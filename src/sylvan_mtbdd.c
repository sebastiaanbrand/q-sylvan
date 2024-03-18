/*
 * Copyright 2011-2016 Formal Methods and Tools, University of Twente
 * Copyright 2016-2017 Tom van Dijk, Johannes Kepler University Linz
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <sylvan_int.h>

#include <inttypes.h>
#include <math.h>
#include <string.h>

#include <sylvan_refs.h>
#include <sylvan_sl.h>
#include <sha2.h>

#include <sylvan_cache.h>

#include <sylvan_mpc.h>

/* Primitives */
int
mtbdd_isleaf(MTBDD bdd)
{
    if (bdd == mtbdd_true || bdd == mtbdd_false) return 1;
    return mtbddnode_isleaf(MTBDD_GETNODE(bdd));
}

// for nodes
uint32_t
mtbdd_getvar(MTBDD node)
{
    return mtbddnode_getvariable(MTBDD_GETNODE(node));
}

MTBDD
mtbdd_getlow(MTBDD mtbdd)
{
    return node_getlow(mtbdd, MTBDD_GETNODE(mtbdd));
}

MTBDD
mtbdd_gethigh(MTBDD mtbdd)
{
    return node_gethigh(mtbdd, MTBDD_GETNODE(mtbdd));
}

// for leaves
uint32_t
mtbdd_gettype(MTBDD leaf)
{
    return mtbddnode_gettype(MTBDD_GETNODE(leaf));
}

uint64_t
mtbdd_getvalue(MTBDD leaf)
{
    return mtbddnode_getvalue(MTBDD_GETNODE(leaf));
}

// for leaf type 0 (integer)
int64_t
mtbdd_getint64(MTBDD leaf)
{
    uint64_t value = mtbdd_getvalue(leaf);
    return *(int64_t*)&value;
}

// for leaf type 1 (double)
double
mtbdd_getdouble(MTBDD leaf)
{
    uint64_t value = mtbdd_getvalue(leaf);
    return *(double*)&value;
}

/**
 * Implementation of garbage collection
 */

/* Recursively mark MDD nodes as 'in use' */
VOID_TASK_IMPL_1(mtbdd_gc_mark_rec, MDD, mtbdd)
{
    if (mtbdd == mtbdd_true) return;
    if (mtbdd == mtbdd_false) return;

    if (llmsset_mark(nodes, MTBDD_STRIPMARK(mtbdd))) {
        mtbddnode_t n = MTBDD_GETNODE(mtbdd);
        if (!mtbddnode_isleaf(n)) {
            SPAWN(mtbdd_gc_mark_rec, mtbddnode_getlow(n));
            CALL(mtbdd_gc_mark_rec, mtbddnode_gethigh(n));
            SYNC(mtbdd_gc_mark_rec);
        }
    }
}

/**
 * External references
 */

refs_table_t mtbdd_refs;
refs_table_t mtbdd_protected;
static int mtbdd_protected_created = 0;

MDD
mtbdd_ref(MDD a)
{
    if (a == mtbdd_true || a == mtbdd_false) return a;
    refs_up(&mtbdd_refs, MTBDD_STRIPMARK(a));
    return a;
}

void
mtbdd_deref(MDD a)
{
    if (a == mtbdd_true || a == mtbdd_false) return;
    refs_down(&mtbdd_refs, MTBDD_STRIPMARK(a));
}

size_t
mtbdd_count_refs()
{
    return refs_count(&mtbdd_refs);
}

void
mtbdd_protect(MTBDD *a)
{
    if (!mtbdd_protected_created) {
        // In C++, sometimes mtbdd_protect is called before Sylvan is initialized. Just create a table.
        protect_create(&mtbdd_protected, 4096);
        mtbdd_protected_created = 1;
    }
    protect_up(&mtbdd_protected, (size_t)a);
}

void
mtbdd_unprotect(MTBDD *a)
{
    if (mtbdd_protected.refs_table != NULL) protect_down(&mtbdd_protected, (size_t)a);
}

size_t
mtbdd_count_protected()
{
    return protect_count(&mtbdd_protected);
}

/* Called during garbage collection */
VOID_TASK_0(mtbdd_gc_mark_external_refs)
{
    // iterate through refs hash table, mark all found
    size_t count=0;
    uint64_t *it = refs_iter(&mtbdd_refs, 0, mtbdd_refs.refs_size);
    while (it != NULL) {
        SPAWN(mtbdd_gc_mark_rec, refs_next(&mtbdd_refs, &it, mtbdd_refs.refs_size));
        count++;
    }
    while (count--) {
        SYNC(mtbdd_gc_mark_rec);
    }
}

VOID_TASK_0(mtbdd_gc_mark_protected)
{
    // iterate through refs hash table, mark all found
    size_t count=0;
    uint64_t *it = protect_iter(&mtbdd_protected, 0, mtbdd_protected.refs_size);
    while (it != NULL) {
        BDD *to_mark = (BDD*)protect_next(&mtbdd_protected, &it, mtbdd_protected.refs_size);
        SPAWN(mtbdd_gc_mark_rec, *to_mark);
        count++;
    }
    while (count--) {
        SYNC(mtbdd_gc_mark_rec);
    }
}

/* Infrastructure for internal markings */
typedef struct mtbdd_refs_task
{
    Task *t;
    void *f;
} *mtbdd_refs_task_t;

typedef struct mtbdd_refs_internal
{
    const MTBDD **pbegin, **pend, **pcur;
    MTBDD *rbegin, *rend, *rcur;
    mtbdd_refs_task_t sbegin, send, scur;
} *mtbdd_refs_internal_t;

DECLARE_THREAD_LOCAL(mtbdd_refs_key, mtbdd_refs_internal_t);

VOID_TASK_2(mtbdd_refs_mark_p_par, const MTBDD**, begin, size_t, count)
{
    if (count < 32) {
        while (count) {
            mtbdd_gc_mark_rec(**(begin++));
            count--;
        }
    } else {
        SPAWN(mtbdd_refs_mark_p_par, begin, count / 2);
        CALL(mtbdd_refs_mark_p_par, begin + (count / 2), count - count / 2);
        SYNC(mtbdd_refs_mark_p_par);
    }
}

VOID_TASK_2(mtbdd_refs_mark_r_par, MTBDD*, begin, size_t, count)
{
    if (count < 32) {
        while (count) {
            mtbdd_gc_mark_rec(*begin++);
            count--;
        }
    } else {
        SPAWN(mtbdd_refs_mark_r_par, begin, count / 2);
        CALL(mtbdd_refs_mark_r_par, begin + (count / 2), count - count / 2);
        SYNC(mtbdd_refs_mark_r_par);
    }
}

VOID_TASK_2(mtbdd_refs_mark_s_par, mtbdd_refs_task_t, begin, size_t, count)
{
    if (count < 32) {
        while (count > 0) {
            Task *t = begin->t;
            if (!TASK_IS_STOLEN(t)) return;
            if (t->f == begin->f && TASK_IS_COMPLETED(t)) {
                mtbdd_gc_mark_rec(*(MTBDD*)TASK_RESULT(t));
            }
            begin += 1;
            count -= 1;
        }
    } else {
        if (!TASK_IS_STOLEN(begin->t)) return;
        SPAWN(mtbdd_refs_mark_s_par, begin, count / 2);
        CALL(mtbdd_refs_mark_s_par, begin + (count / 2), count - count / 2);
        SYNC(mtbdd_refs_mark_s_par);
    }
}

VOID_TASK_0(mtbdd_refs_mark_task)
{
    LOCALIZE_THREAD_LOCAL(mtbdd_refs_key, mtbdd_refs_internal_t);
    SPAWN(mtbdd_refs_mark_p_par, mtbdd_refs_key->pbegin, mtbdd_refs_key->pcur-mtbdd_refs_key->pbegin);
    SPAWN(mtbdd_refs_mark_r_par, mtbdd_refs_key->rbegin, mtbdd_refs_key->rcur-mtbdd_refs_key->rbegin);
    CALL(mtbdd_refs_mark_s_par, mtbdd_refs_key->sbegin, mtbdd_refs_key->scur-mtbdd_refs_key->sbegin);
    SYNC(mtbdd_refs_mark_r_par);
    SYNC(mtbdd_refs_mark_p_par);
}

VOID_TASK_0(mtbdd_refs_mark)
{
    TOGETHER(mtbdd_refs_mark_task);
}

void
mtbdd_refs_init_key(void)
{
    assert(lace_is_worker()); // only use inside Lace workers
    mtbdd_refs_internal_t s = (mtbdd_refs_internal_t)malloc(sizeof(struct mtbdd_refs_internal));
    s->pcur = s->pbegin = (const MTBDD**)malloc(sizeof(MTBDD*) * 1024);
    s->pend = s->pbegin + 1024;
    s->rcur = s->rbegin = (MTBDD*)malloc(sizeof(MTBDD) * 1024);
    s->rend = s->rbegin + 1024;
    s->scur = s->sbegin = (mtbdd_refs_task_t)malloc(sizeof(struct mtbdd_refs_task) * 1024);
    s->send = s->sbegin + 1024;
    SET_THREAD_LOCAL(mtbdd_refs_key, s);
}

VOID_TASK_0(mtbdd_refs_init_task)
{
    mtbdd_refs_init_key();
}

VOID_TASK_0(mtbdd_refs_init)
{
    INIT_THREAD_LOCAL(mtbdd_refs_key);
    TOGETHER(mtbdd_refs_init_task);
    sylvan_gc_add_mark(TASK(mtbdd_refs_mark));
}

void
mtbdd_refs_ptrs_up(mtbdd_refs_internal_t mtbdd_refs_key)
{
    size_t cur = mtbdd_refs_key->pcur - mtbdd_refs_key->pbegin;
    size_t size = mtbdd_refs_key->pend - mtbdd_refs_key->pbegin;
    mtbdd_refs_key->pbegin = (const MTBDD**)realloc(mtbdd_refs_key->pbegin, sizeof(MTBDD*) * size * 2);
    mtbdd_refs_key->pcur = mtbdd_refs_key->pbegin + cur;
    mtbdd_refs_key->pend = mtbdd_refs_key->pbegin + (size * 2);
}

MTBDD __attribute__((noinline))
mtbdd_refs_refs_up(mtbdd_refs_internal_t mtbdd_refs_key, MTBDD res)
{
    long size = mtbdd_refs_key->rend - mtbdd_refs_key->rbegin;
    mtbdd_refs_key->rbegin = (MTBDD*)realloc(mtbdd_refs_key->rbegin, sizeof(MTBDD) * size * 2);
    mtbdd_refs_key->rcur = mtbdd_refs_key->rbegin + size;
    mtbdd_refs_key->rend = mtbdd_refs_key->rbegin + (size * 2);
    return res;
}

void __attribute__((noinline))
mtbdd_refs_tasks_up(mtbdd_refs_internal_t mtbdd_refs_key)
{
    long size = mtbdd_refs_key->send - mtbdd_refs_key->sbegin;
    mtbdd_refs_key->sbegin = (mtbdd_refs_task_t)realloc(mtbdd_refs_key->sbegin, sizeof(struct mtbdd_refs_task) * size * 2);
    mtbdd_refs_key->scur = mtbdd_refs_key->sbegin + size;
    mtbdd_refs_key->send = mtbdd_refs_key->sbegin + (size * 2);
}

void __attribute__((unused))
mtbdd_refs_pushptr(const MTBDD *ptr)
{
    LOCALIZE_THREAD_LOCAL(mtbdd_refs_key, mtbdd_refs_internal_t);
    if (mtbdd_refs_key == 0) {
        mtbdd_refs_init_key();
        mtbdd_refs_pushptr(ptr);
    } else {
        *mtbdd_refs_key->pcur++ = ptr;
        if (mtbdd_refs_key->pcur == mtbdd_refs_key->pend) mtbdd_refs_ptrs_up(mtbdd_refs_key);
    }
}

void __attribute__((unused))
mtbdd_refs_popptr(size_t amount)
{
    LOCALIZE_THREAD_LOCAL(mtbdd_refs_key, mtbdd_refs_internal_t);
    mtbdd_refs_key->pcur -= amount;
}

MTBDD __attribute__((unused))
mtbdd_refs_push(MTBDD mtbdd)
{
    LOCALIZE_THREAD_LOCAL(mtbdd_refs_key, mtbdd_refs_internal_t);
    if (mtbdd_refs_key == 0) {
        mtbdd_refs_init_key();
        return mtbdd_refs_push(mtbdd);
    } else {
        *(mtbdd_refs_key->rcur++) = mtbdd;
        if (mtbdd_refs_key->rcur == mtbdd_refs_key->rend) return mtbdd_refs_refs_up(mtbdd_refs_key, mtbdd);
        else return mtbdd;
    }
}

void __attribute__((unused))
mtbdd_refs_pop(long amount)
{
    LOCALIZE_THREAD_LOCAL(mtbdd_refs_key, mtbdd_refs_internal_t);
    mtbdd_refs_key->rcur -= amount;
}

void
mtbdd_refs_spawn(Task *t)
{
    LOCALIZE_THREAD_LOCAL(mtbdd_refs_key, mtbdd_refs_internal_t);
    mtbdd_refs_key->scur->t = t;
    mtbdd_refs_key->scur->f = t->f;
    mtbdd_refs_key->scur += 1;
    if (mtbdd_refs_key->scur == mtbdd_refs_key->send) mtbdd_refs_tasks_up(mtbdd_refs_key);
}

MTBDD
mtbdd_refs_sync(MTBDD result)
{
    LOCALIZE_THREAD_LOCAL(mtbdd_refs_key, mtbdd_refs_internal_t);
    mtbdd_refs_key->scur -= 1;
    return result;
}

/**
 * Initialize and quit functions
 */

static int mtbdd_initialized = 0;

static void
mtbdd_quit()
{
    refs_free(&mtbdd_refs);
    if (mtbdd_protected_created) {
        protect_free(&mtbdd_protected);
        mtbdd_protected_created = 0;
    }

    mtbdd_initialized = 0;
}

void
sylvan_init_mtbdd()
{
    sylvan_init_mt();

    if (mtbdd_initialized) return;
    mtbdd_initialized = 1;

    sylvan_register_quit(mtbdd_quit);
    sylvan_gc_add_mark(TASK(mtbdd_gc_mark_external_refs));
    sylvan_gc_add_mark(TASK(mtbdd_gc_mark_protected));

    refs_create(&mtbdd_refs, 1024);
    if (!mtbdd_protected_created) {
        protect_create(&mtbdd_protected, 4096);
        mtbdd_protected_created = 1;
    }

    RUN(mtbdd_refs_init);
}

/**
 *  Primitives
 */

/** 
 *  mtbdd_makeleaf()
 * 
 *  Create an endnode with value of given type.
 * 
 *  type:   type of the terminal value, 1 = integer, 2 = ...
 *  value:  value of the terminal node
 * 
 *  return: index to node
 */
MTBDD
mtbdd_makeleaf(uint32_t type, uint64_t value)
{
    struct mtbddnode n;
    mtbddnode_makeleaf(&n, type, value);

    int custom = sylvan_mt_has_custom_hash(type);

    int created;    
    uint64_t index = custom ? llmsset_lookupc(nodes, n.a, n.b, &created) : llmsset_lookup(nodes, n.a, n.b, &created);

    if (index == 0) {
        RUN(sylvan_gc);

        index = custom ? llmsset_lookupc(nodes, n.a, n.b, &created) : llmsset_lookup(nodes, n.a, n.b, &created);
        if (index == 0) {
            fprintf(stderr, "BDD Unique table full, %zu of %zu buckets filled!\n", llmsset_count_marked(nodes), llmsset_get_size(nodes));
            exit(1);
        }
    }

    if (created) sylvan_stats_count(BDD_NODES_CREATED);
    else sylvan_stats_count(BDD_NODES_REUSED);

    return (MTBDD)index;
}

void
__attribute__ ((noinline))
_mtbdd_makenode_gc(MTBDD low, MTBDD high)
{
    mtbdd_refs_push(low);
    mtbdd_refs_push(high);
    RUN(sylvan_gc);
    mtbdd_refs_pop(2);
}

void
__attribute__ ((noinline))
_mtbdd_makenode_exit(void)
{
    fprintf(stderr, "BDD Unique table full, %zu of %zu buckets filled!\n", llmsset_count_marked(nodes), llmsset_get_size(nodes));
    exit(1);
}

/** 
 *  mtbdd_makenode()
 * 
 *  Create a new node with type mtbddnode.
 *
 *  var:    integer, index i of xi
 *  low:    index of low node
 *  high:   index of high node
 * 
 *  return: index to node
 */
MTBDD
_mtbdd_makenode(uint32_t var, MTBDD low, MTBDD high)
{
    // Normalization to keep canonicity
    // low will have no mark

    // TODO: check if selfreference happens (== cyclic graph): 
    // v_low = mtbdd_getvar(low)
    // v_high = mtbdd_getvar(high)
    // v < v_low and v < v_high
    // -> rise assert()? No, always release mode, assert will not work.
    // -> printf();
    //    _mtbdd_makenode_exit();
    //

    MTBDD result = low & mtbdd_complement;
    low ^= result;
    high ^= result;

    struct mtbddnode n;
    mtbddnode_makenode(&n, var, low, high);

    int created;
    uint64_t index = llmsset_lookup(nodes, n.a, n.b, &created);
    if (index == 0) {
        _mtbdd_makenode_gc(low, high);
        index = llmsset_lookup(nodes, n.a, n.b, &created);
        if (index == 0) _mtbdd_makenode_exit();
    }

    if (created) sylvan_stats_count(BDD_NODES_CREATED);
    else sylvan_stats_count(BDD_NODES_REUSED);

    result |= index;
    return result;
}

MTBDD
mtbdd_makemapnode(uint32_t var, MTBDD low, MTBDD high)
{
    struct mtbddnode n;
    uint64_t index;
    int created;

    // in an MTBDDMAP, the low edges eventually lead to 0 and cannot have a low mark
    assert(!MTBDD_HASMARK(low));

    mtbddnode_makemapnode(&n, var, low, high);
    index = llmsset_lookup(nodes, n.a, n.b, &created);
    if (index == 0) {
        mtbdd_refs_push(low);
        mtbdd_refs_push(high);
        RUN(sylvan_gc);
        mtbdd_refs_pop(2);

        index = llmsset_lookup(nodes, n.a, n.b, &created);
        if (index == 0) {
            fprintf(stderr, "BDD Unique table full, %zu of %zu buckets filled!\n", llmsset_count_marked(nodes), llmsset_get_size(nodes));
            exit(1);
        }
    }

    if (created) sylvan_stats_count(BDD_NODES_CREATED);
    else sylvan_stats_count(BDD_NODES_REUSED);

    return index;
}

MTBDD
mtbdd_ithvar(uint32_t var)
{
    return mtbdd_makenode(var, mtbdd_false, mtbdd_true);
}

/* Operations */

/**
 * Calculate greatest common divisor
 * Source: http://lemire.me/blog/archives/2013/12/26/fastest-way-to-compute-the-greatest-common-divisor/
 */
uint32_t
gcd(uint32_t u, uint32_t v)
{
    int shift;
    if (u == 0) return v;
    if (v == 0) return u;
    shift = __builtin_ctz(u | v);
    u >>= __builtin_ctz(u);
    do {
        v >>= __builtin_ctz(v);
        if (u > v) {
            unsigned int t = v;
            v = u;
            u = t;
        }
        v = v - u;
    } while (v != 0);
    return u << shift;
}

/**
 * Create leaves of unsigned/signed integers and doubles
 */

MTBDD
mtbdd_int64(int64_t value)
{
    return mtbdd_makeleaf(0, *(uint64_t*)&value);
}

MTBDD
mtbdd_double(double value)
{
    // normalize all 0.0 to 0.0
    if (value == 0.0) value = 0.0;
    return mtbdd_makeleaf(1, *(uint64_t*)&value);
}

MTBDD
mtbdd_fraction(int64_t nom, uint64_t denom)
{
    if (nom == 0) return mtbdd_makeleaf(2, 1);
    uint32_t c = gcd(nom < 0 ? -nom : nom, denom);
    nom /= c;
    denom /= c;
    if (nom > 2147483647 || nom < -2147483647 || denom > 4294967295) fprintf(stderr, "mtbdd_fraction: fraction overflow\n");
    return mtbdd_makeleaf(2, (nom<<32)|denom);
}

MTBDD
mtbdd_complex_mpc(complex_mpc_t value)
{
    uint32_t terminal_type = MPC_TYPE; // Custom type for complex mpc
    return mtbdd_makeleaf(terminal_type, *(uint64_t*)&value);
}

MTBDD
mtbdd_complex_double(complex_double_t value)
{
    uint32_t terminal_type = 4; // Custom type for complex double
    return mtbdd_makeleaf(terminal_type, *(uint64_t*)&value);
}

/**
 * Create a MTBDD cube representing the conjunction of variables in their positive or negative
 * form depending on whether the cube[idx] equals 0 (negative), 1 (positive) or 2 (any).
 * Use cube[idx]==3 for "s=s'" in interleaved variables (matches with next variable)
 * <variables> is the cube of variables
 */
MTBDD
mtbdd_cube(MTBDD variables, uint8_t *cube, MTBDD terminal)
{
    if (variables == mtbdd_true) return terminal;
    mtbddnode_t n = MTBDD_GETNODE(variables);

    BDD result;
    switch (*cube) {
    case 0:
        result = mtbdd_cube(node_gethigh(variables, n), cube+1, terminal);
        result = mtbdd_makenode(mtbddnode_getvariable(n), result, mtbdd_false);
        return result;
    case 1:
        result = mtbdd_cube(node_gethigh(variables, n), cube+1, terminal);
        result = mtbdd_makenode(mtbddnode_getvariable(n), mtbdd_false, result);
        return result;
    case 2:
        return mtbdd_cube(node_gethigh(variables, n), cube+1, terminal);
    case 3:
    {
        MTBDD variables2 = node_gethigh(variables, n);
        mtbddnode_t n2 = MTBDD_GETNODE(variables2);
        uint32_t var2 = mtbddnode_getvariable(n2);
        result = mtbdd_cube(node_gethigh(variables2, n2), cube+2, terminal);
        BDD low = mtbdd_makenode(var2, result, mtbdd_false);
        mtbdd_refs_push(low);
        BDD high = mtbdd_makenode(var2, mtbdd_false, result);
        mtbdd_refs_pop(1);
        result = mtbdd_makenode(mtbddnode_getvariable(n), low, high);
        return result;
    }
    default:
        return mtbdd_false; // ?
    }
}

/**
 * Same as mtbdd_cube, but also performs "or" with existing MTBDD,
 * effectively adding an item to the set
 */
TASK_IMPL_4(MTBDD, mtbdd_union_cube, MTBDD, mtbdd, MTBDD, vars, uint8_t*, cube, MTBDD, terminal)
{
    /* Terminal cases */
    if (mtbdd == terminal) return terminal;
    if (mtbdd == mtbdd_false) return mtbdd_cube(vars, cube, terminal);
    if (vars == mtbdd_true) return terminal;

    sylvan_gc_test();

    mtbddnode_t nv = MTBDD_GETNODE(vars);
    uint32_t v = mtbddnode_getvariable(nv);

    mtbddnode_t na = MTBDD_GETNODE(mtbdd);
    uint32_t va = mtbddnode_getvariable(na);

    if (va < v) {
        MTBDD low = node_getlow(mtbdd, na);
        MTBDD high = node_gethigh(mtbdd, na);
        mtbdd_refs_spawn(SPAWN(mtbdd_union_cube, high, vars, cube, terminal));
        BDD new_low = mtbdd_union_cube(low, vars, cube, terminal);
        mtbdd_refs_push(new_low);
        BDD new_high = mtbdd_refs_sync(SYNC(mtbdd_union_cube));
        mtbdd_refs_pop(1);
        if (new_low != low || new_high != high) return mtbdd_makenode(va, new_low, new_high);
        else return mtbdd;
    } else if (va == v) {
        MTBDD low = node_getlow(mtbdd, na);
        MTBDD high = node_gethigh(mtbdd, na);
        switch (*cube) {
        case 0:
        {
            MTBDD new_low = mtbdd_union_cube(low, node_gethigh(vars, nv), cube+1, terminal);
            if (new_low != low) return mtbdd_makenode(v, new_low, high);
            else return mtbdd;
        }
        case 1:
        {
            MTBDD new_high = mtbdd_union_cube(high, node_gethigh(vars, nv), cube+1, terminal);
            if (new_high != high) return mtbdd_makenode(v, low, new_high);
            return mtbdd;
        }
        case 2:
        {
            mtbdd_refs_spawn(SPAWN(mtbdd_union_cube, high, node_gethigh(vars, nv), cube+1, terminal));
            MTBDD new_low = mtbdd_union_cube(low, node_gethigh(vars, nv), cube+1, terminal);
            mtbdd_refs_push(new_low);
            MTBDD new_high = mtbdd_refs_sync(SYNC(mtbdd_union_cube));
            mtbdd_refs_pop(1);
            if (new_low != low || new_high != high) return mtbdd_makenode(v, new_low, new_high);
            return mtbdd;
        }
        case 3:
        {
            return mtbdd_false; // currently not implemented
        }
        default:
            return mtbdd_false;
        }
    } else /* va > v */ {
        switch (*cube) {
        case 0:
        {
            MTBDD new_low = mtbdd_union_cube(mtbdd, node_gethigh(vars, nv), cube+1, terminal);
            return mtbdd_makenode(v, new_low, mtbdd_false);
        }
        case 1:
        {
            MTBDD new_high = mtbdd_union_cube(mtbdd, node_gethigh(vars, nv), cube+1, terminal);
            return mtbdd_makenode(v, mtbdd_false, new_high);
        }
        case 2:
        {
            mtbdd_refs_spawn(SPAWN(mtbdd_union_cube, mtbdd, node_gethigh(vars, nv), cube+1, terminal));
            MTBDD new_low = mtbdd_union_cube(mtbdd, node_gethigh(vars, nv), cube+1, terminal);
            mtbdd_refs_push(new_low);
            MTBDD new_high = mtbdd_refs_sync(SYNC(mtbdd_union_cube));
            mtbdd_refs_pop(1);
            return mtbdd_makenode(v, new_low, new_high);
        }
        case 3:
        {
            return mtbdd_false; // currently not implemented
        }
        default:
            return mtbdd_false;
        }
    }
}

/**
 * Apply a binary operation <op> to <a> and <b>.
 */
TASK_IMPL_3(MTBDD, mtbdd_apply, MTBDD, a, MTBDD, b, mtbdd_apply_op, op)
{
    /* Check terminal case: a and b are leaves except for a=0 or b=0 */
    MTBDD result = WRAP(op, &a, &b);
    if (result != mtbdd_invalid) return result;

    /* Maybe perform garbage collection */
    sylvan_gc_test();

    /* Count operation */
    sylvan_stats_count(MTBDD_APPLY);

    /* Check cache */
    if (cache_get3(CACHE_MTBDD_APPLY, a, b, (size_t)op, &result)) {
        sylvan_stats_count(MTBDD_APPLY_CACHED);
        return result;
    }

    /* Get top variable */
    int la = mtbdd_isleaf(a);
    int lb = mtbdd_isleaf(b);
    assert(!la || !lb);
    mtbddnode_t na, nb;
    uint32_t va, vb;
    if (!la) {
        na = MTBDD_GETNODE(a);
        va = mtbddnode_getvariable(na);
    } else {
        na = 0;
        va = 0xffffffff;
    }
    if (!lb) {
        nb = MTBDD_GETNODE(b);
        vb = mtbddnode_getvariable(nb);
    } else {
        nb = 0;
        vb = 0xffffffff;
    }
    uint32_t v = va < vb ? va : vb;

    /* Get cofactors */
    MTBDD alow, ahigh, blow, bhigh;
    if (!la && va == v) {
        alow = node_getlow(a, na);
        ahigh = node_gethigh(a, na);
    } else {
        alow = a;
        ahigh = a;
    }
    if (!lb && vb == v) {
        blow = node_getlow(b, nb);
        bhigh = node_gethigh(b, nb);
    } else {
        blow = b;
        bhigh = b;
    }

    /* Recursive */
    mtbdd_refs_spawn(SPAWN(mtbdd_apply, ahigh, bhigh, op));
    MTBDD low = mtbdd_refs_push(CALL(mtbdd_apply, alow, blow, op));
    MTBDD high = mtbdd_refs_sync(SYNC(mtbdd_apply));
    mtbdd_refs_pop(1);
    result = mtbdd_makenode(v, low, high);

    /* Store in cache */
    if (cache_put3(CACHE_MTBDD_APPLY, a, b, (size_t)op, result)) {
        sylvan_stats_count(MTBDD_APPLY_CACHEDPUT);
    }

    return result;
}

/**
 * Apply a binary operation <op> to <a> and <b> with parameter <p>
 */
TASK_IMPL_5(MTBDD, mtbdd_applyp, MTBDD, a, MTBDD, b, size_t, p, mtbdd_applyp_op, op, uint64_t, opid)
{
    /* Check terminal case */
    MTBDD result = WRAP(op, &a, &b, p);
    if (result != mtbdd_invalid) return result;

    /* Maybe perform garbage collection */
    sylvan_gc_test();

    /* Count operation */
    sylvan_stats_count(MTBDD_APPLY);

    /* Check cache */
    if (cache_get3(opid, a, b, p, &result)) {
        sylvan_stats_count(MTBDD_APPLY_CACHED);
        return result;
    }

    /* Get top variable */
    int la = mtbdd_isleaf(a);
    int lb = mtbdd_isleaf(b);
    mtbddnode_t na, nb;
    uint32_t va, vb;
    if (!la) {
        na = MTBDD_GETNODE(a);
        va = mtbddnode_getvariable(na);
    } else {
        na = 0;
        va = 0xffffffff;
    }
    if (!lb) {
        nb = MTBDD_GETNODE(b);
        vb = mtbddnode_getvariable(nb);
    } else {
        nb = 0;
        vb = 0xffffffff;
    }
    uint32_t v = va < vb ? va : vb;

    /* Get cofactors */
    MTBDD alow, ahigh, blow, bhigh;
    if (!la && va == v) {
        alow = node_getlow(a, na);
        ahigh = node_gethigh(a, na);
    } else {
        alow = a;
        ahigh = a;
    }
    if (!lb && vb == v) {
        blow = node_getlow(b, nb);
        bhigh = node_gethigh(b, nb);
    } else {
        blow = b;
        bhigh = b;
    }

    /* Recursive */
    mtbdd_refs_spawn(SPAWN(mtbdd_applyp, ahigh, bhigh, p, op, opid));
    MTBDD low = mtbdd_refs_push(CALL(mtbdd_applyp, alow, blow, p, op, opid));
    MTBDD high = mtbdd_refs_sync(SYNC(mtbdd_applyp));
    mtbdd_refs_pop(1);
    result = mtbdd_makenode(v, low, high);

    /* Store in cache */
    if (cache_put3(opid, a, b, p, result)) {
        sylvan_stats_count(MTBDD_APPLY_CACHEDPUT);
    }

    return result;
}

/**
 * Apply a unary operation <op> to <dd>.
 */
TASK_IMPL_3(MTBDD, mtbdd_uapply, MTBDD, dd, mtbdd_uapply_op, op, size_t, param)
{
    /* Maybe perform garbage collection */
    sylvan_gc_test();

    /* Count operation */
    sylvan_stats_count(MTBDD_UAPPLY);

    /* Check cache */
    MTBDD result;
    if (cache_get3(CACHE_MTBDD_UAPPLY, dd, (size_t)op, param, &result)) {
        sylvan_stats_count(MTBDD_UAPPLY_CACHED);
        return result;
    }

    /* Check terminal case */
    result = WRAP(op, dd, param);
    if (result != mtbdd_invalid) {
        /* Store in cache */
        if (cache_put3(CACHE_MTBDD_UAPPLY, dd, (size_t)op, param, result)) {
            sylvan_stats_count(MTBDD_UAPPLY_CACHEDPUT);
        }

        return result;
    }

    /* Get cofactors */
    mtbddnode_t ndd = MTBDD_GETNODE(dd);
    MTBDD ddlow = node_getlow(dd, ndd);
    MTBDD ddhigh = node_gethigh(dd, ndd);

    /* Recursive */
    mtbdd_refs_spawn(SPAWN(mtbdd_uapply, ddhigh, op, param));
    MTBDD low = mtbdd_refs_push(CALL(mtbdd_uapply, ddlow, op, param));
    MTBDD high = mtbdd_refs_sync(SYNC(mtbdd_uapply));
    mtbdd_refs_pop(1);
    result = mtbdd_makenode(mtbddnode_getvariable(ndd), low, high);

    /* Store in cache */
    if (cache_put3(CACHE_MTBDD_UAPPLY, dd, (size_t)op, param, result)) {
        sylvan_stats_count(MTBDD_UAPPLY_CACHEDPUT);
    }

    return result;
}

TASK_2(MTBDD, mtbdd_uop_times_uint, MTBDD, a, size_t, k)
{
    if (a == mtbdd_false) return mtbdd_false;
    if (a == mtbdd_true) return mtbdd_true;

    // a != constant
    mtbddnode_t na = MTBDD_GETNODE(a);

    if (mtbddnode_isleaf(na)) {
        if (mtbddnode_gettype(na) == 0) {
            int64_t v = mtbdd_getint64(a);
            return mtbdd_int64(v*k);
        } else if (mtbddnode_gettype(na) == 1) {
            double d = mtbdd_getdouble(a);
            return mtbdd_double(d*k);
        } else if (mtbddnode_gettype(na) == 2) {
            uint64_t v = mtbddnode_getvalue(na);
            int64_t n = (int32_t)(v>>32);
            uint32_t d = v;
            uint32_t c = gcd(d, (uint32_t)k);
            return mtbdd_fraction(n*(k/c), d/c);
        } else {
            assert(0); // failure
        }
    }

    return mtbdd_invalid;
}

TASK_2(MTBDD, mtbdd_uop_pow_uint, MTBDD, a, size_t, k)
{
    if (a == mtbdd_false) return mtbdd_false;
    if (a == mtbdd_true) return mtbdd_true;

    // a != constant
    mtbddnode_t na = MTBDD_GETNODE(a);

    if (mtbddnode_isleaf(na)) {
        if (mtbddnode_gettype(na) == 0) {
            int64_t v = mtbdd_getint64(a);
            return mtbdd_int64(pow(v, k));
        } else if (mtbddnode_gettype(na) == 1) {
            double d = mtbdd_getdouble(a);
            return mtbdd_double(pow(d, k));
        } else if (mtbddnode_gettype(na) == 2) {
            uint64_t v = mtbddnode_getvalue(na);
            return mtbdd_fraction(pow((int32_t)(v>>32), k), (uint32_t)v);
        } else {
            assert(0); // failure
        }
    }

    return mtbdd_invalid;
}

TASK_IMPL_3(MTBDD, mtbdd_abstract_op_plus, MTBDD, a, MTBDD, b, int, k)
{
    if (k==0) {
        return mtbdd_apply(a, b, TASK(mtbdd_op_plus));
    } else {
        uint64_t factor = 1ULL<<k; // skip 1,2,3,4: times 2,4,8,16
        return mtbdd_uapply(a, TASK(mtbdd_uop_times_uint), factor);
    }
}

TASK_IMPL_3(MTBDD, mtbdd_abstract_op_times, MTBDD, a, MTBDD, b, int, k)
{
    if (k==0) {
        return mtbdd_apply(a, b, TASK(mtbdd_op_times));
    } else {
        uint64_t squares = 1ULL<<k; // square k times, ie res^(2^k): 2,4,8,16
        return mtbdd_uapply(a, TASK(mtbdd_uop_pow_uint), squares);
    }
}

TASK_IMPL_3(MTBDD, mtbdd_abstract_op_min, MTBDD, a, MTBDD, b, int, k)
{
    return k == 0 ? mtbdd_apply(a, b, TASK(mtbdd_op_min)) : a;
}

TASK_IMPL_3(MTBDD, mtbdd_abstract_op_max, MTBDD, a, MTBDD, b, int, k)
{
    return k == 0 ? mtbdd_apply(a, b, TASK(mtbdd_op_max)) : a;
}

/**
 * Abstract the variables in <v> from <a> using the operation <op>
 */
TASK_IMPL_3(MTBDD, mtbdd_abstract, MTBDD, a, MTBDD, v, mtbdd_abstract_op, op)
{
    /* Check terminal case */
    if (a == mtbdd_false) return mtbdd_false;
    if (a == mtbdd_true) return mtbdd_true;
    if (v == mtbdd_true) return a;

    /* Maybe perform garbage collection */
    sylvan_gc_test();

    /* Count operation */
    sylvan_stats_count(MTBDD_ABSTRACT);

    /* a != constant, v != constant */
    mtbddnode_t na = MTBDD_GETNODE(a);

    if (mtbddnode_isleaf(na)) {
        /* Count number of variables */
        uint64_t k = 0;
        while (v != mtbdd_true) {
            k++;
            v = node_gethigh(v, MTBDD_GETNODE(v));
        }

        /* Check cache */
        MTBDD result;
        if (cache_get3(CACHE_MTBDD_ABSTRACT, a, v | (k << 40), (size_t)op, &result)) {
            sylvan_stats_count(MTBDD_ABSTRACT_CACHED);
            return result;
        }

        /* Compute result */
        result = WRAP(op, a, a, k);

        /* Store in cache */
        if (cache_put3(CACHE_MTBDD_ABSTRACT, a, v | (k << 40), (size_t)op, result)) {
            sylvan_stats_count(MTBDD_ABSTRACT_CACHEDPUT);
        }

        return result;
    }

    /* Possibly skip k variables */
    mtbddnode_t nv = MTBDD_GETNODE(v);
    uint32_t var_a = mtbddnode_getvariable(na);
    uint32_t var_v = mtbddnode_getvariable(nv);
    uint64_t k = 0;
    while (var_v < var_a) {
        k++;
        v = node_gethigh(v, nv);
        if (v == mtbdd_true) break;
        nv = MTBDD_GETNODE(v);
        var_v = mtbddnode_getvariable(nv);
    }

    /* Check cache */
    MTBDD result;
    if (cache_get3(CACHE_MTBDD_ABSTRACT, a, v | (k << 40), (size_t)op, &result)) {
        sylvan_stats_count(MTBDD_ABSTRACT_CACHED);
        return result;
    }

    /* Recursive */
    if (v == mtbdd_true) {
        result = a;
    } else if (var_a < var_v) {
        mtbdd_refs_spawn(SPAWN(mtbdd_abstract, node_gethigh(a, na), v, op));
        MTBDD low = mtbdd_refs_push(CALL(mtbdd_abstract, node_getlow(a, na), v, op));
        MTBDD high = mtbdd_refs_sync(SYNC(mtbdd_abstract));
        mtbdd_refs_pop(1);
        result = mtbdd_makenode(var_a, low, high);
    } else /* var_a == var_v */ {
        mtbdd_refs_spawn(SPAWN(mtbdd_abstract, node_gethigh(a, na), node_gethigh(v, nv), op));
        MTBDD low = mtbdd_refs_push(CALL(mtbdd_abstract, node_getlow(a, na), node_gethigh(v, nv), op));
        MTBDD high = mtbdd_refs_push(mtbdd_refs_sync(SYNC(mtbdd_abstract)));
        result = WRAP(op, low, high, 0);
        mtbdd_refs_pop(2);
    }

    if (k) {
        mtbdd_refs_push(result);
        result = WRAP(op, result, result, k);
        mtbdd_refs_pop(1);
    }

    /* Store in cache */
    if (cache_put3(CACHE_MTBDD_ABSTRACT, a, v | (k << 40), (size_t)op, result)) {
        sylvan_stats_count(MTBDD_ABSTRACT_CACHEDPUT);
    }

    return result;
}

/**
 * Binary operation Plus (for MTBDDs of same type)
 * Only for MTBDDs where either all leaves are Boolean, or Integer, or Double.
 * For Integer/Double MTBDDs, mtbdd_false is interpreted as "0" or "0.0".
 */
TASK_IMPL_2(MTBDD, mtbdd_op_plus, MTBDD*, pa, MTBDD*, pb)
{
    MTBDD a = *pa, b = *pb;

    if (a == mtbdd_false) return b;
    if (b == mtbdd_false) return a;

    // Handle Boolean MTBDDs: interpret as Or
    if (a == mtbdd_true) return mtbdd_true;
    if (b == mtbdd_true) return mtbdd_true;

    mtbddnode_t na = MTBDD_GETNODE(a);
    mtbddnode_t nb = MTBDD_GETNODE(b);

    if (mtbddnode_isleaf(na) && mtbddnode_isleaf(nb)) {

        uint64_t val_a = mtbddnode_getvalue(na);
        uint64_t val_b = mtbddnode_getvalue(nb);
        
        if (mtbddnode_gettype(na) == 0 && mtbddnode_gettype(nb) == 0) {
            // both integer
            return mtbdd_int64(*(int64_t*)(&val_a) + *(int64_t*)(&val_b));
        
        } else if (mtbddnode_gettype(na) == 1 && mtbddnode_gettype(nb) == 1) {
            // both double
            return mtbdd_double(*(double*)(&val_a) + *(double*)(&val_b));
        
        } else if (mtbddnode_gettype(na) == 2 && mtbddnode_gettype(nb) == 2) {
        
            // both fraction
            int64_t nom_a = (int32_t)(val_a>>32);
            int64_t nom_b = (int32_t)(val_b>>32);
            uint64_t denom_a = val_a&0xffffffff;
            uint64_t denom_b = val_b&0xffffffff;
        
            // common cases
            if (nom_a == 0) return b;
            if (nom_b == 0) return a;
        
            // equalize denominators
            uint32_t c = gcd(denom_a, denom_b);
            nom_a *= denom_b/c;
            nom_b *= denom_a/c;
            denom_a *= denom_b/c;
        
            // add
            return mtbdd_fraction(nom_a + nom_b, denom_a);
        
        } else {

            assert(mtbdd_gettype(a) == MPC_TYPE && mtbdd_gettype(b) == MPC_TYPE);

            mpc_ptr ma = (mpc_ptr)val_a;
            mpc_ptr mb = (mpc_ptr)val_b;

            mpc_t x;
            mpc_init2(x, MPC_PRECISION);
            mpc_add(x, ma, mb, MPC_ROUNDING);

            MTBDD result = mtbdd_makeleaf(MPC_TYPE, (size_t)x);
            mpc_clear(x);
            return result;        
        }
    }

    if (a < b) {
        *pa = b;
        *pb = a;
    }

    return mtbdd_invalid;
}

/**
 * Binary operation Minus (for MTBDDs of same type)
 * Only for MTBDDs where either all leaves are Boolean, or Integer, or Double.
 * For Integer/Double MTBDDs, mtbdd_false is interpreted as "0" or "0.0".
 */
TASK_IMPL_2(MTBDD, mtbdd_op_minus, MTBDD*, pa, MTBDD*, pb)
{
    MTBDD a = *pa, b = *pb;
    if (a == mtbdd_false) return mtbdd_negate(b);
    if (b == mtbdd_false) return a;

    mtbddnode_t na = MTBDD_GETNODE(a);
    mtbddnode_t nb = MTBDD_GETNODE(b);

    if (mtbddnode_isleaf(na) && mtbddnode_isleaf(nb)) {
        uint64_t val_a = mtbddnode_getvalue(na);
        uint64_t val_b = mtbddnode_getvalue(nb);
        if (mtbddnode_gettype(na) == 0 && mtbddnode_gettype(nb) == 0) {
            // both integer
            return mtbdd_int64(*(int64_t*)(&val_a) - *(int64_t*)(&val_b));
        } else if (mtbddnode_gettype(na) == 1 && mtbddnode_gettype(nb) == 1) {
            // both double
            return mtbdd_double(*(double*)(&val_a) - *(double*)(&val_b));
        } else if (mtbddnode_gettype(na) == 2 && mtbddnode_gettype(nb) == 2) {
            // both fraction
            int64_t nom_a = (int32_t)(val_a>>32);
            int64_t nom_b = (int32_t)(val_b>>32);
            uint64_t denom_a = val_a&0xffffffff;
            uint64_t denom_b = val_b&0xffffffff;
            // common cases
            if (nom_b == 0) return a;
            // equalize denominators
            uint32_t c = gcd(denom_a, denom_b);
            nom_a *= denom_b/c;
            nom_b *= denom_a/c;
            denom_a *= denom_b/c;
            // subtract
            return mtbdd_fraction(nom_a - nom_b, denom_a);
        } else {
            assert(0); // failure
        }
    }

    return mtbdd_invalid;
}

/**
 * Binary operation Times (for MTBDDs of same type)
 * Only for MTBDDs where either all leaves are Boolean, or Integer, or Double.
 * For Integer/Double MTBDD, if either operand is mtbdd_false (not defined),
 * then the result is mtbdd_false (i.e. not defined).
 */
TASK_IMPL_2(MTBDD, mtbdd_op_times, MTBDD*, pa, MTBDD*, pb)
{
    MTBDD a = *pa, b = *pb;
    if (a == mtbdd_false || b == mtbdd_false) return mtbdd_false;

    // Handle Boolean MTBDDs: interpret as And
    if (a == mtbdd_true) return b;
    if (b == mtbdd_true) return a;

    mtbddnode_t na = MTBDD_GETNODE(a);
    mtbddnode_t nb = MTBDD_GETNODE(b);

    if (mtbddnode_isleaf(na) && mtbddnode_isleaf(nb)) {
        uint64_t val_a = mtbddnode_getvalue(na);
        uint64_t val_b = mtbddnode_getvalue(nb);
        if (mtbddnode_gettype(na) == 0 && mtbddnode_gettype(nb) == 0) {
            // both integer
            int64_t i_a = *(int64_t*)(&val_a);
            int64_t i_b = *(int64_t*)(&val_b);
            if (i_a == 0) return a;
            if (i_b == 0) return b;
            if (i_a == 1) return b;
            if (i_b == 1) return a;
            return mtbdd_int64(i_a * i_b);
        } else if (mtbddnode_gettype(na) == 1 && mtbddnode_gettype(nb) == 1) {
            // both double
            double d_a = *(double*)(&val_a);
            double d_b = *(double*)(&val_b);
            if (d_a == 0.0) return a;
            if (d_a == 1.0) return b;
            if (d_b == 0.0) return b;
            if (d_b == 1.0) return a;
            return mtbdd_double(d_a * d_b);
        } else if (mtbddnode_gettype(na) == 2 && mtbddnode_gettype(nb) == 2) {
            // both fraction
            int64_t nom_a = (int32_t)(val_a>>32);
            int64_t nom_b = (int32_t)(val_b>>32);
            uint64_t denom_a = val_a&0xffffffff;
            uint64_t denom_b = val_b&0xffffffff;
            if (nom_a == 0) return a;
            if (nom_b == 0) return b;
            // multiply!
            uint32_t c = gcd(nom_b < 0 ? -nom_b : nom_b, denom_a);
            uint32_t d = gcd(nom_a < 0 ? -nom_a : nom_a, denom_b);
            nom_a /= d;
            denom_a /= c;
            nom_a *= (nom_b/c);
            denom_a *= (denom_b/d);
            return mtbdd_fraction(nom_a, denom_a);
        } else {
            assert(0); // failure
        }
    }

    if (a < b) {
        *pa = b;
        *pb = a;
    }

    return mtbdd_invalid;
}

/**
 * Binary operation Minimum (for MTBDDs of same type)
 * Only for MTBDDs where either all leaves are Boolean, or Integer, or Double.
 * For Integer/Double MTBDD, if either operand is mtbdd_false (not defined),
 * then the result is the other operand.
 */
TASK_IMPL_2(MTBDD, mtbdd_op_min, MTBDD*, pa, MTBDD*, pb)
{
    printf("1\n");

    MTBDD a = *pa, b = *pb;
    if (a == mtbdd_true) return b;
    if (b == mtbdd_true) return a;
    if (a == b) return a;

    // Special case where "false" indicates a partial function
    if (a == mtbdd_false && b != mtbdd_false && mtbddnode_isleaf(MTBDD_GETNODE(b))) return b;
    if (b == mtbdd_false && a != mtbdd_false && mtbddnode_isleaf(MTBDD_GETNODE(a))) return a;

    mtbddnode_t na = MTBDD_GETNODE(a);
    mtbddnode_t nb = MTBDD_GETNODE(b);

    printf("2\n");

    if (mtbddnode_isleaf(na) && mtbddnode_isleaf(nb)) {
        uint64_t val_a = mtbddnode_getvalue(na);
        uint64_t val_b = mtbddnode_getvalue(nb);
        if (mtbddnode_gettype(na) == 0 && mtbddnode_gettype(nb) == 0) {
            // both integer
            int64_t va = *(int64_t*)(&val_a);
            int64_t vb = *(int64_t*)(&val_b);
            return va < vb ? a : b;
        } else if (mtbddnode_gettype(na) == 1 && mtbddnode_gettype(nb) == 1) {
            // both double
            double va = *(double*)&val_a;
            double vb = *(double*)&val_b;

            printf("3  %lf  %lf \n", va, vb);

            return va < vb ? a : b;
        } else if (mtbddnode_gettype(na) == 2 && mtbddnode_gettype(nb) == 2) {
            // both fraction
            int64_t nom_a = (int32_t)(val_a>>32);
            int64_t nom_b = (int32_t)(val_b>>32);
            uint64_t denom_a = val_a&0xffffffff;
            uint64_t denom_b = val_b&0xffffffff;
            // equalize denominators
            uint32_t c = gcd(denom_a, denom_b);
            nom_a *= denom_b/c;
            nom_b *= denom_a/c;
            // compute lowest
            return nom_a < nom_b ? a : b;
        } else {
            assert(0); // failure
        }
    }

    if (a < b) {
        *pa = b;
        *pb = a;
    }

    return mtbdd_invalid;
}

/**
 * Binary operation Maximum (for MTBDDs of same type)
 * Only for MTBDDs where either all leaves are Boolean, or Integer, or Double.
 * For Integer/Double MTBDD, if either operand is mtbdd_false (not defined),
 * then the result is the other operand.
 */
TASK_IMPL_2(MTBDD, mtbdd_op_max, MTBDD*, pa, MTBDD*, pb)
{
    MTBDD a = *pa, b = *pb;
    if (a == mtbdd_true) return a;
    if (b == mtbdd_true) return b;
    if (a == mtbdd_false) return b;
    if (b == mtbdd_false) return a;
    if (a == b) return a;

    mtbddnode_t na = MTBDD_GETNODE(a);
    mtbddnode_t nb = MTBDD_GETNODE(b);

    if (mtbddnode_isleaf(na) && mtbddnode_isleaf(nb)) {
        uint64_t val_a = mtbddnode_getvalue(na);
        uint64_t val_b = mtbddnode_getvalue(nb);
        if (mtbddnode_gettype(na) == 0 && mtbddnode_gettype(nb) == 0) {
            // both integer
            int64_t va = *(int64_t*)(&val_a);
            int64_t vb = *(int64_t*)(&val_b);
            return va > vb ? a : b;
        } else if (mtbddnode_gettype(na) == 1 && mtbddnode_gettype(nb) == 1) {
            // both double
            double vval_a = *(double*)&val_a;
            double vval_b = *(double*)&val_b;
            return vval_a > vval_b ? a : b;
        } else if (mtbddnode_gettype(na) == 2 && mtbddnode_gettype(nb) == 2) {
            // both fraction
            int64_t nom_a = (int32_t)(val_a>>32);
            int64_t nom_b = (int32_t)(val_b>>32);
            uint64_t denom_a = val_a&0xffffffff;
            uint64_t denom_b = val_b&0xffffffff;
            // equalize denominators
            uint32_t c = gcd(denom_a, denom_b);
            nom_a *= denom_b/c;
            nom_b *= denom_a/c;
            // compute highest
            return nom_a > nom_b ? a : b;
        } else {
            assert(0); // failure
        }
    }

    if (a < b) {
        *pa = b;
        *pb = a;
    }

    return mtbdd_invalid;
}

TASK_IMPL_2(MTBDD, mtbdd_op_cmpl, MTBDD, a, size_t, k)
{
    // if a is false, then it is a partial function. Keep partial!
    if (a == mtbdd_false) return mtbdd_false;

    // a != constant
    mtbddnode_t na = MTBDD_GETNODE(a);

    if (mtbddnode_isleaf(na)) {
        if (mtbddnode_gettype(na) == 0) {
            int64_t v = mtbdd_getint64(a);
            if (v == 0) return mtbdd_int64(1);
            else return mtbdd_int64(0);
        } else if (mtbddnode_gettype(na) == 1) {
            double d = mtbdd_getdouble(a);
            if (d == 0.0) return mtbdd_double(1.0);
            else return mtbdd_double(0.0);
        } else if (mtbddnode_gettype(na) == 2) {
            uint64_t v = mtbddnode_getvalue(na);
            if (v == 1) return mtbdd_fraction(1, 1);
            else return mtbdd_fraction(0, 1);
        } else {
            assert(0); // failure
        }
    }

    return mtbdd_invalid;
    (void)k; // unused variable
}

TASK_IMPL_2(MTBDD, mtbdd_op_negate, MTBDD, a, size_t, k)
{
    // if a is false, then it is a partial function. Keep partial!
    if (a == mtbdd_false) return mtbdd_false;

    // a != constant
    mtbddnode_t na = MTBDD_GETNODE(a);

    if (mtbddnode_isleaf(na)) {
        if (mtbddnode_gettype(na) == 0) {
            int64_t v = mtbdd_getint64(a);
            return mtbdd_int64(-v);
        } else if (mtbddnode_gettype(na) == 1) {
            double d = mtbdd_getdouble(a);
            return mtbdd_double(-d);
        } else if (mtbddnode_gettype(na) == 2) {
            uint64_t v = mtbddnode_getvalue(na);
            return mtbdd_fraction(-(int32_t)(v>>32), (uint32_t)v);
        } else {
            assert(0); // failure
        }
    }

    return mtbdd_invalid;
    (void)k; // unused variable
}

/**
 * Compute IF <f> THEN <g> ELSE <h>.
 * <f> must be a Boolean MTBDD (or standard BDD).
 */
TASK_IMPL_3(MTBDD, mtbdd_ite, MTBDD, f, MTBDD, g, MTBDD, h)
{
    /* Terminal cases */
    if (f == mtbdd_true) return g;
    if (f == mtbdd_false) return h;
    if (g == h) return g;
    if (g == mtbdd_true && h == mtbdd_false) return f;
    if (h == mtbdd_true && g == mtbdd_false) return MTBDD_TOGGLEMARK(f);

    // If all MTBDD's are Boolean, then there could be further optimizations (see sylvan_bdd.c)

    /* Maybe perform garbage collection */
    sylvan_gc_test();

    /* Count operation */
    sylvan_stats_count(MTBDD_ITE);

    /* Check cache */
    MTBDD result;
    if (cache_get3(CACHE_MTBDD_ITE, f, g, h, &result)) {
        sylvan_stats_count(MTBDD_ITE_CACHED);
        return result;
    }

    /* Get top variable */
    int lg = mtbdd_isleaf(g);
    int lh = mtbdd_isleaf(h);
    mtbddnode_t nf = MTBDD_GETNODE(f);
    mtbddnode_t ng = lg ? 0 : MTBDD_GETNODE(g);
    mtbddnode_t nh = lh ? 0 : MTBDD_GETNODE(h);
    uint32_t vf = mtbddnode_getvariable(nf);
    uint32_t vg = lg ? 0 : mtbddnode_getvariable(ng);
    uint32_t vh = lh ? 0 : mtbddnode_getvariable(nh);
    uint32_t v = vf;
    if (!lg && vg < v) v = vg;
    if (!lh && vh < v) v = vh;

    /* Get cofactors */
    MTBDD flow, fhigh, glow, ghigh, hlow, hhigh;
    flow = (vf == v) ? node_getlow(f, nf) : f;
    fhigh = (vf == v) ? node_gethigh(f, nf) : f;
    glow = (!lg && vg == v) ? node_getlow(g, ng) : g;
    ghigh = (!lg && vg == v) ? node_gethigh(g, ng) : g;
    hlow = (!lh && vh == v) ? node_getlow(h, nh) : h;
    hhigh = (!lh && vh == v) ? node_gethigh(h, nh) : h;

    /* Recursive calls */
    mtbdd_refs_spawn(SPAWN(mtbdd_ite, fhigh, ghigh, hhigh));
    MTBDD low = mtbdd_refs_push(CALL(mtbdd_ite, flow, glow, hlow));
    MTBDD high = mtbdd_refs_sync(SYNC(mtbdd_ite));
    mtbdd_refs_pop(1);
    result = mtbdd_makenode(v, low, high);

    /* Store in cache */
    if (cache_put3(CACHE_MTBDD_ITE, f, g, h, result)) {
        sylvan_stats_count(MTBDD_ITE_CACHEDPUT);
    }

    return result;
}

/**
 * Monad that converts double/fraction to a Boolean MTBDD, translate terminals >= value to 1 and to 0 otherwise;
 */
TASK_IMPL_2(MTBDD, mtbdd_op_threshold_double, MTBDD, a, size_t, svalue)
{
    /* We only expect "double" terminals, or false */
    if (a == mtbdd_false) return mtbdd_false;
    if (a == mtbdd_true) return mtbdd_invalid;

    // a != constant
    mtbddnode_t na = MTBDD_GETNODE(a);

    if (mtbddnode_isleaf(na)) {
        double value = *(double*)&svalue;
        if (mtbddnode_gettype(na) == 1) {
            return mtbdd_getdouble(a) >= value ? mtbdd_true : mtbdd_false;
        } else if (mtbddnode_gettype(na) == 2) {
            double d = (double)mtbdd_getnumer(a);
            d /= mtbdd_getdenom(a);
            return d >= value ? mtbdd_true : mtbdd_false;
        } else {
            assert(0); // failure
        }
    }

    return mtbdd_invalid;
}

/**
 * Monad that converts double/fraction to a Boolean BDD, translate terminals > value to 1 and to 0 otherwise;
 */
TASK_IMPL_2(MTBDD, mtbdd_op_strict_threshold_double, MTBDD, a, size_t, svalue)
{
    /* We only expect "double" terminals, or false */
    if (a == mtbdd_false) return mtbdd_false;
    if (a == mtbdd_true) return mtbdd_invalid;

    // a != constant
    mtbddnode_t na = MTBDD_GETNODE(a);

    if (mtbddnode_isleaf(na)) {
        double value = *(double*)&svalue;
        if (mtbddnode_gettype(na) == 1) {
            return mtbdd_getdouble(a) > value ? mtbdd_true : mtbdd_false;
        } else if (mtbddnode_gettype(na) == 2) {
            double d = (double)mtbdd_getnumer(a);
            d /= mtbdd_getdenom(a);
            return d > value ? mtbdd_true : mtbdd_false;
        } else {
            assert(0); // failure
        }
    }

    return mtbdd_invalid;
}

TASK_IMPL_2(MTBDD, mtbdd_threshold_double, MTBDD, dd, double, d)
{
    return mtbdd_uapply(dd, TASK(mtbdd_op_threshold_double), *(size_t*)&d);
}

TASK_IMPL_2(MTBDD, mtbdd_strict_threshold_double, MTBDD, dd, double, d)
{
    return mtbdd_uapply(dd, TASK(mtbdd_op_strict_threshold_double), *(size_t*)&d);
}

/**
 * Compare two Double MTBDDs, returns Boolean True if they are equal within some value epsilon
 */
TASK_4(MTBDD, mtbdd_equal_norm_d2, MTBDD, a, MTBDD, b, size_t, svalue, int*, shortcircuit)
{
    /* Check short circuit */
    if (*shortcircuit) return mtbdd_false;

    /* Check terminal case */
    if (a == b) return mtbdd_true;
    if (a == mtbdd_false) return mtbdd_false;
    if (b == mtbdd_false) return mtbdd_false;

    mtbddnode_t na = MTBDD_GETNODE(a);
    mtbddnode_t nb = MTBDD_GETNODE(b);
    int la = mtbddnode_isleaf(na);
    int lb = mtbddnode_isleaf(nb);

    if (la && lb) {
        // assume Double MTBDD
        double va = mtbdd_getdouble(a);
        double vb = mtbdd_getdouble(b);
        va -= vb;
        if (va < 0) va = -va;
        return (va < *(double*)&svalue) ? mtbdd_true : mtbdd_false;
    }

    if (b < a) {
        MTBDD t = a;
        a = b;
        b = t;
    }

    /* Maybe perform garbage collection */
    sylvan_gc_test();

    /* Count operation */
    sylvan_stats_count(MTBDD_EQUAL_NORM);

    /* Check cache */
    MTBDD result;
    if (cache_get3(CACHE_MTBDD_EQUAL_NORM, a, b, svalue, &result)) {
        sylvan_stats_count(MTBDD_EQUAL_NORM_CACHED);
        return result;
    }

    /* Get top variable */
    uint32_t va = la ? 0xffffffff : mtbddnode_getvariable(na);
    uint32_t vb = lb ? 0xffffffff : mtbddnode_getvariable(nb);
    uint32_t var = va < vb ? va : vb;

    /* Get cofactors */
    MTBDD alow, ahigh, blow, bhigh;
    alow  = va == var ? node_getlow(a, na)  : a;
    ahigh = va == var ? node_gethigh(a, na) : a;
    blow  = vb == var ? node_getlow(b, nb)  : b;
    bhigh = vb == var ? node_gethigh(b, nb) : b;

    SPAWN(mtbdd_equal_norm_d2, ahigh, bhigh, svalue, shortcircuit);
    result = CALL(mtbdd_equal_norm_d2, alow, blow, svalue, shortcircuit);
    if (result == mtbdd_false) *shortcircuit = 1;
    if (result != SYNC(mtbdd_equal_norm_d2)) result = mtbdd_false;
    if (result == mtbdd_false) *shortcircuit = 1;

    /* Store in cache */
    if (cache_put3(CACHE_MTBDD_EQUAL_NORM, a, b, svalue, result)) {
        sylvan_stats_count(MTBDD_EQUAL_NORM_CACHEDPUT);
    }

    return result;
}

TASK_IMPL_3(MTBDD, mtbdd_equal_norm_d, MTBDD, a, MTBDD, b, double, d)
{
    /* the implementation checks shortcircuit in every task and if the two
       MTBDDs are not equal module epsilon, then the computation tree quickly aborts */
    int shortcircuit = 0;
    return CALL(mtbdd_equal_norm_d2, a, b, *(size_t*)&d, &shortcircuit);
}

/**
 * Compare two Double MTBDDs, returns Boolean True if they are equal within some value epsilon
 * This version computes the relative difference vs the value in a.
 */
TASK_4(MTBDD, mtbdd_equal_norm_rel_d2, MTBDD, a, MTBDD, b, size_t, svalue, int*, shortcircuit)
{
    /* Check short circuit */
    if (*shortcircuit) return mtbdd_false;

    /* Check terminal case */
    if (a == b) return mtbdd_true;
    if (a == mtbdd_false) return mtbdd_false;
    if (b == mtbdd_false) return mtbdd_false;

    mtbddnode_t na = MTBDD_GETNODE(a);
    mtbddnode_t nb = MTBDD_GETNODE(b);
    int la = mtbddnode_isleaf(na);
    int lb = mtbddnode_isleaf(nb);

    if (la && lb) {
        // assume Double MTBDD
        double va = mtbdd_getdouble(a);
        double vb = mtbdd_getdouble(b);
        if (va == 0) return mtbdd_false;
        va = (va - vb) / va;
        if (va < 0) va = -va;
        return (va < *(double*)&svalue) ? mtbdd_true : mtbdd_false;
    }

    /* Maybe perform garbage collection */
    sylvan_gc_test();

    /* Count operation */
    sylvan_stats_count(MTBDD_EQUAL_NORM_REL);

    /* Check cache */
    MTBDD result;
    if (cache_get3(CACHE_MTBDD_EQUAL_NORM_REL, a, b, svalue, &result)) {
        sylvan_stats_count(MTBDD_EQUAL_NORM_REL_CACHED);
        return result;
    }

    /* Get top variable */
    uint32_t va = la ? 0xffffffff : mtbddnode_getvariable(na);
    uint32_t vb = lb ? 0xffffffff : mtbddnode_getvariable(nb);
    uint32_t var = va < vb ? va : vb;

    /* Get cofactors */
    MTBDD alow, ahigh, blow, bhigh;
    alow  = va == var ? node_getlow(a, na)  : a;
    ahigh = va == var ? node_gethigh(a, na) : a;
    blow  = vb == var ? node_getlow(b, nb)  : b;
    bhigh = vb == var ? node_gethigh(b, nb) : b;

    SPAWN(mtbdd_equal_norm_rel_d2, ahigh, bhigh, svalue, shortcircuit);
    result = CALL(mtbdd_equal_norm_rel_d2, alow, blow, svalue, shortcircuit);
    if (result == mtbdd_false) *shortcircuit = 1;
    if (result != SYNC(mtbdd_equal_norm_rel_d2)) result = mtbdd_false;
    if (result == mtbdd_false) *shortcircuit = 1;

    /* Store in cache */
    if (cache_put3(CACHE_MTBDD_EQUAL_NORM_REL, a, b, svalue, result)) {
        sylvan_stats_count(MTBDD_EQUAL_NORM_REL_CACHEDPUT);
    }

    return result;
}

TASK_IMPL_3(MTBDD, mtbdd_equal_norm_rel_d, MTBDD, a, MTBDD, b, double, d)
{
    /* the implementation checks shortcircuit in every task and if the two
       MTBDDs are not equal module epsilon, then the computation tree quickly aborts */
    int shortcircuit = 0;
    return CALL(mtbdd_equal_norm_rel_d2, a, b, *(size_t*)&d, &shortcircuit);
}

/**
 * For two MTBDDs a, b, return mtbdd_true if all common assignments a(s) <= b(s), mtbdd_false otherwise.
 * For domains not in a / b, assume True.
 */
TASK_3(MTBDD, mtbdd_leq_rec, MTBDD, a, MTBDD, b, int*, shortcircuit)
{
    /* Check short circuit */
    if (*shortcircuit) return mtbdd_false;

    /* Check terminal case */
    if (a == b) return mtbdd_true;

    /* For partial functions, just return true */
    if (a == mtbdd_false) return mtbdd_true;
    if (b == mtbdd_false) return mtbdd_true;

    /* Maybe perform garbage collection */
    sylvan_gc_test();

    /* Count operation */
    sylvan_stats_count(MTBDD_LEQ);

    /* Check cache */
    MTBDD result;
    if (cache_get3(CACHE_MTBDD_LEQ, a, b, 0, &result)) {
        sylvan_stats_count(MTBDD_LEQ_CACHED);
        return result;
    }

    mtbddnode_t na = MTBDD_GETNODE(a);
    mtbddnode_t nb = MTBDD_GETNODE(b);
    int la = mtbddnode_isleaf(na);
    int lb = mtbddnode_isleaf(nb);

    if (la && lb) {
        uint64_t va = mtbddnode_getvalue(na);
        uint64_t vb = mtbddnode_getvalue(nb);

        if (mtbddnode_gettype(na) == 0 && mtbddnode_gettype(nb) == 0) {
            // type 0 = integer
            result = *(int64_t*)(&va) <= *(int64_t*)(&vb) ? mtbdd_true : mtbdd_false;
        } else if (mtbddnode_gettype(na) == 1 && mtbddnode_gettype(nb) == 1) {
            // type 1 = double
            double vva = *(double*)&va;
            double vvb = *(double*)&vb;
            result = vva <= vvb ? mtbdd_true : mtbdd_false;
        } else if (mtbddnode_gettype(na) == 2 && mtbddnode_gettype(nb) == 2) {
            // type 2 = fraction
            int64_t nom_a = (int32_t)(va>>32);
            int64_t nom_b = (int32_t)(vb>>32);
            uint64_t da = va&0xffffffff;
            uint64_t db = vb&0xffffffff;
            // equalize denominators
            uint32_t c = gcd(da, db);
            nom_a *= db/c;
            nom_b *= da/c;
            result = nom_a <= nom_b ? mtbdd_true : mtbdd_false;
        } else {
            assert(0); // failure
        }
    } else {
        /* Get top variable */
        uint32_t va = la ? 0xffffffff : mtbddnode_getvariable(na);
        uint32_t vb = lb ? 0xffffffff : mtbddnode_getvariable(nb);
        uint32_t var = va < vb ? va : vb;

        /* Get cofactors */
        MTBDD alow, ahigh, blow, bhigh;
        alow  = va == var ? node_getlow(a, na)  : a;
        ahigh = va == var ? node_gethigh(a, na) : a;
        blow  = vb == var ? node_getlow(b, nb)  : b;
        bhigh = vb == var ? node_gethigh(b, nb) : b;

        SPAWN(mtbdd_leq_rec, ahigh, bhigh, shortcircuit);
        result = CALL(mtbdd_leq_rec, alow, blow, shortcircuit);
        if (result != SYNC(mtbdd_leq_rec)) result = mtbdd_false;
    }

    if (result == mtbdd_false) *shortcircuit = 1;

    /* Store in cache */
    if (cache_put3(CACHE_MTBDD_LEQ, a, b, 0, result)) {
        sylvan_stats_count(MTBDD_LEQ_CACHEDPUT);
    }

    return result;
}

TASK_IMPL_2(MTBDD, mtbdd_leq, MTBDD, a, MTBDD, b)
{
    /* the implementation checks shortcircuit in every task and if the two
       MTBDDs are not equal module epsilon, then the computation tree quickly aborts */
    int shortcircuit = 0;
    return CALL(mtbdd_leq_rec, a, b, &shortcircuit);
}

/**
 * For two MTBDDs a, b, return mtbdd_true if all common assignments a(s) < b(s), mtbdd_false otherwise.
 * For domains not in a / b, assume True.
 */
TASK_3(MTBDD, mtbdd_less_rec, MTBDD, a, MTBDD, b, int*, shortcircuit)
{
    /* Check short circuit */
    if (*shortcircuit) return mtbdd_false;

    /* Check terminal case */
    if (a == b) return mtbdd_false;

    /* For partial functions, just return true */
    if (a == mtbdd_false) return mtbdd_true;
    if (b == mtbdd_false) return mtbdd_true;

    /* Maybe perform garbage collection */
    sylvan_gc_test();

    /* Count operation */
    sylvan_stats_count(MTBDD_LESS);

    /* Check cache */
    MTBDD result;
    if (cache_get3(CACHE_MTBDD_LESS, a, b, 0, &result)) {
        sylvan_stats_count(MTBDD_LESS_CACHED);
        return result;
    }

    mtbddnode_t na = MTBDD_GETNODE(a);
    mtbddnode_t nb = MTBDD_GETNODE(b);
    int la = mtbddnode_isleaf(na);
    int lb = mtbddnode_isleaf(nb);

    if (la && lb) {
        uint64_t va = mtbddnode_getvalue(na);
        uint64_t vb = mtbddnode_getvalue(nb);

        if (mtbddnode_gettype(na) == 0 && mtbddnode_gettype(nb) == 0) {
            // type 0 = integer
            result = *(int64_t*)(&va) < *(int64_t*)(&vb) ? mtbdd_true : mtbdd_false;
        } else if (mtbddnode_gettype(na) == 1 && mtbddnode_gettype(nb) == 1) {
            // type 1 = double
            double vva = *(double*)&va;
            double vvb = *(double*)&vb;
            result = vva < vvb ? mtbdd_true : mtbdd_false;
        } else if (mtbddnode_gettype(na) == 2 && mtbddnode_gettype(nb) == 2) {
            // type 2 = fraction
            int64_t nom_a = (int32_t)(va>>32);
            int64_t nom_b = (int32_t)(vb>>32);
            uint64_t da = va&0xffffffff;
            uint64_t db = vb&0xffffffff;
            // equalize denominators
            uint32_t c = gcd(da, db);
            nom_a *= db/c;
            nom_b *= da/c;
            result = nom_a < nom_b ? mtbdd_true : mtbdd_false;
        } else {
            assert(0); // failure
        }
    } else {
        /* Get top variable */
        uint32_t va = la ? 0xffffffff : mtbddnode_getvariable(na);
        uint32_t vb = lb ? 0xffffffff : mtbddnode_getvariable(nb);
        uint32_t var = va < vb ? va : vb;

        /* Get cofactors */
        MTBDD alow, ahigh, blow, bhigh;
        alow  = va == var ? node_getlow(a, na)  : a;
        ahigh = va == var ? node_gethigh(a, na) : a;
        blow  = vb == var ? node_getlow(b, nb)  : b;
        bhigh = vb == var ? node_gethigh(b, nb) : b;

        SPAWN(mtbdd_less_rec, ahigh, bhigh, shortcircuit);
        result = CALL(mtbdd_less_rec, alow, blow, shortcircuit);
        if (result != SYNC(mtbdd_less_rec)) result = mtbdd_false;
    }

    if (result == mtbdd_false) *shortcircuit = 1;

    /* Store in cache */
    if (cache_put3(CACHE_MTBDD_LESS, a, b, 0, result)) {
        sylvan_stats_count(MTBDD_LESS_CACHEDPUT);
    }

    return result;
}

TASK_IMPL_2(MTBDD, mtbdd_less, MTBDD, a, MTBDD, b)
{
    /* the implementation checks shortcircuit in every task and if the two
       MTBDDs are not equal module epsilon, then the computation tree quickly aborts */
    int shortcircuit = 0;
    return CALL(mtbdd_less_rec, a, b, &shortcircuit);
}

/**
 * For two MTBDDs a, b, return mtbdd_true if all common assignments a(s) >= b(s), mtbdd_false otherwise.
 * For domains not in a / b, assume True.
 */
TASK_3(MTBDD, mtbdd_geq_rec, MTBDD, a, MTBDD, b, int*, shortcircuit)
{
    /* Check short circuit */
    if (*shortcircuit) return mtbdd_false;

    /* Check terminal case */
    if (a == b) return mtbdd_true;

    /* For partial functions, just return true */
    if (a == mtbdd_false) return mtbdd_true;
    if (b == mtbdd_false) return mtbdd_true;

    /* Maybe perform garbage collection */
    sylvan_gc_test();

    /* Count operation */
    sylvan_stats_count(MTBDD_GEQ);

    /* Check cache */
    MTBDD result;
    if (cache_get3(CACHE_MTBDD_GEQ, a, b, 0, &result)) {
        sylvan_stats_count(MTBDD_GEQ_CACHED);
        return result;
    }

    mtbddnode_t na = MTBDD_GETNODE(a);
    mtbddnode_t nb = MTBDD_GETNODE(b);
    int la = mtbddnode_isleaf(na);
    int lb = mtbddnode_isleaf(nb);

    if (la && lb) {
        uint64_t va = mtbddnode_getvalue(na);
        uint64_t vb = mtbddnode_getvalue(nb);

        if (mtbddnode_gettype(na) == 0 && mtbddnode_gettype(nb) == 0) {
            // type 0 = integer
            result = *(int64_t*)(&va) >= *(int64_t*)(&vb) ? mtbdd_true : mtbdd_false;
        } else if (mtbddnode_gettype(na) == 1 && mtbddnode_gettype(nb) == 1) {
            // type 1 = double
            double vva = *(double*)&va;
            double vvb = *(double*)&vb;
            result = vva >= vvb ? mtbdd_true : mtbdd_false;
        } else if (mtbddnode_gettype(na) == 2 && mtbddnode_gettype(nb) == 2) {
            // type 2 = fraction
            int64_t nom_a = (int32_t)(va>>32);
            int64_t nom_b = (int32_t)(vb>>32);
            uint64_t da = va&0xffffffff;
            uint64_t db = vb&0xffffffff;
            // equalize denominators
            uint32_t c = gcd(da, db);
            nom_a *= db/c;
            nom_b *= da/c;
            result = nom_a >= nom_b ? mtbdd_true : mtbdd_false;
        } else {
            assert(0); // failure
        }
    } else {
        /* Get top variable */
        uint32_t va = la ? 0xffffffff : mtbddnode_getvariable(na);
        uint32_t vb = lb ? 0xffffffff : mtbddnode_getvariable(nb);
        uint32_t var = va < vb ? va : vb;

        /* Get cofactors */
        MTBDD alow, ahigh, blow, bhigh;
        alow  = va == var ? node_getlow(a, na)  : a;
        ahigh = va == var ? node_gethigh(a, na) : a;
        blow  = vb == var ? node_getlow(b, nb)  : b;
        bhigh = vb == var ? node_gethigh(b, nb) : b;

        SPAWN(mtbdd_geq_rec, ahigh, bhigh, shortcircuit);
        result = CALL(mtbdd_geq_rec, alow, blow, shortcircuit);
        if (result != SYNC(mtbdd_geq_rec)) result = mtbdd_false;
    }

    if (result == mtbdd_false) *shortcircuit = 1;

    /* Store in cache */
    if (cache_put3(CACHE_MTBDD_GEQ, a, b, 0, result)) {
        sylvan_stats_count(MTBDD_GEQ_CACHEDPUT);
    }

    return result;
}

TASK_IMPL_2(MTBDD, mtbdd_geq, MTBDD, a, MTBDD, b)
{
    /* the implementation checks shortcircuit in every task and if the two
       MTBDDs are not equal module epsilon, then the computation tree quickly aborts */
    int shortcircuit = 0;
    return CALL(mtbdd_geq_rec, a, b, &shortcircuit);
}

/**
 * For two MTBDDs a, b, return mtbdd_true if all common assignments a(s) > b(s), mtbdd_false otherwise.
 * For domains not in a / b, assume True.
 */
TASK_3(MTBDD, mtbdd_greater_rec, MTBDD, a, MTBDD, b, int*, shortcircuit)
{
    /* Check short circuit */
    if (*shortcircuit) return mtbdd_false;

    /* Check terminal case */
    if (a == b) return mtbdd_false;

    /* For partial functions, just return true */
    if (a == mtbdd_false) return mtbdd_true;
    if (b == mtbdd_false) return mtbdd_true;

    /* Maybe perform garbage collection */
    sylvan_gc_test();

    /* Count operation */
    sylvan_stats_count(MTBDD_GREATER);

    /* Check cache */
    MTBDD result;
    if (cache_get3(CACHE_MTBDD_GREATER, a, b, 0, &result)) {
        sylvan_stats_count(MTBDD_GREATER_CACHED);
        return result;
    }

    mtbddnode_t na = MTBDD_GETNODE(a);
    mtbddnode_t nb = MTBDD_GETNODE(b);
    int la = mtbddnode_isleaf(na);
    int lb = mtbddnode_isleaf(nb);

    if (la && lb) {
        uint64_t va = mtbddnode_getvalue(na);
        uint64_t vb = mtbddnode_getvalue(nb);

        if (mtbddnode_gettype(na) == 0 && mtbddnode_gettype(nb) == 0) {
            // type 0 = integer
            result = *(int64_t*)(&va) > *(int64_t*)(&vb) ? mtbdd_true : mtbdd_false;
        } else if (mtbddnode_gettype(na) == 1 && mtbddnode_gettype(nb) == 1) {
            // type 1 = double
            double vva = *(double*)&va;
            double vvb = *(double*)&vb;
            result = vva > vvb ? mtbdd_true : mtbdd_false;
        } else if (mtbddnode_gettype(na) == 2 && mtbddnode_gettype(nb) == 2) {
            // type 2 = fraction
            int64_t nom_a = (int32_t)(va>>32);
            int64_t nom_b = (int32_t)(vb>>32);
            uint64_t da = va&0xffffffff;
            uint64_t db = vb&0xffffffff;
            // equalize denominators
            uint32_t c = gcd(da, db);
            nom_a *= db/c;
            nom_b *= da/c;
            result = nom_a > nom_b ? mtbdd_true : mtbdd_false;
        } else {
            assert(0); // failure
        }
    } else {
        /* Get top variable */
        uint32_t va = la ? 0xffffffff : mtbddnode_getvariable(na);
        uint32_t vb = lb ? 0xffffffff : mtbddnode_getvariable(nb);
        uint32_t var = va < vb ? va : vb;

        /* Get cofactors */
        MTBDD alow, ahigh, blow, bhigh;
        alow  = va == var ? node_getlow(a, na)  : a;
        ahigh = va == var ? node_gethigh(a, na) : a;
        blow  = vb == var ? node_getlow(b, nb)  : b;
        bhigh = vb == var ? node_gethigh(b, nb) : b;

        SPAWN(mtbdd_greater_rec, ahigh, bhigh, shortcircuit);
        result = CALL(mtbdd_greater_rec, alow, blow, shortcircuit);
        if (result != SYNC(mtbdd_greater_rec)) result = mtbdd_false;
    }

    if (result == mtbdd_false) *shortcircuit = 1;

    /* Store in cache */
    if (cache_put3(CACHE_MTBDD_GREATER, a, b, 0, result)) {
        sylvan_stats_count(MTBDD_GREATER_CACHEDPUT);
    }

    return result;
}

TASK_IMPL_2(MTBDD, mtbdd_greater, MTBDD, a, MTBDD, b)
{
    /* the implementation checks shortcircuit in every task and if the two
       MTBDDs are not equal module epsilon, then the computation tree quickly aborts */
    int shortcircuit = 0;
    return CALL(mtbdd_greater_rec, a, b, &shortcircuit);
}

/**
 * Multiply <a> and <b>, and abstract variables <vars> using summation.
 * This is similar to the "and_exists" operation in BDDs.
 */
TASK_IMPL_3(MTBDD, mtbdd_and_abstract_plus, MTBDD, a, MTBDD, b, MTBDD, v)
{
    /* Check terminal case */
    if (v == mtbdd_true) return mtbdd_apply(a, b, TASK(mtbdd_op_times));
    MTBDD result = CALL(mtbdd_op_times, &a, &b);
    if (result != mtbdd_invalid) {
        mtbdd_refs_push(result);
        result = mtbdd_abstract(result, v, TASK(mtbdd_abstract_op_plus));
        mtbdd_refs_pop(1);
        return result;
    }

    /* Maybe perform garbage collection */
    sylvan_gc_test();

    /* Count operation */
    sylvan_stats_count(MTBDD_AND_ABSTRACT_PLUS);

    /* Check cache */
    if (cache_get3(CACHE_MTBDD_AND_ABSTRACT_PLUS, a, b, v, &result)) {
        sylvan_stats_count(MTBDD_AND_ABSTRACT_PLUS_CACHED);
        return result;
    }

    /* Now, v is not a constant, and either a or b is not a constant */

    /* Get top variable */
    int la = mtbdd_isleaf(a);
    int lb = mtbdd_isleaf(b);
    mtbddnode_t na = la ? 0 : MTBDD_GETNODE(a);
    mtbddnode_t nb = lb ? 0 : MTBDD_GETNODE(b);
    uint32_t va = la ? 0xffffffff : mtbddnode_getvariable(na);
    uint32_t vb = lb ? 0xffffffff : mtbddnode_getvariable(nb);
    uint32_t var = va < vb ? va : vb;

    mtbddnode_t nv = MTBDD_GETNODE(v);
    uint32_t vv = mtbddnode_getvariable(nv);

    if (vv < var) {
        /* Recursive, then abstract result */
        result = CALL(mtbdd_and_abstract_plus, a, b, node_gethigh(v, nv));
        mtbdd_refs_push(result);
        result = mtbdd_apply(result, result, TASK(mtbdd_op_plus));
        mtbdd_refs_pop(1);
    } else {
        /* Get cofactors */
        MTBDD alow, ahigh, blow, bhigh;
        alow  = (!la && va == var) ? node_getlow(a, na)  : a;
        ahigh = (!la && va == var) ? node_gethigh(a, na) : a;
        blow  = (!lb && vb == var) ? node_getlow(b, nb)  : b;
        bhigh = (!lb && vb == var) ? node_gethigh(b, nb) : b;

        if (vv == var) {
            /* Recursive, then abstract result */
            mtbdd_refs_spawn(SPAWN(mtbdd_and_abstract_plus, ahigh, bhigh, node_gethigh(v, nv)));
            MTBDD low = mtbdd_refs_push(CALL(mtbdd_and_abstract_plus, alow, blow, node_gethigh(v, nv)));
            MTBDD high = mtbdd_refs_push(mtbdd_refs_sync(SYNC(mtbdd_and_abstract_plus)));
            result = CALL(mtbdd_apply, low, high, TASK(mtbdd_op_plus));
            mtbdd_refs_pop(2);
        } else /* vv > v */ {
            /* Recursive, then create node */
            mtbdd_refs_spawn(SPAWN(mtbdd_and_abstract_plus, ahigh, bhigh, v));
            MTBDD low = mtbdd_refs_push(CALL(mtbdd_and_abstract_plus, alow, blow, v));
            MTBDD high = mtbdd_refs_sync(SYNC(mtbdd_and_abstract_plus));
            mtbdd_refs_pop(1);
            result = mtbdd_makenode(var, low, high);
        }
    }

    /* Store in cache */
    if (cache_put3(CACHE_MTBDD_AND_ABSTRACT_PLUS, a, b, v, result)) {
        sylvan_stats_count(MTBDD_AND_ABSTRACT_PLUS_CACHEDPUT);
    }

    return result;
}

/**
 * Multiply <a> and <b>, and abstract variables <vars> by taking the maximum.
 */
TASK_IMPL_3(MTBDD, mtbdd_and_abstract_max, MTBDD, a, MTBDD, b, MTBDD, v)
{
    /* Check terminal case */
    if (v == mtbdd_true) return mtbdd_apply(a, b, TASK(mtbdd_op_times));
    MTBDD result = CALL(mtbdd_op_times, &a, &b);
    if (result != mtbdd_invalid) {
        mtbdd_refs_push(result);
        result = mtbdd_abstract(result, v, TASK(mtbdd_abstract_op_max));
        mtbdd_refs_pop(1);
        return result;
    }

    /* Now, v is not a constant, and either a or b is not a constant */

    /* Get top variable */
    int la = mtbdd_isleaf(a);
    int lb = mtbdd_isleaf(b);
    mtbddnode_t na = la ? 0 : MTBDD_GETNODE(a);
    mtbddnode_t nb = lb ? 0 : MTBDD_GETNODE(b);
    uint32_t va = la ? 0xffffffff : mtbddnode_getvariable(na);
    uint32_t vb = lb ? 0xffffffff : mtbddnode_getvariable(nb);
    uint32_t var = va < vb ? va : vb;

    mtbddnode_t nv = MTBDD_GETNODE(v);
    uint32_t vv = mtbddnode_getvariable(nv);

    while (vv < var) {
        /* we can skip variables, because max(r,r) = r */
        v = node_gethigh(v, nv);
        if (v == mtbdd_true) return mtbdd_apply(a, b, TASK(mtbdd_op_times));
        nv = MTBDD_GETNODE(v);
        vv = mtbddnode_getvariable(nv);
    }

    /* Maybe perform garbage collection */
    sylvan_gc_test();

    /* Count operation */
    sylvan_stats_count(MTBDD_AND_ABSTRACT_MAX);

    /* Check cache */
    if (cache_get3(CACHE_MTBDD_AND_ABSTRACT_MAX, a, b, v, &result)) {
        sylvan_stats_count(MTBDD_AND_ABSTRACT_MAX_CACHED);
        return result;
    }

    /* Get cofactors */
    MTBDD alow, ahigh, blow, bhigh;
    alow  = (!la && va == var) ? node_getlow(a, na)  : a;
    ahigh = (!la && va == var) ? node_gethigh(a, na) : a;
    blow  = (!lb && vb == var) ? node_getlow(b, nb)  : b;
    bhigh = (!lb && vb == var) ? node_gethigh(b, nb) : b;

    if (vv == var) {
        /* Recursive, then abstract result */
        mtbdd_refs_spawn(SPAWN(mtbdd_and_abstract_max, ahigh, bhigh, node_gethigh(v, nv)));
        MTBDD low = mtbdd_refs_push(CALL(mtbdd_and_abstract_max, alow, blow, node_gethigh(v, nv)));
        MTBDD high = mtbdd_refs_push(mtbdd_refs_sync(SYNC(mtbdd_and_abstract_max)));
        result = CALL(mtbdd_apply, low, high, TASK(mtbdd_op_max));
        mtbdd_refs_pop(2);
    } else /* vv > v */ {
        /* Recursive, then create node */
        mtbdd_refs_spawn(SPAWN(mtbdd_and_abstract_max, ahigh, bhigh, v));
        MTBDD low = mtbdd_refs_push(CALL(mtbdd_and_abstract_max, alow, blow, v));
        MTBDD high = mtbdd_refs_sync(SYNC(mtbdd_and_abstract_max));
        mtbdd_refs_pop(1);
        result = mtbdd_makenode(var, low, high);
    }

    /* Store in cache */
    if (cache_put3(CACHE_MTBDD_AND_ABSTRACT_MAX, a, b, v, result)) {
        sylvan_stats_count(MTBDD_AND_ABSTRACT_MAX_CACHEDPUT);
    }

    return result;
}

/**
 * Calculate the support of a MTBDD, i.e. the cube of all variables that appear in the MTBDD nodes.
 */
TASK_IMPL_1(MTBDD, mtbdd_support, MTBDD, dd)
{
    /* Terminal case */
    if (mtbdd_isleaf(dd)) return mtbdd_true;

    /* Maybe perform garbage collection */
    sylvan_gc_test();

    /* Count operation */
    sylvan_stats_count(BDD_SUPPORT);

    /* Check cache */
    MTBDD result;
    if (cache_get3(CACHE_MTBDD_SUPPORT, dd, 0, 0, &result)) {
        sylvan_stats_count(BDD_SUPPORT_CACHED);
        return result;
    }

    /* Recursive calls */
    mtbddnode_t n = MTBDD_GETNODE(dd);
    mtbdd_refs_spawn(SPAWN(mtbdd_support, node_getlow(dd, n)));
    MTBDD high = mtbdd_refs_push(CALL(mtbdd_support, node_gethigh(dd, n)));
    MTBDD low = mtbdd_refs_push(mtbdd_refs_sync(SYNC(mtbdd_support)));

    /* Compute result */
    result = mtbdd_makenode(mtbddnode_getvariable(n), mtbdd_false, sylvan_and(low, high));
    mtbdd_refs_pop(2);

    /* Write to cache */
    if (cache_put3(CACHE_MTBDD_SUPPORT, dd, 0, 0, result)) {
        sylvan_stats_count(BDD_SUPPORT_CACHEDPUT);
    }

    return result;
}

/**
 * Function composition, for each node with variable <key> which has a <key,value> pair in <map>,
 * replace the node by the result of mtbdd_ite(<value>, <high>, <low>).
 * Each <value> in <map> must be a Boolean MTBDD.
 */
TASK_IMPL_2(MTBDD, mtbdd_compose, MTBDD, a, MTBDDMAP, map)
{
    /* Terminal case */
    if (mtbdd_isleaf(a) || mtbdd_map_isempty(map)) return a;

    /* Determine top level */
    mtbddnode_t n = MTBDD_GETNODE(a);
    uint32_t v = mtbddnode_getvariable(n);

    /* Find in map */
    while (mtbdd_map_key(map) < v) {
        map = mtbdd_map_next(map);
        if (mtbdd_map_isempty(map)) return a;
    }

    /* Perhaps execute garbage collection */
    sylvan_gc_test();

    /* Count operation */
    sylvan_stats_count(MTBDD_COMPOSE);

    /* Check cache */
    MTBDD result;
    if (cache_get3(CACHE_MTBDD_COMPOSE, a, map, 0, &result)) {
        sylvan_stats_count(MTBDD_COMPOSE_CACHED);
        return result;
    }

    /* Recursive calls */
    mtbdd_refs_spawn(SPAWN(mtbdd_compose, node_getlow(a, n), map));
    MTBDD high = mtbdd_refs_push(CALL(mtbdd_compose, node_gethigh(a, n), map));
    MTBDD low = mtbdd_refs_push(mtbdd_refs_sync(SYNC(mtbdd_compose)));

    /* Calculate result */
    MTBDD r = mtbdd_map_key(map) == v ? mtbdd_map_value(map) : mtbdd_makenode(v, mtbdd_false, mtbdd_true);
    mtbdd_refs_push(r);
    result = CALL(mtbdd_ite, r, high, low);
    mtbdd_refs_pop(3);

    /* Store in cache */
    if (cache_put3(CACHE_MTBDD_COMPOSE, a, map, 0, result)) {
        sylvan_stats_count(MTBDD_COMPOSE_CACHEDPUT);
    }

    return result;
}

/**
 * Compute minimum leaf in the MTBDD (for Integer, Double, Rational MTBDDs)
 */
TASK_IMPL_1(MTBDD, mtbdd_minimum, MTBDD, a)
{
    /* Check terminal case */
    if (a == mtbdd_false) return mtbdd_false;
    mtbddnode_t na = MTBDD_GETNODE(a);
    if (mtbddnode_isleaf(na)) return a;

    /* Maybe perform garbage collection */
    sylvan_gc_test();

    /* Count operation */
    sylvan_stats_count(MTBDD_MINIMUM);

    /* Check cache */
    MTBDD result;
    if (cache_get3(CACHE_MTBDD_MINIMUM, a, 0, 0, &result)) {
        sylvan_stats_count(MTBDD_MINIMUM_CACHED);
        return result;
    }

    /* Call recursive */
    SPAWN(mtbdd_minimum, node_getlow(a, na));
    MTBDD high = CALL(mtbdd_minimum, node_gethigh(a, na));
    MTBDD low = SYNC(mtbdd_minimum);

    /* Determine lowest */
    mtbddnode_t nl = MTBDD_GETNODE(low);
    mtbddnode_t nh = MTBDD_GETNODE(high);

    if (mtbddnode_gettype(nl) == 0 && mtbddnode_gettype(nh) == 0) {
        result = mtbdd_getint64(low) < mtbdd_getint64(high) ? low : high;
    } else if (mtbddnode_gettype(nl) == 1 && mtbddnode_gettype(nh) == 1) {
        result = mtbdd_getdouble(low) < mtbdd_getdouble(high) ? low : high;
    } else if (mtbddnode_gettype(nl) == 2 && mtbddnode_gettype(nh) == 2) {
        // type 2 = fraction
        int64_t nom_l = mtbdd_getnumer(low);
        int64_t nom_h = mtbdd_getnumer(high);
        uint64_t denom_l = mtbdd_getdenom(low);
        uint64_t denom_h = mtbdd_getdenom(high);
        // equalize denominators
        uint32_t c = gcd(denom_l, denom_h);
        nom_l *= denom_h/c;
        nom_h *= denom_l/c;
        result = nom_l < nom_h ? low : high;
    } else {
        assert(0); // failure TODO: how to handle custom types?
    }

    /* Store in cache */
    if (cache_put3(CACHE_MTBDD_MINIMUM, a, 0, 0, result)) {
        sylvan_stats_count(MTBDD_MINIMUM_CACHEDPUT);
    }

    return result;
}

/**
 * Compute maximum leaf in the MTBDD (for Integer, Double, Rational MTBDDs)
 */
TASK_IMPL_1(MTBDD, mtbdd_maximum, MTBDD, a)
{
    /* Check terminal case */
    if (a == mtbdd_false) return mtbdd_false;
    mtbddnode_t na = MTBDD_GETNODE(a);
    if (mtbddnode_isleaf(na)) return a;

    /* Maybe perform garbage collection */
    sylvan_gc_test();

    /* Count operation */
    sylvan_stats_count(MTBDD_MAXIMUM);

    /* Check cache */
    MTBDD result;
    if (cache_get3(CACHE_MTBDD_MAXIMUM, a, 0, 0, &result)) {
        sylvan_stats_count(MTBDD_MAXIMUM_CACHED);
        return result;
    }

    /* Call recursive */
    SPAWN(mtbdd_maximum, node_getlow(a, na));
    MTBDD high = CALL(mtbdd_maximum, node_gethigh(a, na));
    MTBDD low = SYNC(mtbdd_maximum);

    /* Determine highest */
    mtbddnode_t nl = MTBDD_GETNODE(low);
    mtbddnode_t nh = MTBDD_GETNODE(high);

    if (mtbddnode_gettype(nl) == 0 && mtbddnode_gettype(nh) == 0) {
        result = mtbdd_getint64(low) > mtbdd_getint64(high) ? low : high;
    } else if (mtbddnode_gettype(nl) == 1 && mtbddnode_gettype(nh) == 1) {
        result = mtbdd_getdouble(low) > mtbdd_getdouble(high) ? low : high;
    } else if (mtbddnode_gettype(nl) == 2 && mtbddnode_gettype(nh) == 2) {
        // type 2 = fraction
        int64_t nom_l = mtbdd_getnumer(low);
        int64_t nom_h = mtbdd_getnumer(high);
        uint64_t denom_l = mtbdd_getdenom(low);
        uint64_t denom_h = mtbdd_getdenom(high);
        // equalize denominators
        uint32_t c = gcd(denom_l, denom_h);
        nom_l *= denom_h/c;
        nom_h *= denom_l/c;
        result = nom_l > nom_h ? low : high;
    } else {
        assert(0); // failure
    }

    /* Store in cache */
    if (cache_put3(CACHE_MTBDD_MAXIMUM, a, 0, 0, result)) {
        sylvan_stats_count(MTBDD_MAXIMUM_CACHEDPUT);
    }

    return result;
}

/**
 * Calculate the number of satisfying variable assignments according to <variables>.
 */
TASK_IMPL_2(double, mtbdd_satcount, MTBDD, dd, size_t, nvars)
{
    /* Trivial cases */
    if (dd == mtbdd_false) return 0.0;

    if (mtbdd_isleaf(dd)) {
        // test if 0
        mtbddnode_t dd_node = MTBDD_GETNODE(dd);
        if (dd != mtbdd_true) {
            if (mtbddnode_gettype(dd_node) == 0 && mtbdd_getint64(dd) == 0) return 0.0;
            else if (mtbddnode_gettype(dd_node) == 1 && mtbdd_getdouble(dd) == 0.0) return 0.0;
            else if (mtbddnode_gettype(dd_node) == 2 && mtbdd_getvalue(dd) == 1) return 0.0;
        }
        return powl(2.0L, nvars);
    }

    /* Perhaps execute garbage collection */
    sylvan_gc_test();

    union {
        double d;
        uint64_t s;
    } hack;

    /* Consult cache */
    if (cache_get3(CACHE_BDD_SATCOUNT, dd, 0, nvars, &hack.s)) {
        sylvan_stats_count(BDD_SATCOUNT_CACHED);
        return hack.d;
    }

    SPAWN(mtbdd_satcount, mtbdd_gethigh(dd), nvars-1);
    double low = CALL(mtbdd_satcount, mtbdd_getlow(dd), nvars-1);
    hack.d = low + SYNC(mtbdd_satcount);

    if (cache_put3(CACHE_BDD_SATCOUNT, dd, 0, nvars, hack.s)) {
        sylvan_stats_count(BDD_SATCOUNT_CACHEDPUT);
    }

    return hack.d;
}

MTBDD
mtbdd_enum_first(MTBDD dd, MTBDD variables, uint8_t *arr, mtbdd_enum_filter_cb filter_cb)
{
    if (dd == mtbdd_false) {
        // the leaf dd is skipped
        return mtbdd_false;
    } else if (mtbdd_isleaf(dd)) {
        // a leaf for which the filter returns 0 is skipped
        if (filter_cb != NULL && filter_cb(dd) == 0) return mtbdd_false;
        // ok, we have a leaf that is not skipped, go for it!
        while (variables != mtbdd_true) {
            *arr++ = 2;
            variables = mtbdd_gethigh(variables);
        }
        return dd;
    } else if (variables == mtbdd_true) {
        // in the case of partial evaluation... treat like a leaf
        if (filter_cb != NULL && filter_cb(dd) == 0) return mtbdd_false;
        return dd;
    } else {
        // if variables == true, then dd must be a leaf. But then this line is unreachable.
        // if this assertion fails, then <variables> is not the support of <dd>.
        assert(variables != mtbdd_true);

        // get next variable from <variables>
        uint32_t v = mtbdd_getvar(variables);
        variables = mtbdd_gethigh(variables);

        // check if MTBDD is on this variable
        mtbddnode_t n = MTBDD_GETNODE(dd);
        if (mtbddnode_getvariable(n) != v) {
            *arr = 2;
            return mtbdd_enum_first(dd, variables, arr+1, filter_cb);
        }

        // first maybe follow low
        MTBDD res = mtbdd_enum_first(node_getlow(dd, n), variables, arr+1, filter_cb);
        if (res != mtbdd_false) {
            *arr = 0;
            return res;
        }

        // if not low, try following high
        res = mtbdd_enum_first(node_gethigh(dd, n), variables, arr+1, filter_cb);
        if (res != mtbdd_false) {
            *arr = 1;
            return res;
        }
        
        // we've tried low and high, return false
        return mtbdd_false;
    }
}

MTBDD
mtbdd_enum_next(MTBDD dd, MTBDD variables, uint8_t *arr, mtbdd_enum_filter_cb filter_cb)
{
    if (mtbdd_isleaf(dd)) {
        // we find the leaf in 'enum_next', then we've seen it before...
        return mtbdd_false;
    } else if (variables == mtbdd_true) {
        // in the case of partial evaluation... treat like a leaf
        return mtbdd_false;
    } else {
        // if variables == true, then dd must be a leaf. But then this line is unreachable.
        // if this assertion fails, then <variables> is not the support of <dd>.
        assert(variables != mtbdd_true);

        variables = mtbdd_gethigh(variables);

        if (*arr == 0) {
            // previous was low
            mtbddnode_t n = MTBDD_GETNODE(dd);
            MTBDD res = mtbdd_enum_next(node_getlow(dd, n), variables, arr+1, filter_cb);
            if (res != mtbdd_false) {
                return res;
            } else {
                // try to find new in high branch
                res = mtbdd_enum_first(node_gethigh(dd, n), variables, arr+1, filter_cb);
                if (res != mtbdd_false) {
                    *arr = 1;
                    return res;
                } else {
                    return mtbdd_false;
                }
            }
        } else if (*arr == 1) {
            // previous was high
            mtbddnode_t n = MTBDD_GETNODE(dd);
            return mtbdd_enum_next(node_gethigh(dd, n), variables, arr+1, filter_cb);
        } else {
            // previous was either
            return mtbdd_enum_next(dd, variables, arr+1, filter_cb);
        }
    }
}

MTBDD
mtbdd_enum_all_first(MTBDD dd, MTBDD variables, uint8_t *arr, mtbdd_enum_filter_cb filter_cb)
{
    if (dd == mtbdd_false) {
        // the leaf False is skipped
        return mtbdd_false;
    } else if (variables == mtbdd_true) {
        // if this assertion fails, then <variables> is not the support of <dd>.
        // actually, remove this check to allow for "partial" enumeration
        // assert(mtbdd_isleaf(dd));
        // for _first, just return the leaf; there is nothing to set, though.
        if (filter_cb != NULL && filter_cb(dd) == 0) return mtbdd_false;
        return dd;
    } else if (mtbdd_isleaf(dd)) {
        // a leaf for which the filter returns 0 is skipped
        if (filter_cb != NULL && filter_cb(dd) == 0) return mtbdd_false;
        // for all remaining variables, set to 0
        while (variables != mtbdd_true) {
            *arr++ = 0;
            variables = mtbdd_gethigh(variables);
        }
        return dd;
    } else {
        // get next variable from <variables>
        mtbddnode_t nv = MTBDD_GETNODE(variables);
        variables = node_gethigh(variables, nv);

        // get cofactors
        mtbddnode_t ndd = MTBDD_GETNODE(dd);
        MTBDD low, high;
        if (mtbddnode_getvariable(ndd) == mtbddnode_getvariable(nv)) {
            low = node_getlow(dd, ndd);
            high = node_gethigh(dd, ndd);
        } else {
            low = high = dd;
        }

        // first maybe follow low
        MTBDD res = mtbdd_enum_all_first(low, variables, arr+1, filter_cb);
        if (res != mtbdd_false) {
            *arr = 0;
            return res;
        }

        // if not low, try following high
        res = mtbdd_enum_all_first(high, variables, arr+1, filter_cb);
        if (res != mtbdd_false) {
            *arr = 1;
            return res;
        }

        // we've tried low and high, return false
        return mtbdd_false;
    }
}

MTBDD
mtbdd_enum_all_next(MTBDD dd, MTBDD variables, uint8_t *arr, mtbdd_enum_filter_cb filter_cb)
{
    if (dd == mtbdd_false) {
        // the leaf False is skipped
        return mtbdd_false;
    } else if (variables == mtbdd_true) {
        // if this assertion fails, then <variables> is not the support of <dd>.
        // actually, remove this check to allow for "partial" enumeration
        // assert(mtbdd_isleaf(dd));
        // no next if there are no variables
        return mtbdd_false;
    } else {
        // get next variable from <variables>
        mtbddnode_t nv = MTBDD_GETNODE(variables);
        variables = node_gethigh(variables, nv);

        // filter leaf (if leaf) or get cofactors (if not leaf)
        mtbddnode_t ndd = MTBDD_GETNODE(dd);
        MTBDD low, high;
        if (mtbdd_isleaf(dd)) {
            // a leaf for which the filter returns 0 is skipped
            if (filter_cb != NULL && filter_cb(dd) == 0) return mtbdd_false;
            low = high = dd;
        } else {
            // get cofactors
            if (mtbddnode_getvariable(ndd) == mtbddnode_getvariable(nv)) {
                low = node_getlow(dd, ndd);
                high = node_gethigh(dd, ndd);
            } else {
                low = high = dd;
            }
        }

        // try recursive next first
        if (*arr == 0) {
            MTBDD res = mtbdd_enum_all_next(low, variables, arr+1, filter_cb);
            if (res != mtbdd_false) return res;
        } else if (*arr == 1) {
            return mtbdd_enum_all_next(high, variables, arr+1, filter_cb);
            // if *arr was 1 and _next returns False, return False
        } else {
            // the array is invalid...
            assert(*arr == 0 || *arr == 1);
            return mtbdd_invalid;  // in Release mode, the assertion is empty code
        }

        // previous was low, try following high
        MTBDD res = mtbdd_enum_all_first(high, variables, arr+1, filter_cb);
        if (res == mtbdd_false) return mtbdd_false;

        // succesful, set arr
        *arr = 1;
        return res;
    }
}

/**
 * Given a MTBDD <dd>, call <cb> with context <context> for every unique path in <dd> ending in leaf <leaf>.
 *
 * Usage:
 * VOID_TASK_3(cb, mtbdd_enum_trace_t, trace, MTBDD, leaf, void*, context) { ... do something ... }
 * mtbdd_enum_par(dd, cb, context);
 */
VOID_TASK_4(mtbdd_enum_par_do, MTBDD, dd, mtbdd_enum_cb, cb, void*, context, mtbdd_enum_trace_t, trace)
{
    if (mtbdd_isleaf(dd)) {
        WRAP(cb, trace, dd, context);
        return;
    }

    mtbddnode_t ndd = MTBDD_GETNODE(dd);
    uint32_t var = mtbddnode_getvariable(ndd);

    struct mtbdd_enum_trace t0 = (struct mtbdd_enum_trace){trace, var, 0};
    struct mtbdd_enum_trace t1 = (struct mtbdd_enum_trace){trace, var, 1};
    SPAWN(mtbdd_enum_par_do, node_getlow(dd, ndd), cb, context, &t0);
    CALL(mtbdd_enum_par_do, node_gethigh(dd, ndd), cb, context, &t1);
    SYNC(mtbdd_enum_par_do);
}

VOID_TASK_IMPL_3(mtbdd_enum_par, MTBDD, dd, mtbdd_enum_cb, cb, void*, context)
{
    CALL(mtbdd_enum_par_do, dd, cb, context, NULL);
}

/**
 * Function composition after partial evaluation.
 *
 * Given a function F(X) = f, compute the composition F'(X) = g(f) for every assignment to X.
 * All variables X in <vars> must appear before all variables in f and g(f).
 *
 * Usage:
 * TASK_2(MTBDD, g, MTBDD, in) { ... return g of <in> ... }
 * MTBDD x_vars = ...;  // the cube of variables x
 * MTBDD result = mtbdd_eval_compose(dd, x_vars, TASK(g));
 */
TASK_IMPL_3(MTBDD, mtbdd_eval_compose, MTBDD, dd, MTBDD, vars, mtbdd_eval_compose_cb, cb)
{
    /* Maybe perform garbage collection */
    sylvan_gc_test();

    /* Count operation */
    sylvan_stats_count(MTBDD_EVAL_COMPOSE);

    /* Check cache */
    MTBDD result;
    if (cache_get3(CACHE_MTBDD_EVAL_COMPOSE, dd, vars, (size_t)cb, &result)) {
        sylvan_stats_count(MTBDD_EVAL_COMPOSE_CACHED);
        return result;
    }

    if (mtbdd_isleaf(dd) || vars == mtbdd_true) {
        /* Apply */
        result = WRAP(cb, dd);
    } else {
        /* Get top variable in dd */
        mtbddnode_t ndd = MTBDD_GETNODE(dd);
        uint32_t var = mtbddnode_getvariable(ndd);

        /* Check if <var> is in <vars> */
        mtbddnode_t nvars = MTBDD_GETNODE(vars);
        uint32_t vv = mtbddnode_getvariable(nvars);

        /* Search/forward <vars> */
        MTBDD _vars = vars;
        while (vv < var) {
            _vars = node_gethigh(_vars, nvars);
            if (_vars == mtbdd_true) break;
            nvars = MTBDD_GETNODE(_vars);
            vv = mtbddnode_getvariable(nvars);
        }

        if (_vars == mtbdd_true) {
            /* Apply */
            result = WRAP(cb, dd);
        } else {
            /* If this fails, then there are variables in f/g BEFORE vars, which breaks functionality. */
            assert(vv == var);
            if (vv != var) return mtbdd_invalid;

            /* Get cofactors */
            MTBDD ddlow = node_getlow(dd, ndd);
            MTBDD ddhigh = node_gethigh(dd, ndd);

            /* Recursive */
            _vars = node_gethigh(_vars, nvars);
            mtbdd_refs_spawn(SPAWN(mtbdd_eval_compose, ddhigh, _vars, cb));
            MTBDD low = mtbdd_refs_push(CALL(mtbdd_eval_compose, ddlow, _vars, cb));
            MTBDD high = mtbdd_refs_sync(SYNC(mtbdd_eval_compose));
            mtbdd_refs_pop(1);
            result = mtbdd_makenode(var, low, high);
        }
    }

    /* Store in cache */
    if (cache_put3(CACHE_MTBDD_EVAL_COMPOSE, dd, vars, (size_t)cb, result)) {
        sylvan_stats_count(MTBDD_EVAL_COMPOSE_CACHEDPUT);
    }

    return result;
}

/**
 * Helper function for recursive unmarking
 */
static void
mtbdd_unmark_rec(MTBDD mtbdd)
{
    mtbddnode_t n = MTBDD_GETNODE(mtbdd);
    if (!mtbddnode_getmark(n)) return;
    mtbddnode_setmark(n, 0);
    if (mtbddnode_isleaf(n)) return;
    mtbdd_unmark_rec(mtbddnode_getlow(n));
    mtbdd_unmark_rec(mtbddnode_gethigh(n));
}

/**
 * Count number of leaves in MTBDD
 */

static size_t
mtbdd_leafcount_mark(MTBDD mtbdd)
{
    if (mtbdd == mtbdd_true) return 0; // do not count true/false leaf
    if (mtbdd == mtbdd_false) return 0; // do not count true/false leaf
    mtbddnode_t n = MTBDD_GETNODE(mtbdd);
    if (mtbddnode_getmark(n)) return 0;
    mtbddnode_setmark(n, 1);
    if (mtbddnode_isleaf(n)) return 1; // count leaf as 1
    return mtbdd_leafcount_mark(mtbddnode_getlow(n)) + mtbdd_leafcount_mark(mtbddnode_gethigh(n));
}

size_t
mtbdd_leafcount_more(const MTBDD *mtbdds, size_t count)
{
    size_t result = 0, i;
    for (i=0; i<count; i++) result += mtbdd_leafcount_mark(mtbdds[i]);
    for (i=0; i<count; i++) mtbdd_unmark_rec(mtbdds[i]);
    return result;
}

/**
 * Count number of nodes in MTBDD
 */

static size_t
mtbdd_nodecount_mark(MTBDD mtbdd)
{
    mtbddnode_t n = MTBDD_GETNODE(mtbdd);
    if (mtbddnode_getmark(n)) return 0;
    mtbddnode_setmark(n, 1);
    if (mtbddnode_isleaf(n)) return 1; // count leaf as 1
    return 1 + mtbdd_nodecount_mark(mtbddnode_getlow(n)) + mtbdd_nodecount_mark(mtbddnode_gethigh(n));
}

size_t
mtbdd_nodecount_more(const MTBDD *mtbdds, size_t count)
{
    size_t result = 0, i;
    for (i=0; i<count; i++) result += mtbdd_nodecount_mark(mtbdds[i]);
    for (i=0; i<count; i++) mtbdd_unmark_rec(mtbdds[i]);
    return result;
}

TASK_2(int, mtbdd_test_isvalid_rec, MTBDD, dd, uint32_t, parent_var)
{
    // check if True/False leaf
    if (dd == mtbdd_true || dd == mtbdd_false) return 1;

    // check if index is in array
    uint64_t index = dd & (~mtbdd_complement);
    assert(index > 1 && index < nodes->table_size);
    if (index <= 1 || index >= nodes->table_size) return 0;

    // check if marked
    int marked = llmsset_is_marked(nodes, index);
    assert(marked);
    if (marked == 0) return 0;

    // check if leaf
    mtbddnode_t n = MTBDD_GETNODE(dd);
    if (mtbddnode_isleaf(n)) return 1; // we're fine

    // check variable order
    uint32_t var = mtbddnode_getvariable(n);
    assert(var > parent_var);
    if (var <= parent_var) return 0;

    // check cache
    uint64_t result;
    if (cache_get3(CACHE_BDD_ISBDD, dd, 0, 0, &result)) {
        sylvan_stats_count(BDD_ISBDD_CACHED);
        return result;
    }

    // check recursively
    SPAWN(mtbdd_test_isvalid_rec, node_getlow(dd, n), var);
    result = (uint64_t)CALL(mtbdd_test_isvalid_rec, node_gethigh(dd, n), var);
    if (!SYNC(mtbdd_test_isvalid_rec)) result = 0;

    // put in cache and return result
    if (cache_put3(CACHE_BDD_ISBDD, dd, 0, 0, result)) {
        sylvan_stats_count(BDD_ISBDD_CACHEDPUT);
    }

    return result;
}

TASK_IMPL_1(int, mtbdd_test_isvalid, MTBDD, dd)
{
    // check if True/False leaf
    if (dd == mtbdd_true || dd == mtbdd_false) return 1;

    // check if index is in array
    uint64_t index = dd & (~mtbdd_complement);
    assert(index > 1 && index < nodes->table_size);
    if (index <= 1 || index >= nodes->table_size) return 0;

    // check if marked
    int marked = llmsset_is_marked(nodes, index);
    assert(marked);
    if (marked == 0) return 0;

    // check if leaf
    mtbddnode_t n = MTBDD_GETNODE(dd);
    if (mtbddnode_isleaf(n)) return 1; // we're fine

    // check recursively
    uint32_t var = mtbddnode_getvariable(n);
    SPAWN(mtbdd_test_isvalid_rec, node_getlow(dd, n), var);
    int result = CALL(mtbdd_test_isvalid_rec, node_gethigh(dd, n), var);
    if (!SYNC(mtbdd_test_isvalid_rec)) result = 0;
    return result;
}

/**
 * Write a text representation of a leaf to the given file.
 */
void
mtbdd_fprint_leaf(FILE *out, MTBDD leaf)
{
    char buf[64];
    char *ptr = mtbdd_leaf_to_str(leaf, buf, 64);
    if (ptr != NULL) {
        fputs(ptr, out);
        if (ptr != buf) free(ptr);
    }
}

/**
 * Write a text representation of a leaf to stdout.
 */
void
mtbdd_print_leaf(MTBDD leaf)
{
    mtbdd_fprint_leaf(stdout, leaf);
}

/**
 * Obtain the textual representation of a leaf.
 * The returned result is either equal to the given <buf> (if the results fits)
 * or to a newly allocated array (with malloc).
 */
char *
mtbdd_leaf_to_str(MTBDD leaf, char *buf, size_t buflen)
{
    mtbddnode_t n = MTBDD_GETNODE(leaf);
    uint32_t type = mtbddnode_gettype(n);
    uint64_t value = mtbddnode_getvalue(n);
    int complement = MTBDD_HASMARK(leaf) ? 1 : 0;

    return sylvan_mt_to_str(complement, type, value, buf, buflen);
}

/**
 * Export to .dot file
 */

static void
mtbdd_fprintdot_rec(FILE *out, MTBDD mtbdd)
{
    mtbddnode_t n = MTBDD_GETNODE(mtbdd); // also works for mtbdd_false
    if (mtbddnode_getmark(n)) return;
    mtbddnode_setmark(n, 1);

    if (mtbdd == mtbdd_true || mtbdd == mtbdd_false) {
        fprintf(out, "0 [shape=box, style=filled, label=\"F\"];\n");
    } else if (mtbddnode_isleaf(n)) {
        fprintf(out, "%" PRIu64 " [shape=box, style=filled, label=\"", MTBDD_STRIPMARK(mtbdd));
        mtbdd_fprint_leaf(out, mtbdd);
        fprintf(out, "\"];\n");
    } else {
        fprintf(out, "%" PRIu64 " [label=\"%" PRIu32 "\"];\n",
                MTBDD_STRIPMARK(mtbdd), mtbddnode_getvariable(n));

        mtbdd_fprintdot_rec(out, mtbddnode_getlow(n));
        mtbdd_fprintdot_rec(out, mtbddnode_gethigh(n));

        fprintf(out, "%" PRIu64 " -> %" PRIu64 " [style=dashed];\n",
                MTBDD_STRIPMARK(mtbdd), mtbddnode_getlow(n));
        fprintf(out, "%" PRIu64 " -> %" PRIu64 " [style=solid dir=both arrowtail=%s];\n",
                MTBDD_STRIPMARK(mtbdd), MTBDD_STRIPMARK(mtbddnode_gethigh(n)),
                mtbddnode_getcomp(n) ? "dot" : "none");
    }
}

void
mtbdd_fprintdot(FILE *out, MTBDD mtbdd)
{
    fprintf(out, "digraph \"DD\" {\n");
    fprintf(out, "graph [dpi = 300];\n");
    fprintf(out, "center = true;\n");
    fprintf(out, "edge [dir = forward];\n");
    fprintf(out, "root [style=invis];\n");
    fprintf(out, "root -> %" PRIu64 " [style=solid dir=both arrowtail=%s];\n",
            MTBDD_STRIPMARK(mtbdd), MTBDD_HASMARK(mtbdd) ? "dot" : "none");

    mtbdd_fprintdot_rec(out, mtbdd);
    mtbdd_unmark_rec(mtbdd);

    fprintf(out, "}\n");
}

/**
 * Export to .dot file, but do not display complement edges. Expand instead.
 */

static void
mtbdd_fprintdot_nc_rec(FILE *out, MTBDD mtbdd)
{
    mtbddnode_t n = MTBDD_GETNODE(mtbdd); // also works for mtbdd_false
    if (mtbddnode_getmark(n)) return;
    mtbddnode_setmark(n, 1);

    if (mtbdd == mtbdd_true) {
        fprintf(out, "%" PRIu64 " [shape=box, style=filled, label=\"T\"];\n", mtbdd);
    } else if (mtbdd == mtbdd_false) {
        fprintf(out, "0 [shape=box, style=filled, label=\"F\"];\n");
    } else if (mtbddnode_isleaf(n)) {
        fprintf(out, "%" PRIu64 " [shape=box, style=filled, label=\"", mtbdd);
        mtbdd_fprint_leaf(out, mtbdd);
        fprintf(out, "\"];\n");
    } else {
        fprintf(out, "%" PRIu64 " [label=\"%" PRIu32 "\"];\n", mtbdd, mtbddnode_getvariable(n));

        mtbdd_fprintdot_nc_rec(out, mtbddnode_getlow(n));
        mtbdd_fprintdot_nc_rec(out, mtbddnode_gethigh(n));

        fprintf(out, "%" PRIu64 " -> %" PRIu64 " [style=dashed];\n", mtbdd, node_getlow(mtbdd, n));
        fprintf(out, "%" PRIu64 " -> %" PRIu64 " [style=solid];\n", mtbdd, node_gethigh(mtbdd, n));
    }
}

void
mtbdd_fprintdot_nc(FILE *out, MTBDD mtbdd)
{
    fprintf(out, "digraph \"DD\" {\n");
    fprintf(out, "graph [dpi = 300];\n");
    fprintf(out, "center = true;\n");
    fprintf(out, "edge [dir = forward];\n");
    fprintf(out, "root [style=invis];\n");
    fprintf(out, "root -> %" PRIu64 " [style=solid];\n", mtbdd);

    mtbdd_fprintdot_nc_rec(out, mtbdd);
    mtbdd_unmark_rec(mtbdd);

    fprintf(out, "}\n");
}

/**
 * Generate SHA2 structural hashes.
 * Hashes are independent of location.
 * Mainly useful for debugging purposes.
 */
static void
mtbdd_sha2_rec(MTBDD dd, SHA256_CTX *ctx)
{
    if (dd == mtbdd_true || dd == mtbdd_false) {
        SHA256_Update(ctx, (void*)&dd, sizeof(MTBDD));
        return;
    }

    mtbddnode_t node = MTBDD_GETNODE(dd);
    if (mtbddnode_getmark(node) == 0) {
        mtbddnode_setmark(node, 1);
        if (mtbddnode_isleaf(node)) {
            uint32_t type = mtbddnode_gettype(node);
            SHA256_Update(ctx, (void*)&type, sizeof(uint32_t));
            uint64_t value = mtbddnode_getvalue(node);
            value = sylvan_mt_hash(type, value, value);
            SHA256_Update(ctx, (void*)&value, sizeof(uint64_t));
        } else {
            uint32_t level = mtbddnode_getvariable(node);
            if (MTBDD_STRIPMARK(mtbddnode_gethigh(node))) level |= 0x80000000;
            SHA256_Update(ctx, (void*)&level, sizeof(uint32_t));
            mtbdd_sha2_rec(mtbddnode_gethigh(node), ctx);
            mtbdd_sha2_rec(mtbddnode_getlow(node), ctx);
        }
    }
}

void
mtbdd_printsha(MTBDD dd)
{
    mtbdd_fprintsha(stdout, dd);
}

void
mtbdd_fprintsha(FILE *f, MTBDD dd)
{
    char buf[80];
    mtbdd_getsha(dd, buf);
    fprintf(f, "%s", buf);
}

void
mtbdd_getsha(MTBDD dd, char *target)
{
    SHA256_CTX ctx;
    SHA256_Init(&ctx);
    mtbdd_sha2_rec(dd, &ctx);
    if (dd != mtbdd_true && dd != mtbdd_false) mtbdd_unmark_rec(dd);
    SHA256_End(&ctx, target);
}

/**
 * Implementation of visitor operations
 */

VOID_TASK_IMPL_4(mtbdd_visit_seq, MTBDD, dd, mtbdd_visit_pre_cb, pre_cb, mtbdd_visit_post_cb, post_cb, void*, ctx)
{
    int children = 1;
    if (pre_cb != NULL) children = WRAP(pre_cb, dd, ctx);
    if (children && !mtbdd_isleaf(dd)) {
        CALL(mtbdd_visit_seq, mtbdd_getlow(dd), pre_cb, post_cb, ctx);
        CALL(mtbdd_visit_seq, mtbdd_gethigh(dd), pre_cb, post_cb, ctx);
    }
    if (post_cb != NULL) WRAP(post_cb, dd, ctx);
}

VOID_TASK_IMPL_4(mtbdd_visit_par, MTBDD, dd, mtbdd_visit_pre_cb, pre_cb, mtbdd_visit_post_cb, post_cb, void*, ctx)
{
    int children = 1;
    if (pre_cb != NULL) children = WRAP(pre_cb, dd, ctx);
    if (children && !mtbdd_isleaf(dd)) {
        SPAWN(mtbdd_visit_par, mtbdd_getlow(dd), pre_cb, post_cb, ctx);
        CALL(mtbdd_visit_par, mtbdd_gethigh(dd), pre_cb, post_cb, ctx);
        SYNC(mtbdd_visit_par);
    }
    if (post_cb != NULL) WRAP(post_cb, dd, ctx);
}

/**
 * Writing MTBDD files using a skiplist as a backend
 */

TASK_2(int, mtbdd_writer_add_visitor_pre, MTBDD, dd, sylvan_skiplist_t, sl)
{
    if (mtbdd_isleaf(dd)) return 0;
    return sylvan_skiplist_get(sl, MTBDD_STRIPMARK(dd)) == 0 ? 1 : 0;
}

VOID_TASK_2(mtbdd_writer_add_visitor_post, MTBDD, dd, sylvan_skiplist_t, sl)
{
    if (dd == mtbdd_true || dd == mtbdd_false) return;
    sylvan_skiplist_assign_next(sl, MTBDD_STRIPMARK(dd));
}

sylvan_skiplist_t
mtbdd_writer_start()
{
    size_t sl_size = nodes->table_size > 0x7fffffff ? 0x7fffffff : nodes->table_size;
    return sylvan_skiplist_alloc(sl_size);
}

VOID_TASK_IMPL_2(mtbdd_writer_add, sylvan_skiplist_t, sl, MTBDD, dd)
{
    mtbdd_visit_seq(dd, (mtbdd_visit_pre_cb)TASK(mtbdd_writer_add_visitor_pre), (mtbdd_visit_post_cb)TASK(mtbdd_writer_add_visitor_post), (void*)sl);
}

void
mtbdd_writer_writebinary(FILE *out, sylvan_skiplist_t sl)
{
    size_t nodecount = sylvan_skiplist_count(sl);
    fwrite(&nodecount, sizeof(size_t), 1, out);
    for (size_t i=1; i<=nodecount; i++) {
        MTBDD dd = sylvan_skiplist_getr(sl, i);

        mtbddnode_t n = MTBDD_GETNODE(dd);
        if (mtbddnode_isleaf(n)) {
            /* write leaf */
            fwrite(n, sizeof(struct mtbddnode), 1, out);
            uint32_t type = mtbddnode_gettype(n);
            uint64_t value = mtbddnode_getvalue(n);
            sylvan_mt_write_binary(type, value, out);
        } else {
            struct mtbddnode node;
            MTBDD low = sylvan_skiplist_get(sl, mtbddnode_getlow(n));
            MTBDD high = mtbddnode_gethigh(n);
            high = MTBDD_TRANSFERMARK(high, sylvan_skiplist_get(sl, MTBDD_STRIPMARK(high)));
            mtbddnode_makenode(&node, mtbddnode_getvariable(n), low, high);
            fwrite(&node, sizeof(struct mtbddnode), 1, out);
        }
    }
}

uint64_t
mtbdd_writer_get(sylvan_skiplist_t sl, MTBDD dd)
{
    return MTBDD_TRANSFERMARK(dd, sylvan_skiplist_get(sl, MTBDD_STRIPMARK(dd)));
}

void
mtbdd_writer_end(sylvan_skiplist_t sl)
{
    sylvan_skiplist_free(sl);
}

VOID_TASK_IMPL_3(mtbdd_writer_tobinary, FILE *, out, MTBDD *, dds, int, count)
{
    sylvan_skiplist_t sl = mtbdd_writer_start();

    for (int i=0; i<count; i++) {
        CALL(mtbdd_writer_add, sl, dds[i]);
    }

    mtbdd_writer_writebinary(out, sl);

    fwrite(&count, sizeof(int), 1, out);
    
    for (int i=0; i<count; i++) {
        uint64_t v = mtbdd_writer_get(sl, dds[i]);
        fwrite(&v, sizeof(uint64_t), 1, out);
    }

    mtbdd_writer_end(sl);
}

void
mtbdd_writer_writetext(FILE *out, sylvan_skiplist_t sl)
{
    fprintf(out, "[\n");
    size_t nodecount = sylvan_skiplist_count(sl);
    for (size_t i=1; i<=nodecount; i++) {
        MTBDD dd = sylvan_skiplist_getr(sl, i);

        mtbddnode_t n = MTBDD_GETNODE(dd);
        if (mtbddnode_isleaf(n)) {
            /* serialize leaf, does not support customs yet */
            fprintf(out, "  leaf(%zu,%u,\"", i, mtbddnode_gettype(n));
            mtbdd_fprint_leaf(out, MTBDD_STRIPMARK(dd));
            fprintf(out, "\"),\n");
        } else {
            MTBDD low = sylvan_skiplist_get(sl, mtbddnode_getlow(n));
            MTBDD high = mtbddnode_gethigh(n);
            high = MTBDD_TRANSFERMARK(high, sylvan_skiplist_get(sl, MTBDD_STRIPMARK(high)));
            fprintf(out, "  node(%zu,%u,%zu,%s%zu),\n", i, mtbddnode_getvariable(n), (size_t)low, MTBDD_HASMARK(high)?"~":"", (size_t)MTBDD_STRIPMARK(high));
        }
    }

    fprintf(out, "]");
}

VOID_TASK_IMPL_3(mtbdd_writer_totext, FILE *, out, MTBDD *, dds, int, count)
{
    sylvan_skiplist_t sl = mtbdd_writer_start();

    for (int i=0; i<count; i++) {
        CALL(mtbdd_writer_add, sl, dds[i]);
    }

    mtbdd_writer_writetext(out, sl);

    fprintf(out, ",[");
    
    for (int i=0; i<count; i++) {
        uint64_t v = mtbdd_writer_get(sl, dds[i]);
        fprintf(out, "%s%zu,", MTBDD_HASMARK(v)?"~":"", (size_t)MTBDD_STRIPMARK(v));
    }

    fprintf(out, "]\n");

    mtbdd_writer_end(sl);
}

/**
 * Reading a file earlier written with mtbdd_writer_writebinary
 * Returns an array with the conversion from stored identifier to MTBDD
 * This array is allocated with malloc and must be freed afterwards.
 * This method does not support custom leaves.
 */
TASK_IMPL_1(uint64_t*, mtbdd_reader_readbinary, FILE*, in)
{
    size_t nodecount;
    if (fread(&nodecount, sizeof(size_t), 1, in) != 1) {
        return NULL;
    }

    uint64_t *arr = malloc(sizeof(uint64_t)*(nodecount+1));
    arr[0] = 0;
    for (size_t i=1; i<=nodecount; i++) {
        struct mtbddnode node;
        if (fread(&node, sizeof(struct mtbddnode), 1, in) != 1) {
            free(arr);
            return NULL;
        }

        if (mtbddnode_isleaf(&node)) {
            /* serialize leaf */
            uint32_t type = mtbddnode_gettype(&node);
            uint64_t value = mtbddnode_getvalue(&node);
            sylvan_mt_read_binary(type, &value, in);
            arr[i] = mtbdd_makeleaf(type, value);
        } else {
            MTBDD low = arr[mtbddnode_getlow(&node)];
            MTBDD high = mtbddnode_gethigh(&node);
            high = MTBDD_TRANSFERMARK(high, arr[MTBDD_STRIPMARK(high)]);
            arr[i] = mtbdd_makenode(mtbddnode_getvariable(&node), low, high);
        }
    }

    return arr;
}

/**
 * Retrieve the MTBDD of the given stored identifier.
 */
MTBDD
mtbdd_reader_get(uint64_t* arr, uint64_t identifier)
{
    return MTBDD_TRANSFERMARK(identifier, arr[MTBDD_STRIPMARK(identifier)]);
}

/**
 * Free the allocated translation array
 */
void
mtbdd_reader_end(uint64_t *arr)
{
    free(arr);
}

/**
 * Reading a file earlier written with mtbdd_writer_tobinary
 */
TASK_IMPL_3(int, mtbdd_reader_frombinary, FILE*, in, MTBDD*, dds, int, count)
{
    uint64_t *arr = CALL(mtbdd_reader_readbinary, in);
    if (arr == NULL) return -1;

    /* Read stored count */
    int actual_count;
    if (fread(&actual_count, sizeof(int), 1, in) != 1) {
        mtbdd_reader_end(arr);
        return -1;
    }

    /* If actual count does not agree with given count, abort */
    if (actual_count != count) {
        mtbdd_reader_end(arr);
        return -1;
    }
    
    /* Read every stored identifier, and translate to MTBDD */
    for (int i=0; i<count; i++) {
        uint64_t v;
        if (fread(&v, sizeof(uint64_t), 1, in) != 1) {
            mtbdd_reader_end(arr);
            return -1;
        }
        dds[i] = mtbdd_reader_get(arr, v);
    }

    mtbdd_reader_end(arr);
    return 0;
}

/**
 * Implementation of variable sets, i.e., cubes of (positive) variables.
 */

/**
 * Create a set of variables, represented as the conjunction of (positive) variables.
 */
MTBDD
mtbdd_set_from_array(uint32_t* arr, size_t length)
{
    if (length == 0) return mtbdd_true;
    else if (length == 1) return mtbdd_makenode(*arr, mtbdd_false, mtbdd_true);
    else return mtbdd_set_add(mtbdd_fromarray(arr+1, length-1), *arr);
}

/**
 * Write all variables in a variable set to the given array.
 * The array must be sufficiently large.
 */
void
mtbdd_set_to_array(MTBDD set, uint32_t *arr)
{
    while (set != mtbdd_true) {
        mtbddnode_t n = MTBDD_GETNODE(set);
        *arr++ = mtbddnode_getvariable(n);
        set = node_gethigh(set, n);
    }
}

/**
 * Add the variable <var> to <set>.
 */
MTBDD
mtbdd_set_add(MTBDD set, uint32_t var)
{
    if (set == mtbdd_true) return mtbdd_makenode(var, mtbdd_false, mtbdd_true);

    mtbddnode_t set_node = MTBDD_GETNODE(set);
    uint32_t set_var = mtbddnode_getvariable(set_node);
    if (var < set_var) return mtbdd_makenode(var, mtbdd_false, set);
    else if (set_var == var) return set;
    else {
        MTBDD sub = mtbddnode_followhigh(set, set_node);
        MTBDD res = mtbdd_set_add(sub, var);
        res = sub == res ? set : mtbdd_makenode(set_var, mtbdd_false, res);
        return res;
    }
}

/**
 * Remove the variable <var> from <set>.
 */
MTBDD
mtbdd_set_remove(MTBDD set, uint32_t var)
{
    if (set == mtbdd_true) return mtbdd_true;

    mtbddnode_t set_node = MTBDD_GETNODE(set);
    uint32_t set_var = mtbddnode_getvariable(set_node);
    if (var < set_var) return set;
    else if (set_var == var) return mtbddnode_followhigh(set, set_node);
    else {
        MTBDD sub = mtbddnode_followhigh(set, set_node);
        MTBDD res = mtbdd_set_remove(sub, var);
        res = sub == res ? set : mtbdd_makenode(set_var, mtbdd_false, res);
        return res;
    }
}

/**
 * Remove variables in <set2> from <set1>.
 */
TASK_IMPL_2(MTBDD, mtbdd_set_minus, MTBDD, set1, MTBDD, set2)
{
    if (set1 == mtbdd_true) return mtbdd_true;
    if (set2 == mtbdd_true) return set1;
    if (set1 == set2) return mtbdd_true;

    mtbddnode_t set1_node = MTBDD_GETNODE(set1);
    mtbddnode_t set2_node = MTBDD_GETNODE(set2);
    uint32_t set1_var = mtbddnode_getvariable(set1_node);
    uint32_t set2_var = mtbddnode_getvariable(set2_node);

    if (set1_var == set2_var) {
        return mtbdd_set_minus(mtbddnode_followhigh(set1, set1_node), mtbddnode_followhigh(set2, set2_node));
    }

    if (set1_var > set2_var) {
        return mtbdd_set_minus(set1, mtbddnode_followhigh(set2, set2_node));
    }

    /* set1_var < set2_var */
    MTBDD sub = mtbddnode_followhigh(set1, set1_node);
    MTBDD res = mtbdd_set_minus(sub, set2);
    return res == sub ? set1 : mtbdd_makenode(set1_var, mtbdd_false, res);
}

/**
 * Return 1 if <set> contains <var>, 0 otherwise.
 */
int
mtbdd_set_contains(MTBDD set, uint32_t var)
{
    while (set != mtbdd_true) {
        mtbddnode_t n = MTBDD_GETNODE(set);
        uint32_t v = mtbddnode_getvariable(n);
        if (v == var) return 1;
        if (v > var) return 0;
        set = node_gethigh(set, n);
    }
    return 0;
}

/**
 * Compute the number of variables in a given set of variables.
 */
size_t
mtbdd_set_count(MTBDD set)
{
    size_t result = 0;
    while (set != mtbdd_true) {
        result++;
        set = mtbdd_gethigh(set);
    }
    return result;
}

/**
 * Sanity check if the given MTBDD is a conjunction of positive variables,
 * and if all nodes are marked in the nodes table (detects violations after garbage collection).
 */
void
mtbdd_test_isset(MTBDD set)
{
    while (set != mtbdd_true) {
        assert(set != mtbdd_false);
        assert(llmsset_is_marked(nodes, set));
        mtbddnode_t n = MTBDD_GETNODE(set);
        assert(node_getlow(set, n) == mtbdd_false);
        set = node_gethigh(set, n);
    }
}

/**
 * Return 1 if the map contains the key, 0 otherwise.
 */
int
mtbdd_map_contains(MTBDDMAP map, uint32_t key)
{
    while (!mtbdd_map_isempty(map)) {
        mtbddnode_t n = MTBDD_GETNODE(map);
        uint32_t k = mtbddnode_getvariable(n);
        if (k == key) return 1;
        if (k > key) return 0;
        map = node_getlow(map, n);
    }

    return 0;
}

/**
 * Retrieve the number of keys in the map.
 */
size_t
mtbdd_map_count(MTBDDMAP map)
{
    size_t r = 0;

    while (!mtbdd_map_isempty(map)) {
        r++;
        map = mtbdd_map_next(map);
    }

    return r;
}

/**
 * Add the pair <key,value> to the map, overwrites if key already in map.
 */
MTBDDMAP
mtbdd_map_add(MTBDDMAP map, uint32_t key, MTBDD value)
{
    if (mtbdd_map_isempty(map)) {
        return mtbdd_makemapnode(key, mtbdd_map_empty(), value);
    }

    mtbddnode_t n = MTBDD_GETNODE(map);
    uint32_t k = mtbddnode_getvariable(n);

    if (k < key) {
        // add recursively and rebuild tree
        MTBDDMAP low = mtbdd_map_add(node_getlow(map, n), key, value);
        return mtbdd_makemapnode(k, low, node_gethigh(map, n));
    } else if (k > key) {
        return mtbdd_makemapnode(key, map, value);
    } else {
        // replace old
        return mtbdd_makemapnode(key, node_getlow(map, n), value);
    }
}

/**
 * Add all values from map2 to map1, overwrites if key already in map1.
 */
MTBDDMAP
mtbdd_map_update(MTBDDMAP map1, MTBDDMAP map2)
{
    if (mtbdd_map_isempty(map1)) return map2;
    if (mtbdd_map_isempty(map2)) return map1;

    mtbddnode_t n1 = MTBDD_GETNODE(map1);
    mtbddnode_t n2 = MTBDD_GETNODE(map2);
    uint32_t k1 = mtbddnode_getvariable(n1);
    uint32_t k2 = mtbddnode_getvariable(n2);

    MTBDDMAP result;
    if (k1 < k2) {
        MTBDDMAP low = mtbdd_map_update(node_getlow(map1, n1), map2);
        result = mtbdd_makemapnode(k1, low, node_gethigh(map1, n1));
    } else if (k1 > k2) {
        MTBDDMAP low = mtbdd_map_update(map1, node_getlow(map2, n2));
        result = mtbdd_makemapnode(k2, low, node_gethigh(map2, n2));
    } else {
        MTBDDMAP low = mtbdd_map_update(node_getlow(map1, n1), node_getlow(map2, n2));
        result = mtbdd_makemapnode(k2, low, node_gethigh(map2, n2));
    }

    return result;
}

/**
 * Remove the key <key> from the map and return the result
 */
MTBDDMAP
mtbdd_map_remove(MTBDDMAP map, uint32_t key)
{
    if (mtbdd_map_isempty(map)) return map;

    mtbddnode_t n = MTBDD_GETNODE(map);
    uint32_t k = mtbddnode_getvariable(n);

    if (k < key) {
        MTBDDMAP low = mtbdd_map_remove(node_getlow(map, n), key);
        return mtbdd_makemapnode(k, low, node_gethigh(map, n));
    } else if (k > key) {
        return map;
    } else {
        return node_getlow(map, n);
    }
}

/**
 * Remove all keys in the cube <variables> from the map and return the result
 */
MTBDDMAP
mtbdd_map_removeall(MTBDDMAP map, MTBDD variables)
{
    if (mtbdd_map_isempty(map)) return map;
    if (variables == mtbdd_true) return map;

    mtbddnode_t n1 = MTBDD_GETNODE(map);
    mtbddnode_t n2 = MTBDD_GETNODE(variables);
    uint32_t k1 = mtbddnode_getvariable(n1);
    uint32_t k2 = mtbddnode_getvariable(n2);

    if (k1 < k2) {
        MTBDDMAP low = mtbdd_map_removeall(node_getlow(map, n1), variables);
        return mtbdd_makemapnode(k1, low, node_gethigh(map, n1));
    } else if (k1 > k2) {
        return mtbdd_map_removeall(map, node_gethigh(variables, n2));
    } else {
        return mtbdd_map_removeall(node_getlow(map, n1), node_gethigh(variables, n2));
    }
}

// ===================================================================================================
//
//  Extension by Richard van Dijk  / Sebastiaan Brand of existing API for matrix / vector operations 
//
// ===================================================================================================

/**
 * Util functions for dynamically allocation of memory for arrays
 */

int
allocate_matrix_array(MatArr_t ***W_arr, int n)
{
    if(n < 0)
        return 1;

    *W_arr = (MatArr_t **)malloc((1 << n) * sizeof(MatArr_t *)); // row pointers

    if(*W_arr == NULL)
        return 2;

    for(int i=0; i < (1 << n); i++) {

        (*W_arr)[i] = (MatArr_t *)malloc((1 << n) * sizeof(MatArr_t)); // all elements

        if((*W_arr)[i] == NULL)
            return 3;

        for(int j=0; j < (1 << n); j++)
            (*W_arr)[i][j] = 0.0;
    }

    return 0;
}

int
free_matrix_array(MatArr_t **W_arr, int n)
{
    if(n<0)
        return 1;

    // free all mallocs
    for(int i=0; i < (1 << n); i++) {
        free(W_arr[i]);
    }

    free(W_arr);

    return 0;
}

int
print_vector_array(VecArr_t *v_arr, int n)
{
    if(n<0)
        return 1;

    for(int row=0; row < (1 << n); row++)
        printf("v_arr[%d] = %lf\n", row, v_arr[row]);

    return 0;
}

int
print_matrix_array(MatArr_t **W_arr, int n)
{
    if(n<0)
        return 1;

    for(int row=0; row < (1 << n); row++)
        for(int column=0; column < (1 << n); column++)
            printf("W_arr[%d][%d] = %lf\n", row, column, W_arr[row][column]);

    return 0;
}

/**
 * Utility functions:
 * 
 * Convert a matrix array M[row][col] into a MTBDD.
 * 
 * Convert a vector array v[row] into a MTBDD.
 * 
 * Convert a MTBDD into a matrix array.
 * 
 * Convert a MTBDD into a vector array.
 * 
 * The mode of the conversions refers to how the leafs are 
 * filled with the array values.
 * 
 * This can be row wise (transpose mode) or column wise.
 * 
 * Suppose M[row][col] is a matrix:
 * 
 *      M[0][0]     M[0][1]
 *      M[1][0]     M[1][1]
 * 
 * The transpose of M[row][col] is M[col][row]
 * 
 *      M[0][0]     M[1][0]
 *      M[0][1]     M[1][1]
 * 
 * Then MTBDD is column wise:
 * 
 *                          x0
 *                x1                   x1
 *
 *        M[0][0]    M[0][1]   M[1][0]     M[1][1]
 *
 * MTBDD is row wise (= transpose):
 * 
 *                          x0
 *                x1                   x1
 *
 *        M[0][0]    M[1][0]   M[0][1]     M[1][1]
 * 
 */

MTBDD vector_array_to_mtbdd(VecArr_t *v_arr, int n, row_column_mode_t mode)
{
    if(v_arr == NULL)
        return MTBDD_ZERO;

    if(n < 0)
        return MTBDD_ZERO;
    
    if(n == 0)
        return mtbdd_double(v_arr[0]);

    if(n == 1 && (mode == COLUMN_WISE_MODE || mode == ALTERNATE_COLUMN_FIRST_WISE_MODE)) {
        return mtbdd_makenode(1, mtbdd_double(v_arr[0]), mtbdd_double(v_arr[1]));
    }

    if(n == 1 && (mode == ROW_WISE_MODE || mode == ALTERNATE_ROW_FIRST_WISE_MODE)) {

        return mtbdd_makenode(0, mtbdd_double(v_arr[0]), mtbdd_double(v_arr[1]));
    }

    // TODO: if n > 1, recursive, divide in four parts, alternate_column_first

    return MTBDD_ZERO;
}

MTBDD matrix_array_to_mtbdd(MatArr_t **M_arr, int n, row_column_mode_t mode)
{
    if(M_arr == NULL)
        return MTBDD_ZERO;

    if(n < 0)
        return MTBDD_ZERO;

    if(n == 0) 
        return mtbdd_double(M_arr[0][0]);

    if(n == 1 && (mode == COLUMN_WISE_MODE || mode == ALTERNATE_COLUMN_FIRST_WISE_MODE)) {

        MTBDD column0 = mtbdd_makenode(1, mtbdd_double(M_arr[0][0]), mtbdd_double(M_arr[1][0]));
        MTBDD column1 = mtbdd_makenode(1, mtbdd_double(M_arr[0][1]), mtbdd_double(M_arr[1][1]));

        return mtbdd_makenode(0, column0, column1);
    }

    if(n == 1 && (mode == ROW_WISE_MODE || mode == ALTERNATE_ROW_FIRST_WISE_MODE)) {

        MTBDD row0 = mtbdd_makenode(1, mtbdd_double(M_arr[0][0]), mtbdd_double(M_arr[0][1]));
        MTBDD row1 = mtbdd_makenode(1, mtbdd_double(M_arr[1][0]), mtbdd_double(M_arr[1][1]));

        return mtbdd_makenode(0, row0, row1);
    }

    // TODO: if n > 1, recursive, divide in four parts

    return MTBDD_ZERO;
}

void mtbdd_to_vector_array(MTBDD v, int n, row_column_mode_t mode, VecArr_t *w)
{
    MatArr_t **W_arr = NULL;
    allocate_matrix_array(&W_arr, n);

    mtbdd_to_matrix_array(v, n, mode, W_arr);

    for(int row=0; row < (1 << n); row++)
        w[row] = W_arr[row][0];

    free_matrix_array(W_arr, n);

    return;
}

void mtbdd_to_matrix_array(MTBDD M, int n, row_column_mode_t mode, MatArr_t **W) // TODO: refactor to make more compact
{
    //
    // The leaf can be reached by the composite row * 2^n + column
    //
    // Example: 
    //
    //   row = 3, column = 5, n = 3. 
    //
    // That corresponds to a matrix of 2^3 x 2^3 = 8 x 8, or in C: "MatArr_t W[8][8]". 
    //
    // The row and column index = {0,1,..,7} = {0, ..., ((2^3)-1)}, W[row][column].
    //
    //  W[3][5] = W[011][100] = M[011100]
    //
    // By starting traversing through M from the root, first go to node 011. That is 2^3 deep.
    //
    // Under this node you find a mtbdd with depth 2^3 with leafs corresponding to W[][column].
    //
    // So, traversing further on with node[100] = getlow(getlow(gethigh(node))) you reach the leaf
    //
    // with as value equal to W[3][5].
    //
    //  W[3][5] = W[011][100] = leaf(node[100]) = leaf(M[011100])
    //

    if(M == MTBDD_ZERO)
        return;

    if(n <= 0)
        return;

    // f(c0,r0,c1,r1) = W[r0r1][c0c1] or f(r0,c0,r1,c1) = W[r0r1][c0c1]
    if(mode == ALTERNATE_COLUMN_FIRST_WISE_MODE || mode == ALTERNATE_ROW_FIRST_WISE_MODE) {

        for(int index=0; index < (1 << (2 * n)); index++) {
            
            MTBDD node = M;
            
            bool turn_for_column = true;

            int row = 0;
            int column = 0;

            int bit_row = (1 << (n-1)) - 1;
            int bit_column = (1 << (n-1)) - 1;

            for(int bit=((2 * n) - 1); bit >= 0; bit--) {

                if((index & (1 << bit)) == 0) {

                    if(mtbdd_getlow(node) != MTBDD_ZERO && mtbdd_isleaf(node) == 0) 
                        node = mtbdd_getlow(node);

                    if(turn_for_column) {
                        column += (0 << bit_column); // TODO: can be removed
                        bit_column -= 1;
                        turn_for_column = false;
                    }

                    else {
                        row += (0 << bit_row); // TODO: can be removed
                        bit_row -= 1;
                        turn_for_column = true;
                    }

                    //printf("node = %ld, index = %d, getlow -> bit=%d, row=%d, column=%d\n", node, index, bit, row, column);

                }
                else {

                    if(mtbdd_gethigh(node) != MTBDD_ZERO && mtbdd_isleaf(node) == 0) 
                        node = mtbdd_gethigh(node);

                    if(turn_for_column) {
                        column += (1 << bit_column);
                        bit_column -= 1;
                        turn_for_column = false;
                    }

                    else {
                        row += (1 << bit_row);
                        bit_row -= 1;
                        turn_for_column = true;
                    }

                    //printf("node = %ld, index = %d, getlow -> bit=%d, row=%d, column=%d\n", node, index, bit, row, column);

                }
            }

            //printf("row=%d, column=%d\n", row, column);

            if(mode == ALTERNATE_COLUMN_FIRST_WISE_MODE)
                W[row][column] = mtbdd_getdouble(node);

            if(mode == ALTERNATE_ROW_FIRST_WISE_MODE)
                W[column][row] = mtbdd_getdouble(node);
        }

        return;
    }

    // f(r0,r1,c0,c1) or f(c0,c1,r0,r1) 
    for(int row=0; row < (1 << n); row++) {

        MTBDD node = M;

        // Traverse through M from root to row node
        for(int bit=(n-1); bit >= 0; bit--) {

            if((row & (1 << bit)) == 0) {
                if(mtbdd_getlow(node) != MTBDD_ZERO && mtbdd_isleaf(node) == 0) 
                    node = mtbdd_getlow(node);
                
                printf("row=%d getlow()\n", row);
            }
            else {
                if(mtbdd_gethigh(node) != MTBDD_ZERO && mtbdd_isleaf(node) == 0) 
                    node = mtbdd_gethigh(node);
                
                //printf("row=%d gethigh()\n", row);
            }
        }

        MTBDD row_node = node;

        for(int column=0; column < (1 << n); column++) {

            node = row_node;

            // Traverse through M from node to column node
            for(int bit=(n-1); bit >= 0; bit--) {

                if((column & (1 << bit)) == 0) {
                    if(mtbdd_getlow(node) != MTBDD_ZERO && mtbdd_isleaf(node) == 0) 
                        node = mtbdd_getlow(node);
                    
                    //printf("column=%d getlow()\n", column);
                }
                else {
                    if(mtbdd_gethigh(node) != MTBDD_ZERO && mtbdd_isleaf(node) == 0) 
                        node = mtbdd_gethigh(node);
                    
                    //printf("column=%d gethigh()\n", column);
                }
            }

            if(mode == COLUMN_WISE_MODE)
                W[row][column] = mtbdd_getdouble(node);

            if(mode == ROW_WISE_MODE)
                W[column][row] = mtbdd_getdouble(node);
        }
    }

    return;
}

/**
 * Matrix . Vector multiplication in MTBDD domain
 * 
 * Computes M.v = w for a 2^n x 2^n matrix M[row][col], and a 2^n vector v[col].
 * 
 * Results in an array of a 2^n vector w[col].
 * 
 */
void array_matrix_vector_product(MatArr_t **M, VecArr_t *v, int n, VecArr_t *w)
{
    // Convert to MTBDD
    MTBDD M_col = matrix_array_to_mtbdd(M, n, COLUMN_WISE_MODE);
    MTBDD v_row = vector_array_to_mtbdd(v, n, ROW_WISE_MODE);

    // Multiply in MTBDD domain
    MTBDD w_col = mtbdd_matvec_mult_alt(M_col, v_row, n); // TODO: take none alt

    mtbdd_to_vector_array(w_col, n, COLUMN_WISE_MODE, w);

    return;
}

/**
 * Matrix . Matrix multiplication in MTBDD domain
 * 
 * Computes M.V = W for two 2^n x 2^n matrices M[row][col] and V[row][col].
 * 
 * Results in an array of a 2^n x 2^n matrix W[row][col]. 
 * 
 */
void array_matrix_matrix_product(MatArr_t **M1, MatArr_t **M2, int n, MatArr_t **W)
{
    // Convert to MTBDD
    MTBDD M1_col = matrix_array_to_mtbdd(M1, n, COLUMN_WISE_MODE);
    MTBDD M2_row = matrix_array_to_mtbdd(M2, n, ROW_WISE_MODE);

    // Multiply in MTBDD domain
    MTBDD W_col = mtbdd_matmat_mult_alt(M1_col, M2_row, n); // TODO: take none alt

    mtbdd_to_matrix_array(W_col, n, COLUMN_WISE_MODE, W);
    
    return;
}

/**
 * 
 * Util functions for matrix multiplication
 * 
 */

MTBDD mtbdd_is_result_in_cache_3(int function, uint64_t num1, uint64_t num2, uint64_t num3)
{
    MTBDD result = MTBDD_ZERO;

    if(!cache_get3(function, num1, num2, num3, &result)) {
        printf("Error in cache_get4()\n");
        return result; //exit(0);
    }

    return result;
}

void mtbdd_put_result_in_cache_3(int function, uint64_t num1, uint64_t num2, uint64_t num3, uint64_t result)
{
    if(!cache_put3(function, num1, num2, num3, result)) {
        printf("Error in cache_put4()\n");
        return; //exit(0);
    }
    return;
}

MTBDD mtbdd_is_result_in_cache_4(int function, uint64_t num1, uint64_t num2, uint64_t num3, uint64_t num4)
{
    MTBDD result = MTBDD_ZERO;

    if(!cache_get4(function, num1, num2, num3, num4, &result)) {
        printf("Error in cache_get4()\n");
        return result; //exit(0);
    }

    return result;
}

void mtbdd_put_result_in_cache_4(int function, uint64_t num1, uint64_t num2, uint64_t num3, uint64_t num4, uint64_t result)
{
    if(!cache_put4(function, num1, num2, num3, num4, result)) {
        printf("Error in cache_put4()\n");
        return; //exit(0);
    }
    return;
}

/**
 * Renumber vars in a decision diagram to let it start from new_var.
 * 
 * So, the var_set = {2, 3, ..., n} to a decision diagram 
 * 
 * is renumbered to {0, 1, ..., n-2} if new_var = 0.
 */

MTBDD mtbdd_renumber_variables(MTBDD M, uint32_t new_var)
{
    MTBDD result = MTBDD_ZERO;

    if(M == MTBDD_ZERO)
        return result;

    if(mtbdd_isleaf(M))
        return M;

    if(mtbdd_isnode(M) && mtbdd_getvar(M) == new_var)
        return M;

    MTBDD getlow = mtbdd_getlow(M);
    MTBDD gethigh = mtbdd_gethigh(M);

    result = mtbdd_is_result_in_cache_3(CACHE_MTBDD_RENUMBER_VARS, new_var, getlow, gethigh);
    if(result != MTBDD_ZERO)
        return result;

    MTBDD low = mtbdd_renumber_variables(getlow, new_var + 1);
    MTBDD high = mtbdd_renumber_variables(gethigh, new_var + 1);

    result = mtbdd_makenode(new_var, low, high);

    //printf("I was here new var = %d, low = %ld, high = %ld!\n", new_var, low, high);

    mtbdd_put_result_in_cache_3(CACHE_MTBDD_RENUMBER_VARS, new_var, low, high, result);

    return result;
}

void determine_top_var_and_leafcount(MTBDD M, int *botvar, int *topvar, int *leafcount)
{
    //printf("M = %ld, isleaf = %d, getvar=%d\n, ", M, mtbdd_isleaf(M), mtbdd_getvar(M));

    if(M == MTBDD_ZERO)
        return;

    if(mtbdd_isleaf(M)) {
        (*leafcount)++;
        return;
    }

    if(mtbdd_isnode(M)) {
        int var = mtbdd_getvar(M);

        if(*topvar < var)
            *topvar = var;

        if(*botvar > var)
            *botvar = var;

        determine_top_var_and_leafcount(mtbdd_getlow(M), botvar, topvar, leafcount);
        determine_top_var_and_leafcount(mtbdd_gethigh(M), botvar, topvar, leafcount);

        return;
    }

    assert(false); // M should be leaf or node

    return;
}

//
//  Get the children (low and high) on the level of the given var
//
//  var = 1: suppose [x1] is a virtual node
//
//    |         |
//    |        [x1]
//    |    -->  | | <-- return these as M_low and M_high
//    x3        x3 -
//    |         |
//
//  var = 5: create virtual node [x5] and return low and high
//
//    |         |
//    |  -->   [x5] 
//    |         | | <-- return these as M_low and M_high
//    x7        x7 -
//    |         |
//
//
//  var = 4: -> NOTE: This may not happen in the call to this function.
//
//    |                     |
//    x3 -                  x3 -
//    |                     |
//    x4 -                  x4 - M_high
//    |         M_low  -->  | 
//    x7 -                  x7 -
//    |                     |
//
//  var = 9: -> NOTE: This may not happen in the call to this function
//
//    |         |
//    x7 -  --> x7 - M_high
//    |         | <-- M_low
//   [x9] -
//    |
//

void mtbdd_get_children_of_var(MTBDD M, MTBDD *M_low, MTBDD *M_high, uint32_t var)
{
    *M_low = MTBDD_ZERO;
    *M_high = MTBDD_ZERO;

    if(M == MTBDD_ZERO) {
        *M_low = M;
        *M_high = M;
        return;
    }

    //
    // TODO: Introduce mtbddnode_t node = MTBDD_GETNODE(M); That is more efficient.
    //

    if(mtbdd_isleaf(M)) {
        *M_low = M;
        *M_high = M;
        return;
    }

    if(mtbdd_getvar(M) > var) {
        *M_low = M;
        *M_high = M;
        return;
    }

    if(mtbdd_getvar(M) == var) {
        *M_low = mtbdd_getlow(M);
        *M_high = mtbdd_gethigh(M);
        return;
    }

    if(mtbdd_getvar(M) < var) { // NOTE: This may not happen!
        printf("Error in mtbdd_get_children_of_var(): var(M) < var\n");
        // exit(0);
        return; 
    }

    return;
}

void mtbdd_split_mtbdd_into_two_parts(MTBDD v, MTBDD *v0, MTBDD *v1, uint32_t var)
{
    mtbdd_get_children_of_var(v, v0, v1, var);
    return;
}


void mtbdd_split_mtbdd_into_four_parts(MTBDD M, MTBDD *M00, MTBDD *M01, MTBDD *M10, MTBDD *M11, uint32_t var)
{
    MTBDD M_0 = MTBDD_ZERO;
    MTBDD M_1 = MTBDD_ZERO;

    mtbdd_get_children_of_var(M, &M_0, &M_1, var);
    mtbdd_get_children_of_var(M_0, M00, M01, var + 1);
    mtbdd_get_children_of_var(M_1, M10, M11, var + 1);

    return;
}

/**
 * Matrix . Vector multiplication in MTBDD domain
 * 
 * Computes M.v for an 2^n vector v and a 2^n x 2^n matrix M.
 * 
 * The vector is row wise sorted, the matrix is column wise sorted.
 * 
 * The resulting vector w is column wise sorted.
 * 
 */
MTBDD mtbdd_matvec_mult_alt(MTBDD M, MTBDD v, int n)
{
    // 
    // M(r0r1) is a matrix sorted row wise from the top
    //
    // v(c1) is a vector sorted column wise from the top
    //

    MTBDD result = MTBDD_ZERO;

    // Validate input arguments
    if(n < 0 || M == MTBDD_ZERO || v == MTBDD_ZERO)
        return result;

    // Check if result already in cache
    result = mtbdd_is_result_in_cache_3(CACHE_MTBDD_MATVEC_MULT, M, v, n);
    if(result != MTBDD_ZERO)
        return result;

    // Calculate M with size 1 x 1, and v with size 1
    if(n == 0 && mtbdd_isleaf(M) && mtbdd_isleaf(v)) {

        //
        // M should contain only one leaf, no variables
        //
        // v should contain only one leaf, no variables
        //

        return mtbdd_times(M, v); // TODO: take mtbdd_times(M,v), is type independent!
    }

    // Calculate M with size 2^1 x 2^1 and v with size 2^1
    if(n == 1) {

        // 
        // M should be row-wise top to bottom with variables 0 and 1
        //
        // v should be column-first wise with variable 1
        //
        // result is two leafs with variable 0
        //
        //  M = (a b)   v = (A)     Mv = (aA + bC)
        //      (c d)       (C)          (cA + dC)
        //
        //
        //  M =        x0
        //          x1     x1
        //         a  b   c  d
        //
        //  v =        x1
        //          A       C
        //
        //  Mv =       x0
        //      aA + bC   cA + dC
        //
        //  M.v = (x0|.x1|.a + x0|.x1.b + x0.x1|.c + x0.x1.d).(x1|.A + x1.C)
        //      =  x0|.x1|.A.a + x0.x1|.c.A + x0|.x1.b.C + x0.x1.d.C
        //      =  x0|.(A.a + b.C) + x0.(c.A + d.C)
        //

        // Check if var(M) = {0,1} and var(v) = {1} // TODO: extend for sparse mtbdd's
        if(mtbdd_getvar(M) != 0 || mtbdd_getvar(mtbdd_getlow(M)) != 1 || mtbdd_getvar(v) != 1)
            return result;

        // Prepare variable set to be removed from the v
        size_t length_var_set = 1;
        uint32_t var[length_var_set];
        var[0] = n;
        MTBDD var_set = mtbdd_set_from_array(var, length_var_set);

        // Compute multiplications and additions
        result = mtbdd_and_abstract_plus(M, v, var_set);

        // Put result in cache
        mtbdd_put_result_in_cache_3(CACHE_MTBDD_MATVEC_MULT, M, v, n, result);

        return result;
    }

/*
    // Split matrix M in four parts with size 2^(n-1) x 2^(n-1)

    MTBDD M00 = mtbdd_getlow(mtbdd_getlow(M)); // TODO: incorporate vars!
    MTBDD M10 = mtbdd_gethigh(mtbdd_getlow(M));
    MTBDD M01 = mtbdd_getlow(mtbdd_gethigh(M));
    MTBDD M11 = mtbdd_gethigh(mtbdd_gethigh(M));

    // Split vector v in two parts with size 2^(n-1)
    MTBDD v0 = mtbdd_getlow(v);
    MTBDD v1 = mtbdd_gethigh(v);

    // w0 = M00 . v0 + M10 . v1
    MTBDD w0 = mtbdd_plus(mtbdd_matvec_mult(M00, v0, n-1), mtbdd_matvec_mult(M10, v1, n-1));

    // w1 = M01 . v0 + M11 . v1
    MTBDD w1 = mtbdd_plus(mtbdd_matvec_mult(M01, v0, n-1), mtbdd_matvec_mult(M11, v1, n-1));

    result = mtbdd_makenode(0, w0, w1);

    // Put result in cache
    mtbdd_put_result_in_cache(CACHE_MTBDD_MATVEC_MULT, M, v, n, result);
*/
    return result;
}

/*
 * Matrix . Matrix multiplication in MTBDD domain
 * 
 *   Computes A.B for two 2^n x 2^n matrices.
 * 
 */

MTBDD mtbdd_matmat_mult_alt(MTBDD M1, MTBDD M2, int n)
{
    MTBDD result = MTBDD_ZERO;

    // Validate input arguments
    if(n <= 0 || M1 == MTBDD_ZERO || M2 == MTBDD_ZERO)
        return result;

    // Check if result already in cache
    result = mtbdd_is_result_in_cache_3(CACHE_MTBDD_MATMAT_MULT, M1, M2, n);
    if(result != MTBDD_ZERO) {
        printf("Result was cached (M1, M2, n) = (%ld, %ld, %d) !\n", M1, M2, n);
        return result;
    }

    // Calculate M1 with size 2^1 x 2^1 and M2 with identical size
    if(n == 1) {

        // 
        // M1 should be row-wise top to bottom with variables 0 and 1
        //
        // M2 should be column-first wise with variables 0 and 1
        //
        // result is two leafs with variable 0
        //
        //  M1 = (a b)   M2 = (A B)   M1.M2 = (a.A + b.C  a.B + b.D) = (M00 M01)
        //       (c d)        (C D)           (c.A + d.C  c.B + d.D)   (M10 M11)
        //
        //
        //  M1 =       x0           M1_ =        x0
        //          x1     x1                 x1     x1
        //         a  b   c  d               c  d   a  b
        //
        //  M2 =       x0
        //          x1     x1
        //         A  C   B  D
        //
        //  M1.M2 =     x0
        //          x1        x1
        //       M00  M10  M01  M11
        //
        //  Calculation is done in two steps:
        //
        //  1.  M00M11(x0) = M1 (x0x1).M2(x0x1), x1 reduced
        //  2.  M10M01(x0) = M1'(x0x1).M2(x0x1), x1 reduced 
        //  3.  Combine M00M11(x0) and M10M01(x0) to get M1.M2
        //
        //
        //  Proof:
        //
        //  M1.M2  = (x0|.x1|.a + x0|.x1.b + x0.x1|.c + x0.x1.d) . (x0|.x1|.A + x0|.x1.C + x0.x1|.B + x0.x1.D)
        //         => x0|.(x1|+x1).(a.A + b.C) + x0.(x1|+x1).(c.B + d.D), x1 reduced
        //         =  M00M11(x0) -> M00 = M00M11(0), M11 = M00M11(1)
        //
        //  M1'.M2 = (x0|.x1|.c + x0|.x1.d + x0.x1|.a + x0.x1.b) . (x0|.x1|.A + x0|.x1.C + x0.x1|.B + x0.x1.D)
        //         => x0|.(c.A + d.C) + x0.(a.B + b.D), x1 reduced
        //         =  M10M01(x0) -> M10 = M00M11(0), M01 = M00M11(1)
        //

        M1 = mtbdd_renumber_variables(M1, 0);
        M2 = mtbdd_renumber_variables(M2, 0);

        // Compose M1' out of M1 by swapping low and high of x0
        MTBDD M1_ = mtbdd_makenode(0, mtbdd_gethigh(M1), mtbdd_getlow(M1));

        // Prepare variable set to be removed x1 from the matrices 
        size_t length_var_set = 1;
        uint32_t var[length_var_set];
        var[0] = n;
        MTBDD var_set = mtbdd_set_from_array(var, length_var_set);

        // Compute multiplications and additions
        MTBDD M00M11 = mtbdd_and_abstract_plus(M1,  M2, var_set);
        MTBDD M10M01 = mtbdd_and_abstract_plus(M1_, M2, var_set);

        MTBDD low  = mtbdd_makenode(1, mtbdd_getlow(M00M11), mtbdd_getlow(M10M01));
        MTBDD high = mtbdd_makenode(1, mtbdd_gethigh(M10M01), mtbdd_gethigh(M00M11));

        result = mtbdd_makenode(0, low, high);

        // Put result in cache
        mtbdd_put_result_in_cache_3(CACHE_MTBDD_MATMAT_MULT, M1, M2, n, result);

        return result;
    }

    //
    // Split matrix M1 in four parts with size 2^(n-1) x 2^(n-1)
    //
    //  M1 is row-wise sorted
    //
    MTBDD M1_00 = mtbdd_getlow(  mtbdd_getlow( M1) );
    MTBDD M1_01 = mtbdd_gethigh( mtbdd_getlow( M1) );
    MTBDD M1_10 = mtbdd_getlow(  mtbdd_gethigh(M1) );
    MTBDD M1_11 = mtbdd_gethigh( mtbdd_gethigh(M1) );

    M1_00 = mtbdd_renumber_variables(M1_00, 0);
    M1_01 = mtbdd_renumber_variables(M1_01, 0);
    M1_10 = mtbdd_renumber_variables(M1_10, 0);
    M1_11 = mtbdd_renumber_variables(M1_11, 0);

    //
    // Split matrix M2 in four parts with size 2^(n-1) x 2^(n-1)
    //
    //  M2 is column-wise sorted
    //
    MTBDD M2_00 = mtbdd_getlow(  mtbdd_getlow( M2) );
    MTBDD M2_10 = mtbdd_gethigh( mtbdd_getlow( M2) );
    MTBDD M2_01 = mtbdd_getlow(  mtbdd_gethigh(M2) );
    MTBDD M2_11 = mtbdd_gethigh( mtbdd_gethigh(M2) );

    M2_00 = mtbdd_renumber_variables(M2_00, 0);
    M2_01 = mtbdd_renumber_variables(M2_01, 0);
    M2_10 = mtbdd_renumber_variables(M2_10, 0);
    M2_11 = mtbdd_renumber_variables(M2_11, 0);

    //
    // W00 = M1_00 . M2_00 + M1_01 . M2_01    W01 = M1_00 . M2_10 + M1_01 . M2_11
    //
    //  W00 and W01 are column-wise sorted
    //
    MTBDD W00 = mtbdd_plus(mtbdd_matmat_mult_alt(M1_00, M2_00, n-1), mtbdd_matmat_mult_alt(M1_01, M2_01, n-1));
    MTBDD W01 = mtbdd_plus(mtbdd_matmat_mult_alt(M1_00, M2_10, n-1), mtbdd_matmat_mult_alt(M1_01, M2_11, n-1));

    //
    // W10 = M1_10 . M2_00 + M1_11 . M2_01    W11 = M1_10 . M2_10 + M1_11 . M2_11 
    //
    //  W10 and W11 are column-wise sorted
    //
    MTBDD W10 = mtbdd_plus(mtbdd_matmat_mult_alt(M1_10, M2_00, n-1), mtbdd_matmat_mult_alt(M1_11, M2_01, n-1));
    MTBDD W11 = mtbdd_plus(mtbdd_matmat_mult_alt(M1_10, M2_10, n-1), mtbdd_matmat_mult_alt(M1_11, M2_11, n-1));

    // 
    // Compose final matrix as result
    //
    //  Matrix is column-wise sorted
    //
    MTBDD x0_low  = mtbdd_makenode(n-1, W00, W10);
    MTBDD x0_high = mtbdd_makenode(n-1, W01, W11);

    result = mtbdd_makenode(n-2, x0_low, x0_high);

    result = mtbdd_renumber_variables(result, 0);
 
    // Put result in cache
    mtbdd_put_result_in_cache_3(CACHE_MTBDD_MATMAT_MULT, M1, M2, n-2, result);

    return result;
}

/*
 * Matrix . Vector multiplication in MTBDD domain
 *
 *   Computes A.v for one 2^n matrix and one 2^n vector.
 *
 */
MTBDD mtbdd_matvec_mult(MTBDD M, MTBDD v, int nvars, int currentvar)
{
    int maxvar = -1;
    int minvar = 100;
    int leafcount = 0;
    determine_top_var_and_leafcount(M, &minvar, &maxvar, &leafcount);

    maxvar = -1;
    minvar = 100;
    leafcount = 0;
    determine_top_var_and_leafcount(v, &minvar, &maxvar, &leafcount);

    MTBDD result = MTBDD_ZERO;

    // Validate input arguments
    if(M == MTBDD_ZERO || v == MTBDD_ZERO)
        return result;

    // Multiply two leafs
    if(nvars == currentvar) {

        if(!mtbdd_isleaf(M)) {
            printf("Error in mtbdd_matvec_mult, M not a leaf\n");
            //exit(0);
        }

        if(!mtbdd_isleaf(v)) {
            printf("Error in mtbdd_matvec_mult, v not a leaf\n");
            //exit(0);
        }

        return mtbdd_times(M, v);
    }

    // Check if result already in cache
    //result = mtbdd_is_result_in_cache_4(CACHE_MTBDD_MATVEC_MULT, M, v, currentvar, nvars);
    //if(result != MTBDD_ZERO) {
    //    printf("Result was cached based on (M, v, currentvar, nvars) = (%ld, %ld, %d, %d) !\n", M, v, currentvar, nvars);
    //    return result;
    //}

    // Multiply recursive, reduce M1 and M2 with two vars
    MTBDD M_00 = M;
    MTBDD M_01 = M;
    MTBDD M_10 = M;
    MTBDD M_11 = M;

    mtbdd_split_mtbdd_into_four_parts(M, &M_00, &M_01, &M_10, &M_11, currentvar);

    maxvar = -1; minvar = 100; leafcount = 0;
    determine_top_var_and_leafcount(M, &minvar, &maxvar, &leafcount);
    printf("M: maxvar = %d, minvar = %d, leaves = %d\n", maxvar, minvar, leafcount);

    MTBDD v_0 = v;
    MTBDD v_1 = v;

    mtbdd_split_mtbdd_into_two_parts(v, &v_0, &v_1, currentvar);

    currentvar = currentvar+2;

    //
    // W00 = M_00 . v_0 + M_01 . v_1    
    //
    // W01 = M_10 . v_0 + M_11 . v_1
    //
    MTBDD w_0 = mtbdd_plus(
        mtbdd_matvec_mult(M_00, v_0, nvars, currentvar),
        mtbdd_matvec_mult(M_01, v_1, nvars, currentvar));

    MTBDD w_1 = mtbdd_plus(
        mtbdd_matvec_mult(M_10, v_0, nvars, currentvar),
        mtbdd_matvec_mult(M_11, v_1, nvars, currentvar));

    currentvar = currentvar-2;

    result = mtbdd_makenode(currentvar, w_0, w_1);

    // Put result in cache
    //mtbdd_put_result_in_cache_4(CACHE_MTBDD_MATVEC_MULT, M, v, currentvar, nvars, result);

    return result;
}

/*
 * Matrix . Matrix multiplication in MTBDD domain
 * 
 *   Computes A.B for two 2^n x 2^n matrices.
 *  
 */
MTBDD mtbdd_matmat_mult(MTBDD M1, MTBDD M2, int nvars, int currentvar)
{
/*
    printf("M1 = %ld, M2 = %ld, nvars = %d, currentvar = %d\n", M1, M2, nvars, currentvar);

    int maxvar = -1;
    int minvar = 100;
    int leafcount = 0;
    determine_top_var_and_leafcount(M1, &minvar, &maxvar, &leafcount);

    printf("M1: maxvar = %d, minvar = %d, leaves = %d\n", maxvar, minvar, leafcount);

    maxvar = -1;
    minvar = 100;
    leafcount = 0;
    determine_top_var_and_leafcount(M2, &minvar, &maxvar, &leafcount);

    printf("M2: maxvar = %d, minvar = %d, leaves = %d\n", maxvar, minvar, leafcount);
*/

    MTBDD result = MTBDD_ZERO;

    // Validate input arguments
    if(M1 == MTBDD_ZERO || M2 == MTBDD_ZERO)
        return result;

    // Multiply two leafs
    if(nvars == currentvar) {

        if(!mtbdd_isleaf(M1)) {
            printf("Error in mtbdd_matmat_mult, M1 not a leaf\n");
            //exit(0);
        }

        if(!mtbdd_isleaf(M2)) {
            printf("Error in mtbdd_matmat_mult, M2 not a leaf\n");
            //exit(0);
        }

        return mtbdd_times(M1, M2);
    }

    // Check if result already in cache
    result = mtbdd_is_result_in_cache_4(CACHE_MTBDD_MATMAT_MULT, M1, M2, currentvar, nvars);
    if(result != MTBDD_ZERO) {
        printf("Result was cached based on (M1, M2, currentvar, nvars) = (%ld, %ld, %d, %d) !\n", M1, M2, currentvar, nvars);
        return result;
    }

    // Multiply recursive, reduce M1 and M2 with two vars
    MTBDD M1_00 = M1;
    MTBDD M1_01 = M1;
    MTBDD M1_10 = M1;
    MTBDD M1_11 = M1;

    mtbdd_split_mtbdd_into_four_parts(M1, &M1_00, &M1_01, &M1_10, &M1_11, currentvar);
/*
    maxvar = -1; minvar = 100; leafcount = 0;
    determine_top_var_and_leafcount(M1_00, &minvar, &maxvar, &leafcount);
    printf("M1_00: maxvar = %d, minvar = %d, leaves = %d\n", maxvar, minvar, leafcount);
*/
    MTBDD M2_00 = M2;
    MTBDD M2_01 = M2;
    MTBDD M2_10 = M2;
    MTBDD M2_11 = M2;

    mtbdd_split_mtbdd_into_four_parts(M2, &M2_00, &M2_01, &M2_10, &M2_11, currentvar);

    currentvar = currentvar+2;
/*
    MTBDD a = mtbdd_matmat_mult(M1_00, M2_00, nvars, currentvar);
    MTBDD b = mtbdd_matmat_mult(M1_01, M2_10, nvars, currentvar);

    maxvar = -1; minvar = 100; leafcount = 0;
    determine_top_var_and_leafcount(a, &minvar, &maxvar, &leafcount);
    printf("a: maxvar = %d, minvar = %d, leaves = %d\n", maxvar, minvar, leafcount);

    maxvar = -1; minvar = 100; leafcount = 0;
    determine_top_var_and_leafcount(b, &minvar, &maxvar, &leafcount);
    printf("b: maxvar = %d, minvar = %d, leaves = %d\n", maxvar, minvar, leafcount);
*/
    MTBDD W00 = mtbdd_plus(
        mtbdd_matmat_mult(M1_00, M2_00, nvars, currentvar), 
        mtbdd_matmat_mult(M1_01, M2_10, nvars, currentvar));

    MTBDD W01 = mtbdd_plus(
        mtbdd_matmat_mult(M1_00, M2_01, nvars, currentvar), 
        mtbdd_matmat_mult(M1_01, M2_11, nvars, currentvar));

    MTBDD W10 = mtbdd_plus(
        mtbdd_matmat_mult(M1_10, M2_00, nvars, currentvar), 
        mtbdd_matmat_mult(M1_11, M2_10, nvars, currentvar));

    MTBDD W11 = mtbdd_plus(
        mtbdd_matmat_mult(M1_10, M2_01, nvars, currentvar), 
        mtbdd_matmat_mult(M1_11, M2_11, nvars, currentvar));

    currentvar = currentvar-2;

    MTBDD low  = mtbdd_makenode(currentvar+1, W00, W01);
    MTBDD high = mtbdd_makenode(currentvar+1, W10, W11);

    result = mtbdd_makenode(currentvar, low, high);

    // Put result in cache
    mtbdd_put_result_in_cache_4(CACHE_MTBDD_MATMAT_MULT, M1, M2, currentvar, nvars, result);

    return result;
}

/**
 * Multiplication of a and b matrices to get M = a . b, size 2^n x 2^n = 2^2n
 * 
 *  M(r,c) is element in row r and column c. Same for a(r,c) and b(r,c).
 * 
 * In indices:
 * 
 *  M(i,j,N) = sum k=0 to N-1 of ( a(i,k) . b(k,j) ) with N = 2^n = number of elements.
 * 
 *  M holds all elements M(i,j,N) with 0 <= i < N, and 0 <= j < N, N = number of elements in a = b
 * 
 * Derive recursive calculation of M(i,j,N):
 * 
 *  if (N == 0) 
 *      -> M(0,0,0) = NULL
 * 
 *  if (N == 1) 
 *      -> M(i,j,1) = sum k=0 to 0 of ( a(i,k) . b(k,j) ) = a(i,0) . b(0,j) and i,j < n
 * 
 *  if (N == 2)
 *      -> M(i,j,2) = sum k=0 to 1 of ( a(i,k) . b(k,j) ) = a(i,0) . b(0,j) + a(i,1) . b(1,j) and i,j < N
 *                  = M(i,j,1) + a(i,1) . b(1,j) and i,j < N
 * 
 *  if (N > 2) and (0 <= i < N) and (0 <= j < N)    
 *      -> M(i,j,N) = sum k=0 to N-1 of ( a(i,k) . b(k,j) ) and i,j < N
 *                  = sum k=0 to N-2 of ( a(i,k) . b(k,j) ) + a(i,N-1) . b(N-1,j) and i,j < N 
 *                  = M(i,j,N-1) + a(i,N-1) . b(N-1,j) and i,j < N
 *                  = M(i,j,N-2) + a(i,N-2) . b(N-2,j) + a(i,N-1) . b(N-1,j) and i,j < N
 *                  
 * Split in two -> (i = [0,1,2, ..., (N-1)/2-1] and i = [N-1)/2, ..., N-1])
 * 
 *  M, a and b are mtbdd's.
 * 
 * 
*/

/*

//
//
//
MTBDD mtbdd_matmat_indices(uint32_t n, MTBDD a, MTBDD b)
{
    MTBDD result = MTBDD_ZERO;

    if (n==0)
        return mtbdd_M(0,0,0,a,b);

    if (n==1)
        return mtbdd_M(0,0,1,a,b);

    // Fill result column-wise
    //
    // M:       |
    //          x0 -        "row", 2^1 = 2 leaves
    //          |
    //          x1 -        "col"
    //          |
    //          x2 -        "row"
    //          |
    //          x3 -        "col", 2^4 = 16 leaves, var = n-1
    //          |
    //
    // Run through the indices, M has size n x n
    //
    
    uint32_t pair = MTBDD_ZERO;

    for(uint32_t row=0; row<n; row++)
    {
        for(uint32_t col=0; col<n; col+=2)
        {
            value_low = mtbdd_M(row, col, n, a, b);
            value_high = mtbdd_M(row, col+1, n, a, b);

            pair = mtbdd_makenode(n-1, mtbdd_double(value_low), mtbdd_double(value_high));
        }

        M = mtbdd_makenode(n-2, pair, ?);
    }
    return result;
}



MTBDD mtbdd_M(uint32_t row, uint32_t col, uint32_t n, MTBDD a, MTBDD b)
{
    MTBDD result = MTBDD_ZERO;

    if (n==0)
        return result;

    if (n==1)
        return mtbdd_times(a, b);

    if (n==2)
        // ...
        return result;

    return result;
}
*/

/**
 * Kronecker multiplication or Tensor product
 * 
 * Vector (x) Vector product
 * 
 * Matrix (x) Matrix product
 * 
 */
//MTBDD mtbdd_vec_tensor_prod(MTBDD v, MTBDD w, int n)
//{}

//MTBDD mtbdd_mat_tensor_prod(MTBDD M1, MTBDD M2, int n)
//{}

/**
 * Tensor product (same function for matrix or vector MTBDDs)
 * TODO: this function should be a Lace TASK so that it can be parallelized.
 */
MTBDD mtbdd_tensor_prod(MTBDD a, MTBDD b, int leaf_var_a) // leaf_var_a == n
{
    mtbddnode_t na = MTBDD_GETNODE(a);
    mtbddnode_t nb = MTBDD_GETNODE(b);

    // If both are leaves, multiply their values
    if (mtbddnode_isleaf(na) && mtbddnode_isleaf(nb)) {
        return RUN(mtbdd_op_times, &a, &b);
    }

    // Check cache
    MTBDD result;
    if (cache_get3(CACHE_MTBDD_TENSOR, a, b, (uint64_t) leaf_var_a, &result)) {
        return result;
    }

    MTBDD low, high;
    BDDVAR var;
    // Recurse over A first (A \tensor B != B \tensor A)
    if (!mtbddnode_isleaf(na)) {
        low  = mtbdd_tensor_prod(mtbddnode_getlow(na), b, leaf_var_a);
        high = mtbdd_tensor_prod(mtbddnode_gethigh(na), b, leaf_var_a);
        var  = mtbddnode_getvariable(na);
    }
    else { // A is a leaf and B is not
        low  = mtbdd_tensor_prod(a, mtbddnode_getlow(nb), leaf_var_a);
        high = mtbdd_tensor_prod(a, mtbddnode_gethigh(nb), leaf_var_a);
        var  = mtbddnode_getvariable(nb) + leaf_var_a; // vars in B are offset by the depth of A
    }
    result = mtbdd_makenode(var, low, high);

    // Cache put
    cache_put3(CACHE_MTBDD_TENSOR, a, b, (uint64_t) leaf_var_a, result);

    return result;
}

