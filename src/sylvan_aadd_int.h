
/* Do not include this file directly. Instead, include sylvan_int.h */

/**
 * Internals for AADDs
 */

#ifndef SYLVAN_AADD_INT_H
#define SYLVAN_AADD_INT_H

#include <sylvan_edge_weights.h>

/**
 * These golbal variables control some of the AADD behavior, and are set in 
 * sylvan_init_aadd().
 * TODO: Maybe handle this in a cleaner way than with global variables?
 */
// using [wgts,ptr] [33,30] bits if set to true (default [23,40])
extern bool larger_wgt_indices;
extern int weight_norm_strat;
extern AADD_WGT (*normalize_weights)(AADD_WGT *, AADD_WGT *);


/*****************<Bit level manipulation of AADD / aaddnode_t>****************/

/**
 * When edge weight table <= 2^23 (larger_wgt_indices = false)
 * -----------------------------------------------------------------------------
 * AADD edge structure (64 bits)
 *       1 bit:  unused
 *      23 bits: index of edge weight in weight table (AADD_WGT)
 *      40 bits: index of next node in node table (AADD_TARG)
 * 
 * AADD node structure (128 bits)
 * (note: because of normalization of the edge weights, we only need 1 weight
 *  per node, the other will always be 0 or 1 or dependent on the first value.)
 * 
 * 64 bits low:
 *       1 bit:  unused
 *       8 bits: variable/qubit number of this node
 *       1 bit:  if 0 (1) normalized WGT is on low (high)
 *       1 bit:  if 0 (1) normalized WGT is AADD_ZERO (AADD_ONE)
 *      13 bits: unused
 *      40 bits: low edge pointer to next node (AADD_TARG)
 * 64 bits high:
 *       1 bit:  marked/unmarked flag
 *      23 bits: index of edge weight of high edge in ctable (AADD_WGT)
 *      40 bits: high edge pointer to next node (AADD_TARG)
 * -----------------------------------------------------------------------------
 * 
 * 
 * When edge weight table > 2^23 (larger_wgt_indices = true)
 * -----------------------------------------------------------------------------
 * AADD edge structure (64 bits)
 *       1 bit:  unused
 *      23 bits: index of edge weight in weight table (AADD_WGT)
 *      40 bits: index of next node in node table (AADD_TARG)
 * 
 * AADD node structure (128 bits)
 * 64 bits low:
 *       1 bit:  unused
 *       8 bits: variable/qubit number of this node
 *       1 bit:  if 0 (1) normalized WGT is on low (high)
 *       1 bit:  if 0 (1) normalized WGT is AADD_ZERO (AADD_ONE)
 *      23 bits: unused
 *      30 bits: low edge pointer to next node (AADD_TARG)
 * 64 bits high:
 *       1 bit:  marked/unmarked flag
 *      33 bits: index of edge weight of high edge in ctable (AADD_WGT)
 *      30 bits: high edge pointer to next node (AADD_TARG)
 * -----------------------------------------------------------------------------
 */
typedef struct __attribute__((packed)) aaddnode {
    AADD low, high;
} *aaddnode_t; // 16 bytes

static const AADD aadd_marked_mask  = 0x8000000000000000LL;
static const AADD aadd_var_mask_low = 0x7f80000000000000LL;
static const AADD aadd_wgt_pos_mask = 0x0040000000000000LL;
static const AADD aadd_wgt_val_mask = 0x0020000000000000LL;
static const AADD aadd_wgt_mask_23  = 0x7fffff0000000000LL;
static const AADD aadd_wgt_mask_33  = 0x7fffffffc0000000LL;
static const AADD aadd_ptr_mask_30  = 0x000000003fffffffLL;
static const AADD aadd_ptr_mask_40  = 0x000000ffffffffffLL;


/**
 * Gets only the AADD_WGT information of an AADD edge.
 */
static inline AADD_WGT
AADD_WEIGHT(AADD a)
{
    if (larger_wgt_indices) {
        return (a & aadd_wgt_mask_33) >> 30; // 33 bits
    }
    else {
        return (a & aadd_wgt_mask_23) >> 40; // 23 bits
    }
}

/**
 * Gets only the TARGET information of an AADD edge.
 */
static inline AADD_TARG
AADD_TARGET(AADD a)
{
    if (larger_wgt_indices) {
        return a & aadd_ptr_mask_30; // 30 bits
    }
    else {
        return a & aadd_ptr_mask_40; // 40 bits
    }
}

/**
 * Gets the variable number of a given node `n`.
 */
static inline BDDVAR
aaddnode_getvar(aaddnode_t n)
{
    return (BDDVAR) ((n->low & aadd_var_mask_low) >> 55 ); // 8 bits
}

/**
 * Gets only the AADD_TARG of the low edge of <n>.
 */
static inline AADD_TARG
aaddnode_getptrlow(aaddnode_t n)
{
    return (AADD_TARG) AADD_TARGET(n->low);
}

/**
 * Gets only the AADD_TARG of the high edge of <n>.
 */
static inline AADD_TARG
aaddnode_getptrhigh(aaddnode_t n)
{
    return (AADD_TARG) AADD_TARGET(n->high);
}

/**
 * Gets the value of the "marked" flag.
 */
static inline bool
aaddnode_getmark(aaddnode_t n)
{
    return n->high & aadd_marked_mask ? 1 : 0;
}

/**
 * Sets the value of the "marked" flag to `mark`.
 */
static inline void
aaddnode_setmark(aaddnode_t n, bool mark)
{
    if (mark) n->high |=  aadd_marked_mask; // set 1st bit from left to 1
    else      n->high &= ~aadd_marked_mask; // set 1st bit from left to 0
}

/**
 * Gets the node <p> is pointing to.
 * TODO (?) return special node for when p == AADD_TERMINAL
 */
static inline aaddnode_t
AADD_GETNODE(AADD_TARG p)
{
    return (aaddnode_t) llmsset_index_to_ptr(nodes, p);
}

/**
 * Packs a AADD_TARG and AADD_WGT into a single 64 bit AADD.
 */
static inline AADD
aadd_bundle(AADD_TARG p, AADD_WGT a)
{
    if (larger_wgt_indices) {
        assert (p <= 0x000000003ffffffe);   // avoid clash with sylvan_invalid
        assert (a <= (1LL<<33));
        return (a << 30 | p);
    }else {
        assert (p <= 0x000000fffffffffe);   // avoid clash with sylvan_invalid
        assert (a <= (1<<23));
        return (a << 40 | p);
    }
}

static void __attribute__((unused))
aaddnode_unpack(aaddnode_t n, AADD_TARG *low, AADD_TARG *high, AADD_WGT *a, AADD_WGT *b)
{
    *low  = aaddnode_getptrlow(n);
    *high = aaddnode_getptrhigh(n);
    bool norm_pos = (n->low & aadd_wgt_pos_mask) >> 54;
    bool norm_val = (n->low & aadd_wgt_val_mask) >> 53;

    if (weight_norm_strat == NORM_L2) {
        *b = AADD_WEIGHT(n->high);
        *a = wgt_get_low_L2normed(*b);
    }
    else {
        if (norm_pos == 0) { // low WGT is AADD_ZERO or AADD_ONE, high WGT in table
            *a = (norm_val == 0) ? AADD_ZERO : AADD_ONE;
            *b = AADD_WEIGHT(n->high);
        }
        else { // high WGT is AADD_ZERO or AADD_ONE, low WGT in table
            *b = (norm_val == 0) ? AADD_ZERO : AADD_ONE;
            *a = AADD_WEIGHT(n->high);
        }
    }
}

static void __attribute__((unused))
aaddnode_getchilderen(aaddnode_t n, AADD *low, AADD *high)
{
    AADD_TARG l, h;
    AADD_WGT  a, b;
    aaddnode_unpack(n, &l, &h, &a, &b);
    *low  = aadd_bundle(l, a);
    *high = aadd_bundle(h, b);
}

static void __attribute__((unused))
aaddnode_pack(aaddnode_t n, BDDVAR var, AADD_TARG low, AADD_TARG high, AADD_WGT a, AADD_WGT b)
{
    // We only want to store 1 edge weight per node (which has 2 outgoing
    // edges). For NORM_LOW and NORM_MAX this is relatively easy because in
    // both those cases there is at least one edge weight equal to 1 or 0.
    //
    // For NORM_L2 it is a bit more complicated: both edge weights can be
    // outside of {0, 1}, but under the constraint that |low|^2 + |high|^2 = 1,
    // (or both are 0) and that |low| \in R+, we only need to store high, and
    // can derive low.

    // these will be set depending on the normalization strategy
    // (retrieval of edge weights is also dependent on normalization strategy)
    AADD_WGT wgt_high;
    bool norm_pos;
    bool norm_val;
    
    if (weight_norm_strat == NORM_L2) {
        assert(!(a == AADD_ZERO && b == AADD_ZERO)); // redundant node (caught before)
        norm_pos = 0;
        norm_val = 0;
        wgt_high = b; // we can derive a from b
    }
    else {
        /// weight_norm_strat == NORM_LOW or NORM_MAX or NORM_MIN
        assert(a == AADD_ZERO || a == AADD_ONE || b == AADD_ZERO || b == AADD_ONE);
        norm_pos = (a == AADD_ZERO || a == AADD_ONE) ? 0 : 1;
        if (norm_pos == 0) {
            norm_val = (a == AADD_ZERO) ? 0 : 1;
            wgt_high = b;
        }
        else {
            norm_val = (b == AADD_ZERO) ? 0 : 1;
            wgt_high = a;
        }
    }

    // organize the bit structure of low and high
    n->low  = ((uint64_t)var)<<55 | ((uint64_t)norm_pos)<<54 | ((uint64_t)norm_val)<<53 | low;
    if (larger_wgt_indices) {
        n->high = wgt_high<<30 | high;
    }
    else {
        n->high = wgt_high<<40 | high;
    }
}

static AADD_TARG __attribute__((unused))
_aadd_makenode(BDDVAR var, AADD_TARG low, AADD_TARG high, AADD_WGT a, AADD_WGT b)
{
    struct aaddnode n;

    aaddnode_pack(&n, var, low, high, a, b);

    AADD_TARG result;
    int created;
    AADD_TARG index = llmsset_lookup(nodes, n.low, n.high, &created);
    if (index == 0) {
        printf("auto gc of node table triggered\n");

        aadd_refs_push(low);
        aadd_refs_push(high);
        sylvan_gc();
        aadd_refs_pop(2);

        index = llmsset_lookup(nodes, n.low, n.high, &created);
        if (index == 0) {
            fprintf(stderr, "AADD/BDD Unique table full, %zu of %zu buckets filled!\n", llmsset_count_marked(nodes), llmsset_get_size(nodes));
            exit(1);
        }
    }

    if (created) sylvan_stats_count(AADD_NODES_CREATED);
    else sylvan_stats_count(AADD_NODES_REUSED);

    result = index;
    //return mark ? result | aadd_marked_mask : result;
    return result;
}

static AADD __attribute__((unused))
aadd_makenode(BDDVAR var, AADD low, AADD high)
{ 
    AADD_TARG low_trg  = AADD_TARGET(low);
    AADD_WGT  low_wgt  = AADD_WEIGHT(low);
    AADD_TARG high_trg = AADD_TARGET(high);
    AADD_WGT  high_wgt = AADD_WEIGHT(high);

    // Edges with weight 0 should point straight to terminal.
    if (low_wgt  == AADD_ZERO) low_trg  = AADD_TERMINAL;
    if (high_wgt == AADD_ZERO) high_trg = AADD_TERMINAL;

    // If both low and high are the same (both TARG and WGT) return low
    if (low == high) return low;
    else {
        // If the edges are not the same
        AADD_WGT norm  = (*normalize_weights)(&low_wgt, &high_wgt);
        AADD_TARG res  = _aadd_makenode(var, low_trg, high_trg, low_wgt, high_wgt);
        return aadd_bundle(res, norm);
    }
}

/****************</Bit level manipulation of AADD / aaddnode_t>****************/


/**
 * Gets either the top node of the AADD, or the node with var 't' if this
 * variable would otherwise be skipped. The "node" is returned as 
 * (*topvar, *low, *high).
 */
static void __attribute__((unused))
aadd_get_topvar(AADD a, BDDVAR t, BDDVAR *topvar, AADD *low, AADD *high)
{
    bool skipped = false;
    if(AADD_TARGET(a) == AADD_TERMINAL) {
        skipped = true;
    }
    else {
        aaddnode_t node = AADD_GETNODE(AADD_TARGET(a));
        *topvar = aaddnode_getvar(node);
        if (*topvar > t) skipped = true;
    }

    if (skipped) {
        *low  = aadd_bundle(AADD_TARGET(a), AADD_ONE);
        *high = aadd_bundle(AADD_TARGET(a), AADD_ONE);
        *topvar = t;
    }
    else {
        aaddnode_t node = AADD_GETNODE(AADD_TARGET(a));
        aaddnode_getchilderen(node, low, high);
    }
}

#endif