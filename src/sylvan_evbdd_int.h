
/* Do not include this file directly. Instead, include sylvan_int.h */

/**
 * Internals for EVBDDs
 */

#ifndef SYLVAN_EVBDD_INT_H
#define SYLVAN_EVBDD_INT_H

#include <sylvan_edge_weights.h>

/**
 * These golbal variables control some of the EVBDD behavior, and are set in 
 * sylvan_init_evbdd().
 * TODO: Maybe handle this in a cleaner way than with global variables?
 */
// using [wgts,ptr] [33,30] bits if set to true (default [23,40])
extern bool larger_wgt_indices;
extern int weight_norm_strat;
extern EVBDD_WGT (*normalize_weights)(EVBDD_WGT *, EVBDD_WGT *);


/*****************<Bit level manipulation of EVBDD / evbddnode_t>****************/

/**
 * When edge weight table <= 2^23 (larger_wgt_indices = false)
 * -----------------------------------------------------------------------------
 * EVBDD edge structure (64 bits)
 *       1 bit:  unused
 *      23 bits: index of edge weight in weight table (EVBDD_WGT)
 *      40 bits: index of next node in node table (EVBDD_TARG)
 * 
 * EVBDD node structure (128 bits)
 * (note: because of normalization of the edge weights, we only need 1 weight
 *  per node, the other will always be 0 or 1 or dependent on the first value.)
 * 
 * 64 bits low:
 *       1 bit:  unused
 *      16 bits: variable/qubit number of this node
 *       1 bit:  if 0 (1) normalized WGT is on low (high)
 *       1 bit:  if 0 (1) normalized WGT is EVBDD_ZERO (EVBDD_ONE)
 *       5 bits: unused
 *      40 bits: low edge pointer to next node (EVBDD_TARG)
 * 64 bits high:
 *       1 bit:  marked/unmarked flag
 *      23 bits: index of edge weight of high edge in ctable (EVBDD_WGT)
 *      40 bits: high edge pointer to next node (EVBDD_TARG)
 * -----------------------------------------------------------------------------
 * 
 * 
 * When edge weight table > 2^23 (larger_wgt_indices = true)
 * -----------------------------------------------------------------------------
 * EVBDD edge structure (64 bits)
 *       1 bit:  unused
 *      23 bits: index of edge weight in weight table (EVBDD_WGT)
 *      40 bits: index of next node in node table (EVBDD_TARG)
 * 
 * EVBDD node structure (128 bits)
 * 64 bits low:
 *       1 bit:  unused
 *      16 bits: variable/qubit number of this node
 *       1 bit:  if 0 (1) normalized WGT is on low (high)
 *       1 bit:  if 0 (1) normalized WGT is EVBDD_ZERO (EVBDD_ONE)
 *      15 bits: unused
 *      30 bits: low edge pointer to next node (EVBDD_TARG)
 * 64 bits high:
 *       1 bit:  marked/unmarked flag
 *      33 bits: index of edge weight of high edge in ctable (EVBDD_WGT)
 *      30 bits: high edge pointer to next node (EVBDD_TARG)
 * -----------------------------------------------------------------------------
 */
typedef struct __attribute__((packed)) evbddnode {
    EVBDD low, high;
} *evbddnode_t; // 16 bytes

static const EVBDD evbdd_marked_mask  = 0x8000000000000000LL; // 1000000000000000000000000000000000000000000000000000000000000000
static const EVBDD evbdd_var_mask_low = 0x7fff800000000000LL; // 0111111111111111100000000000000000000000000000000000000000000000
static const EVBDD evbdd_wgt_pos_mask = 0x0000400000000000LL; // 0000000000000000010000000000000000000000000000000000000000000000
static const EVBDD evbdd_wgt_val_mask = 0x0000200000000000LL; // 0000000000000000001000000000000000000000000000000000000000000000
static const EVBDD evbdd_wgt_mask_23  = 0x7fffff0000000000LL; // 0111111111111111111111110000000000000000000000000000000000000000
static const EVBDD evbdd_wgt_mask_33  = 0x7fffffffc0000000LL; // 0111111111111111111111111111111111000000000000000000000000000000
static const EVBDD evbdd_ptr_mask_30  = 0x000000003fffffffLL; // 0000000000000000000000000000000000111111111111111111111111111111
static const EVBDD evbdd_ptr_mask_40  = 0x000000ffffffffffLL; // 0000000000000000000000001111111111111111111111111111111111111111


/**
 * Gets only the EVBDD_WGT information of an EVBDD edge.
 */
static inline EVBDD_WGT
EVBDD_WEIGHT(EVBDD a)
{
    if (larger_wgt_indices) {
        return (a & evbdd_wgt_mask_33) >> 30; // 33 bits
    }
    else {
        return (a & evbdd_wgt_mask_23) >> 40; // 23 bits
    }
}

/**
 * Gets only the TARGET information of an EVBDD edge.
 */
static inline EVBDD_TARG
EVBDD_TARGET(EVBDD a)
{
    if (larger_wgt_indices) {
        return a & evbdd_ptr_mask_30; // 30 bits
    }
    else {
        return a & evbdd_ptr_mask_40; // 40 bits
    }
}

/**
 * Gets the variable number of a given node `n`.
 */
static inline BDDVAR
evbddnode_getvar(evbddnode_t n)
{
    return (BDDVAR) ((n->low & evbdd_var_mask_low) >> 47 ); // 16 bits
}

/**
 * Gets only the EVBDD_TARG of the low edge of <n>.
 */
static inline EVBDD_TARG
evbddnode_getptrlow(evbddnode_t n)
{
    return (EVBDD_TARG) EVBDD_TARGET(n->low);
}

/**
 * Gets only the EVBDD_TARG of the high edge of <n>.
 */
static inline EVBDD_TARG
evbddnode_getptrhigh(evbddnode_t n)
{
    return (EVBDD_TARG) EVBDD_TARGET(n->high);
}

/**
 * Gets the value of the "marked" flag.
 */
static inline bool
evbddnode_getmark(evbddnode_t n)
{
    return n->high & evbdd_marked_mask ? 1 : 0;
}

/**
 * Sets the value of the "marked" flag to `mark`.
 */
static inline void
evbddnode_setmark(evbddnode_t n, bool mark)
{
    if (mark) n->high |=  evbdd_marked_mask; // set 1st bit from left to 1
    else      n->high &= ~evbdd_marked_mask; // set 1st bit from left to 0
}

/**
 * Gets the node <p> is pointing to.
 * TODO (?) return special node for when p == EVBDD_TERMINAL
 */
static inline evbddnode_t
EVBDD_GETNODE(EVBDD_TARG p)
{
    return (evbddnode_t) llmsset_index_to_ptr(nodes, p);
}

/**
 * Packs a EVBDD_TARG and EVBDD_WGT into a single 64 bit EVBDD.
 */
static inline EVBDD
evbdd_bundle(EVBDD_TARG p, EVBDD_WGT a)
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
evbddnode_unpack(evbddnode_t n, EVBDD_TARG *low, EVBDD_TARG *high, EVBDD_WGT *a, EVBDD_WGT *b)
{
    *low  = evbddnode_getptrlow(n);
    *high = evbddnode_getptrhigh(n);
    bool norm_pos = (n->low & evbdd_wgt_pos_mask) >> 46;
    bool norm_val = (n->low & evbdd_wgt_val_mask) >> 45;

    if (weight_norm_strat == NORM_L2) {
        *b = EVBDD_WEIGHT(n->high);
        *a = wgt_get_low_L2normed(*b);
    }
    else {
        if (norm_pos == 0) { // low WGT is EVBDD_ZERO or EVBDD_ONE, high WGT in table
            *a = (norm_val == 0) ? EVBDD_ZERO : EVBDD_ONE;
            *b = EVBDD_WEIGHT(n->high);
        }
        else { // high WGT is EVBDD_ZERO or EVBDD_ONE, low WGT in table
            *b = (norm_val == 0) ? EVBDD_ZERO : EVBDD_ONE;
            *a = EVBDD_WEIGHT(n->high);
        }
    }
}

static void __attribute__((unused))
evbddnode_getchilderen(evbddnode_t n, EVBDD *low, EVBDD *high)
{
    EVBDD_TARG l, h;
    EVBDD_WGT  a, b;
    evbddnode_unpack(n, &l, &h, &a, &b);
    *low  = evbdd_bundle(l, a);
    *high = evbdd_bundle(h, b);
}

static void __attribute__((unused))
evbddnode_pack(evbddnode_t n, BDDVAR var, EVBDD_TARG low, EVBDD_TARG high, EVBDD_WGT a, EVBDD_WGT b)
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
    EVBDD_WGT wgt_high;
    bool norm_pos;
    bool norm_val;
    
    if (weight_norm_strat == NORM_L2) {
        assert(!(a == EVBDD_ZERO && b == EVBDD_ZERO)); // redundant node (caught before)
        norm_pos = 0;
        norm_val = 0;
        wgt_high = b; // we can derive a from b
    }
    else {
        /// weight_norm_strat == NORM_LOW or NORM_MAX or NORM_MIN
        assert(a == EVBDD_ZERO || a == EVBDD_ONE || b == EVBDD_ZERO || b == EVBDD_ONE);
        norm_pos = (a == EVBDD_ZERO || a == EVBDD_ONE) ? 0 : 1;
        if (norm_pos == 0) {
            norm_val = (a == EVBDD_ZERO) ? 0 : 1;
            wgt_high = b;
        }
        else {
            norm_val = (b == EVBDD_ZERO) ? 0 : 1;
            wgt_high = a;
        }
    }

    // organize the bit structure of low and high
    n->low  = ((uint64_t)var)<<47 | ((uint64_t)norm_pos)<<46 | ((uint64_t)norm_val)<<45 | low;
    if (larger_wgt_indices) {
        n->high = wgt_high<<30 | high;
    }
    else {
        n->high = wgt_high<<40 | high;
    }
}

static EVBDD_TARG __attribute__((unused))
_evbdd_makenode(BDDVAR var, EVBDD_TARG low, EVBDD_TARG high, EVBDD_WGT a, EVBDD_WGT b)
{
    struct evbddnode n;

    evbddnode_pack(&n, var, low, high, a, b);

    EVBDD_TARG result;
    int created;
    EVBDD_TARG index = llmsset_lookup(nodes, n.low, n.high, &created);
    if (index == 0) {
        //printf("auto gc of node table triggered\n");

        evbdd_refs_push(low);
        evbdd_refs_push(high);
        sylvan_gc();
        evbdd_refs_pop(2);

        index = llmsset_lookup(nodes, n.low, n.high, &created);
        if (index == 0) {
            fprintf(stderr, "EVBDD/BDD Unique table full, %zu of %zu buckets filled!\n", llmsset_count_marked(nodes), llmsset_get_size(nodes));
            exit(1);
        }
    }

    if (created) sylvan_stats_count(EVBDD_NODES_CREATED);
    else sylvan_stats_count(EVBDD_NODES_REUSED);

    result = index;
    //return mark ? result | evbdd_marked_mask : result;
    return result;
}

static EVBDD __attribute__((unused))
evbdd_makenode(BDDVAR var, EVBDD low, EVBDD high)
{ 
    if (var > 1<<16) {
        fprintf(stderr, "ERROR: EVBDDs currently only support up to %d variables.\n", 1<<16);
        exit(EXIT_FAILURE);
    }
    EVBDD_TARG low_trg  = EVBDD_TARGET(low);
    EVBDD_WGT  low_wgt  = EVBDD_WEIGHT(low);
    EVBDD_TARG high_trg = EVBDD_TARGET(high);
    EVBDD_WGT  high_wgt = EVBDD_WEIGHT(high);

    // Edges with weight 0 should point straight to terminal.
    if (low_wgt  == EVBDD_ZERO) low_trg  = EVBDD_TERMINAL;
    if (high_wgt == EVBDD_ZERO) high_trg = EVBDD_TERMINAL;

    // If both low and high are the same (both TARG and WGT) return low
    if (low == high) return low;
    else {
        // If the edges are not the same
        EVBDD_WGT norm  = (*normalize_weights)(&low_wgt, &high_wgt);
        EVBDD_TARG res  = _evbdd_makenode(var, low_trg, high_trg, low_wgt, high_wgt);
        return evbdd_bundle(res, norm);
    }
}

/****************</Bit level manipulation of EVBDD / evbddnode_t>****************/


/**
 * Gets either the top node of the EVBDD, or the node with var 't' if this
 * variable would otherwise be skipped. The "node" is returned as 
 * (*topvar, *low, *high).
 */
static void __attribute__((unused))
evbdd_get_topvar(EVBDD a, BDDVAR t, BDDVAR *topvar, EVBDD *low, EVBDD *high)
{
    bool skipped = false;
    if(EVBDD_TARGET(a) == EVBDD_TERMINAL) {
        skipped = true;
    }
    else {
        evbddnode_t node = EVBDD_GETNODE(EVBDD_TARGET(a));
        *topvar = evbddnode_getvar(node);
        if (*topvar > t) skipped = true;
    }

    if (skipped) {
        *low  = evbdd_bundle(EVBDD_TARGET(a), EVBDD_ONE);
        *high = evbdd_bundle(EVBDD_TARGET(a), EVBDD_ONE);
        *topvar = t;
    }
    else {
        evbddnode_t node = EVBDD_GETNODE(EVBDD_TARGET(a));
        evbddnode_getchilderen(node, low, high);
    }
}

#endif