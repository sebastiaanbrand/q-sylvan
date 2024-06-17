#include <qsylvan_simulator.h>

#include <inttypes.h>


static bool testing_mode = 0; // turns on/off (expensive) sanity checks
static int granularity = 1; // operation cache access granularity


/***************<Helper functions for chaching QMDD operations>****************/

// Pack 2 BDDVARs (assumed max 8 bits each)
static inline uint32_t
QMDD_PARAM_PACK_16(BDDVAR a, BDDVAR b) 
{
    return b<<8 | a;
}

//Pack 24 bit gateid, 2 possible qubit parameters (e.g. control/target)
static inline uint64_t
GATE_OPID_40(uint32_t gateid, BDDVAR a, BDDVAR b)
{
    uint64_t res = ((uint64_t)b)<<32 | ((uint64_t)a)<<24 | gateid;
    return res;
}

// Pack 24 bit gateid w/ 5 possible qubit parameters (e.g. target/control/range)
static inline uint64_t
GATE_OPID_64(uint32_t gateid, BDDVAR a, BDDVAR b, BDDVAR c, BDDVAR d, BDDVAR e)
{
    uint64_t res = ((uint64_t)e)<<56 | 
                   ((uint64_t)d)<<48 | 
                   ((uint64_t)c)<<40 | 
                   ((uint64_t)b)<<32 | 
                   ((uint64_t)a)<<24 | 
                   gateid;
    return res;
}

/**************</Helper functions for chaching QMDD operations>****************/





/******************************<Initialization>********************************/

void
qsylvan_init_simulator(size_t min_tablesize, size_t max_tablesize, double wgt_tab_tolerance, int edge_weigth_backend, int norm_strat)
{
    sylvan_init_aadd(min_tablesize, max_tablesize, wgt_tab_tolerance, edge_weigth_backend, norm_strat, &qmdd_gates_init);
}

void
qsylvan_init_defaults(size_t wgt_tab_size)
{
    qsylvan_init_simulator(wgt_tab_size, wgt_tab_size, -1, COMP_HASHMAP, NORM_LOW);
}

/*****************************</Initialization>********************************/





/***************************<Initial state creation>***************************/

QMDD
qmdd_create_all_zero_state(BDDVAR n)
{
    bool x[n];
    for (BDDVAR k=0; k<n; k++) x[k] = 0;
    return qmdd_create_basis_state(n, x);
}

QMDD
qmdd_create_basis_state(BDDVAR n, bool* x)
{
    // start at terminal, and build backwards
    QMDD low, high, prev = AADD_TERMINAL;

    for (int k = n-1; k >= 0; k--) {
        if (x[k] == 0) {
            low = aadd_bundle(AADD_TARGET(prev), AADD_ONE);
            high = aadd_bundle(AADD_TERMINAL, AADD_ZERO);
        }
        else if (x[k] == 1) {
            low = aadd_bundle(AADD_TERMINAL, AADD_ZERO);
            high = aadd_bundle(AADD_TARGET(prev), AADD_ONE);
        }
        // add node to unique table
        prev = aadd_makenode(k, low, high);
    }
    return prev;
}

QMDD
qmdd_stack_matrix(QMDD below, BDDVAR k, gate_id_t gateid)
{
    // This function effectively does a Kronecker product gate \tensor below
    BDDVAR s, t;
    QMDD u00, u01, u10, u11, low, high, res;

    // Even + uneven variable are used to encode the 4 values
    s = 2*k;
    t = s + 1;

    // Matrix U = [u00 u01
    //             u10 u11] endoded in a small tree
    u00 = aadd_bundle(AADD_TARGET(below), gates[gateid][0]);
    u10 = aadd_bundle(AADD_TARGET(below), gates[gateid][2]);
    u01 = aadd_bundle(AADD_TARGET(below), gates[gateid][1]);
    u11 = aadd_bundle(AADD_TARGET(below), gates[gateid][3]);
    low  = aadd_makenode(t, u00, u10);
    high = aadd_makenode(t, u01, u11);
    res  = aadd_makenode(s, low, high);

    // Propagate common factor on previous root amp to new root amp
    AMP new_root_amp = wgt_mul(AADD_WEIGHT(below), AADD_WEIGHT(res));
    res = aadd_bundle(AADD_TARGET(res), new_root_amp);
    return res;
}

QMDD
qmdd_stack_control(QMDD case0, QMDD case1, BDDVAR k)
{
    // Effectively does |0><0| \tensor case0 + |1><1| \tensor case1
    BDDVAR s, t;
    QMDD u00, u01, u10, u11, low, high, res;

    s = 2*k;
    t = s + 1;

    u00 = case0;
    u10 = aadd_bundle(AADD_TERMINAL, AADD_ZERO);
    u01 = aadd_bundle(AADD_TERMINAL, AADD_ZERO);
    u11 = case1;
    low  = aadd_makenode(t, u00, u10);
    high = aadd_makenode(t, u01, u11);
    res  = aadd_makenode(s, low, high);

    // Weights of case0/case1 already dealt with by aadd_makenode
    return res;
}

QMDD
qmdd_create_all_identity_matrix(BDDVAR n)
{
    // Start at terminal and build backwards
    QMDD prev = aadd_bundle(AADD_TERMINAL, AADD_ONE);
    for (int k = n-1; k >= 0; k--) {
        prev = qmdd_stack_matrix(prev, k, GATEID_I);
    }
    return prev;
}

QMDD
qmdd_create_single_qubit_gate(BDDVAR n, BDDVAR t, gate_id_t gateid)
{
    // Start at terminal and build backwards
    QMDD prev = aadd_bundle(AADD_TERMINAL, AADD_ONE);
    for (int k = n-1; k >= 0; k--) {
        if ((unsigned int)k == t)
            prev = qmdd_stack_matrix(prev, k, gateid);
        else
            prev = qmdd_stack_matrix(prev, k, GATEID_I);
    }
    return prev;
}

QMDD
qmdd_create_single_qubit_gates(BDDVAR n, gate_id_t *gateids)
{
    // Start at terminal and build backwards
    QMDD prev = aadd_bundle(AADD_TERMINAL, AADD_ONE);
    for (int k = n-1; k >= 0; k--) {
        prev = qmdd_stack_matrix(prev, k, gateids[k]);
    }
    return prev;
}

QMDD
qmdd_create_single_qubit_gates_same(BDDVAR n, gate_id_t gateid)
{
    // Start at terminal and build backwards
    QMDD prev = aadd_bundle(AADD_TERMINAL, AADD_ONE);
    for (int k = n-1; k >= 0; k--) {
        prev = qmdd_stack_matrix(prev, k, gateid);
    }
    return prev;
}

QMDD
qmdd_create_controlled_gate(BDDVAR n, BDDVAR c, BDDVAR t, gate_id_t gateid)
{
    // for now, assume t > c
    assert(t > c);
    // Start at terminal and build backwards
    QMDD prev = aadd_bundle(AADD_TERMINAL, AADD_ONE);
    QMDD branch0 = AADD_TERMINAL, branch1 = AADD_TERMINAL;
    for (int k = n-1; k>= 0; k--) {
        if ((unsigned int)k > t || (unsigned int) k < c) {
            prev = qmdd_stack_matrix(prev, k, GATEID_I);
        }
        else if ((unsigned int) k == t) {
            branch0 = qmdd_stack_matrix(prev, k, GATEID_I);
            branch1 = qmdd_stack_matrix(prev, k, gateid);
        }
        else if ((unsigned int) k < t && (unsigned int) k > c) {
            branch0 = qmdd_stack_matrix(branch0, k, GATEID_I);
            branch1 = qmdd_stack_matrix(branch1, k, GATEID_I);
        }
        else if ((unsigned int) k == c) {
            prev = qmdd_stack_control(branch0, branch1, k);
        }
        else {
            assert("all cases should have been covered" && false);
        }
    }
    return prev;
}

QMDD
qmdd_create_multi_cgate(BDDVAR n, int *c_options, gate_id_t gateid)
{
    // c_options[k] = -1 -> ignore qubit k (apply I)
    // c_options[k] =  0 -> control on q_k = |0>
    // c_options[k] =  1 -> control on q_k = |1>
    // c_options[k] =  2 -> target qubit (for now assume 1 target)

    // C_1 C_2 C_3 U_0 =  U \tensor |111><111| +    (U_proj)
    //                    I^{\tensor 4} +           (I)
    //                   -I \tensor |111><111|      (proj)

    // TODO: protect from gc?
    QMDD U_proj = aadd_bundle(AADD_TERMINAL, AADD_ONE);
    QMDD proj   = aadd_bundle(AADD_TERMINAL, AADD_ONE);
    QMDD I      = qmdd_create_all_identity_matrix(n);

    for (int k = n-1; k >= 0; k--) {
        // -1 : Ignore qubit (apply I)
        if (c_options[k] == -1) {
            U_proj = qmdd_stack_matrix(U_proj, k, GATEID_I);
            proj   = qmdd_stack_matrix(proj,   k, GATEID_I);
        }
        // 0 : control on q_k = |0>
        else if (c_options[k] == 0) {
            U_proj = qmdd_stack_matrix(U_proj, k, GATEID_proj0);
            proj   = qmdd_stack_matrix(proj,   k, GATEID_proj0);
        }
        // 1 : control on q_k = |1>
        else if (c_options[k] == 1) {
            U_proj = qmdd_stack_matrix(U_proj, k, GATEID_proj1);
            proj   = qmdd_stack_matrix(proj,   k, GATEID_proj1);
        }
        // 2 : target qubit
        else if (c_options[k] == 2) {
            U_proj = qmdd_stack_matrix(U_proj, k, gateid);
            proj   = qmdd_stack_matrix(proj,   k, GATEID_I);
        }
        else {
            printf("Invalid option %d for qubit %d (options = {-1,0,1,2}\n", c_options[k], k);
            exit(1);
        }
    }

    // multiply root edge of proj with -1
    proj = aadd_bundle(AADD_TARGET(proj), wgt_neg(AADD_WEIGHT(proj)));

    return aadd_plus(U_proj, aadd_plus(I, proj));
}

QMDD
qmdd_create_all_control_phase(BDDVAR n, bool *x)
{
    QMDD identity = aadd_bundle(AADD_TERMINAL, AADD_ONE);
    QMDD ccphase  = aadd_bundle(AADD_TERMINAL, AADD_ONE);

    // Start with (-1)Z gate on last qubit. Z if control on 1 and -Z if 0.
    if (x[n-1] == 1) {
        ccphase = qmdd_stack_matrix(ccphase, n-1, GATEID_Z);
    }
    else if (x[n-1] == 0) {
        ccphase = qmdd_stack_matrix(ccphase, n-1, GATEID_Z);
        ccphase = aadd_bundle(AADD_TARGET(ccphase), wgt_neg(AADD_WEIGHT(ccphase)));
    }

    // Stack remaining controls
    for (int k = n-2; k >= 0; k--) {
        // "Identity stack" for doing nothing on each qubit's non-control branch
        identity = qmdd_stack_matrix(identity, k+1, GATEID_I);

        // Check if this qubit should be controlled on 0 or 1
        if (x[k] == 1)
            ccphase = qmdd_stack_control(identity, ccphase, k);
        else if (x[k] == 0)
            ccphase = qmdd_stack_control(ccphase, identity, k);
    }

    return ccphase;
}

/**************************</Initial state creation>***************************/





/*******************************<Applying gates>*******************************/

static int periodic_gc_nodetable = 0; // trigger for gc of node table
static uint64_t gate_counter = 0;

static void
qmdd_do_before_gate(QMDD* qmdd)
{
    // check if ctable needs gc
    if (aadd_test_gc_wgt_table()) {
        aadd_protect(qmdd);
        aadd_gc_wgt_table();
        aadd_unprotect(qmdd);
    }

    if (periodic_gc_nodetable) {
        gate_counter++;
        if (gate_counter % periodic_gc_nodetable == 0) {
            aadd_protect(qmdd);
            sylvan_gc();
            aadd_unprotect(qmdd);
        }
    }

    // log stuff (if logging is enabled)
    qmdd_stats_log(*qmdd);
}

static void swap(BDDVAR *a, BDDVAR *b)
{
    BDDVAR tmp = *a;
    *a = *b;
    *b = tmp;
}

static bool
check_ctrls_before_targ(BDDVAR *c1, BDDVAR *c2, BDDVAR *c3, BDDVAR t)
{
    // sort controls
    if (*c1 > *c3) swap(c1, c2);
    if (*c1 > *c2) swap(c1, c2);
    if (*c2 > *c3) swap(c2, c3);

    // check if controls before target
    if (*c1 != AADD_INVALID_VAR && *c1 > t) return false;
    if (*c2 != AADD_INVALID_VAR && *c2 > t) return false;
    if (*c3 != AADD_INVALID_VAR && *c3 > t) return false;
    return true;
}

/* Wrapper for applying a single qubit gate. */
TASK_IMPL_3(QMDD, qmdd_gate, QMDD, qmdd, gate_id_t, gate, BDDVAR, target)
{
    qmdd_do_before_gate(&qmdd);
    return qmdd_gate_rec(qmdd, gate, target);
}

/* Wrapper for applying controlled gates with 1, 2, or 3 control qubits. */
QMDD _qmdd_cgate(QMDD state, gate_id_t gate, BDDVAR c1, BDDVAR c2, BDDVAR c3, BDDVAR t, BDDVAR n)
{
    qmdd_do_before_gate(&state);

    if (check_ctrls_before_targ(&c1, &c2, &c3, t)) {
        BDDVAR cs[4] = {c1, c2, c3, AADD_INVALID_VAR}; // last pos is to mark end
        return qmdd_cgate_rec(state, gate, cs, t);
    }
    else {
        assert(n != 0 && "ERROR: when ctrls > targ, nqubits must be passed to cgate() function.");
        int *c_options = malloc(sizeof(int)*(n+1));
        for (uint32_t k = 0; k < n; k++) c_options[k] = -1;
        if (c1 != AADD_INVALID_VAR && c1 < n) c_options[c1] = 1;
        if (c2 != AADD_INVALID_VAR && c2 < n) c_options[c2] = 1;
        if (c3 != AADD_INVALID_VAR && c3 < n) c_options[c3] = 1;
        if (t  != AADD_INVALID_VAR && t  < n) c_options[t] = 2;
        QMDD gate_matrix = qmdd_create_multi_cgate(n, c_options, gate);
        free(c_options);
        return aadd_matvec_mult(gate_matrix, state, n);
    }
}

/* Wrapper for applying a controlled gate where the controls are a range. */
TASK_IMPL_5(QMDD, qmdd_cgate_range, QMDD, qmdd, gate_id_t, gate, BDDVAR, c_first, BDDVAR, c_last, BDDVAR, t)
{
    qmdd_do_before_gate(&qmdd);
    return qmdd_cgate_range_rec(qmdd,gate,c_first,c_last,t);
}

TASK_IMPL_3(QMDD, qmdd_gate_rec, QMDD, q, gate_id_t, gate, BDDVAR, target)
{
    // Trivial cases
    if (AADD_WEIGHT(q) == AADD_ZERO) return q;

    BDDVAR var;
    QMDD res, low, high;
    aadd_get_topvar(q, target, &var, &low, &high);
    assert(var <= target);

    // Check cache
    bool cachenow = ((var % granularity) == 0);
    if (cachenow) {
        if (cache_get3(CACHE_QMDD_GATE, sylvan_false, AADD_TARGET(q), GATE_OPID_40(gate, target, 0), &res)) {
            sylvan_stats_count(QMDD_GATE_CACHED);
            // Multiply root of res with root of input qmdd
            AMP new_root_amp = wgt_mul(AADD_WEIGHT(q), AADD_WEIGHT(res));
            res = aadd_bundle(AADD_TARGET(res), new_root_amp);
            return res;
        }
    }

    if (var == target) {
        AMP a_u00 = wgt_mul(AADD_WEIGHT(low), gates[gate][0]);
        AMP a_u10 = wgt_mul(AADD_WEIGHT(low), gates[gate][2]);
        AMP b_u01 = wgt_mul(AADD_WEIGHT(high), gates[gate][1]);
        AMP b_u11 = wgt_mul(AADD_WEIGHT(high), gates[gate][3]);
        QMDD low1, low2, high1, high2;
        low1  = aadd_bundle(AADD_TARGET(low), a_u00);
        low2  = aadd_bundle(AADD_TARGET(high),b_u01);
        high1 = aadd_bundle(AADD_TARGET(low), a_u10);
        high2 = aadd_bundle(AADD_TARGET(high),b_u11);
        aadd_refs_spawn(SPAWN(aadd_plus, high1, high2));
        low = aadd_refs_push(CALL(aadd_plus, low1, low2));
        high = aadd_refs_sync(SYNC(aadd_plus));
        aadd_refs_pop(1);
        res = aadd_makenode(target, low, high);
    }
    else { // var < target: not at target qubit yet, recursive calls down
        aadd_refs_spawn(SPAWN(qmdd_gate_rec, high, gate, target));
        low = aadd_refs_push(CALL(qmdd_gate_rec, low, gate, target));
        high = aadd_refs_sync(SYNC(qmdd_gate_rec));
        aadd_refs_pop(1);
        res  = aadd_makenode(var, low, high);
    }

    // Store not yet "root normalized" result in cache
    if (cachenow) {
        if (cache_put3(CACHE_QMDD_GATE, sylvan_false, AADD_TARGET(q), GATE_OPID_40(gate, target, 0), res)) 
            sylvan_stats_count(QMDD_GATE_CACHEDPUT);
    }
    // Multiply amp res with amp of input qmdd
    AMP new_root_amp = wgt_mul(AADD_WEIGHT(q), AADD_WEIGHT(res));
    res = aadd_bundle(AADD_TARGET(res), new_root_amp);
    return res;
}

TASK_IMPL_5(QMDD, qmdd_cgate_rec, QMDD, q, gate_id_t, gate, BDDVAR*, cs, uint32_t, ci, BDDVAR, t)
{
    // Get current control qubit. If no more control qubits, apply gate here
    BDDVAR c = cs[ci];
    if (c == AADD_INVALID_VAR || ci > MAX_CONTROLS) {
        return CALL(qmdd_gate_rec, q, gate, t);
    }

    assert(c < t && "ctrl < target required");
    if (ci > 0) 
        assert(cs[ci-1] < cs[ci]  && "order required for multiple controls");

    BDDVAR var;
    QMDD res, low, high;
    aadd_get_topvar(q, c, &var, &low, &high);
    assert(var <= c);

    // Check cache
    bool cachenow = ((var % granularity) == 0);
    if (cachenow) {
        if (cache_get3(CACHE_QMDD_CGATE, sylvan_false, AADD_TARGET(q),
                    GATE_OPID_64(gate, ci, cs[0], cs[1], cs[2], t),
                    &res)) {
            sylvan_stats_count(QMDD_CGATE_CACHED);
            // Multiply root amp of res with input root amp
            AMP new_root_amp = wgt_mul(AADD_WEIGHT(q), AADD_WEIGHT(res));
            res = aadd_bundle(AADD_TARGET(res), new_root_amp);
            return res;
        }
    }

    // If current node is (one of) the control qubit(s), 
    // control on q_c = |1> (high edge)
    if (var == c) {
        high = CALL(qmdd_cgate_rec, high, gate, cs, ci+1, t);
    }
    // Not at control qubit yet, apply to both childeren.
    else {
        aadd_refs_spawn(SPAWN(qmdd_cgate_rec, high, gate, cs, ci, t));
        low = aadd_refs_push(CALL(qmdd_cgate_rec, low, gate, cs, ci, t));
        high = aadd_refs_sync(SYNC(qmdd_cgate_rec));
        aadd_refs_pop(1);
    }
    res = aadd_makenode(var, low, high);

    // Store not yet "root normalized" result in cache
    if (cachenow) {
        if (cache_put3(CACHE_QMDD_CGATE, sylvan_false, AADD_TARGET(q),
                    GATE_OPID_64(gate, ci, cs[0], cs[1], cs[2], t),
                    res)) {
            sylvan_stats_count(QMDD_CGATE_CACHEDPUT);
        }
    }
    // Multiply root amp of res with input root amp
    AMP new_root_amp = wgt_mul(AADD_WEIGHT(q), AADD_WEIGHT(res));
    res = aadd_bundle(AADD_TARGET(res), new_root_amp);
    return res;
}

TASK_IMPL_6(QMDD, qmdd_cgate_range_rec, QMDD, q, gate_id_t, gate, BDDVAR, c_first, BDDVAR, c_last, BDDVAR, t, BDDVAR, k)
{
    // Past last control (done with "control part" of controlled gate)
    if (k > c_last) {
        return CALL(qmdd_gate_rec, q, gate, t);
    }

    assert(c_first <= c_last);
    assert(c_last < t);

    // Check cache
    QMDD res;
    bool cachenow = ((k % granularity) == 0);
    if (cachenow) {
        if (cache_get3(CACHE_QMDD_CGATE_RANGE, sylvan_false, AADD_TARGET(q),
                       GATE_OPID_64(gate, c_first, c_last, k, t, 0),
                       &res)) {
            sylvan_stats_count(QMDD_CGATE_CACHED);
            // Multiply root amp of result with the input root amp
            AMP new_root_amp = wgt_mul(AADD_WEIGHT(q), AADD_WEIGHT(res));
            res = aadd_bundle(AADD_TARGET(res), new_root_amp);
            return res;
        }
    }

    // Get top node
    BDDVAR var, nextvar;
    QMDD low, high;
    if (k < c_first) {
        // possibly skip to c_first if we are before c_first
        aadd_get_topvar(q, c_first, &var, &low, &high);
        assert(var <= c_first);
    }
    else {
        // k is a control qubit, so get node with var = k (insert if skipped)
        assert(c_first <= k && k <= c_last);
        aadd_get_topvar(q, k, &var, &low, &high);
        assert(var == k);
    }
    nextvar = var + 1;
    
    // Not at first control qubit yet, apply to both children
    if (var < c_first) {
        aadd_refs_spawn(SPAWN(qmdd_cgate_range_rec, high, gate, c_first, c_last, t, nextvar));
        low = aadd_refs_push(CALL(qmdd_cgate_range_rec, low, gate, c_first, c_last, t, nextvar));
        high = aadd_refs_sync(SYNC(qmdd_cgate_range_rec));
        aadd_refs_pop(1);
    }
    // Current var is a control qubit, control on q_k = |1> (high edge)
    else {
        high = CALL(qmdd_cgate_range_rec, high, gate, c_first, c_last, t, nextvar);
    }
    res = aadd_makenode(var, low, high);

    // Store not yet "root normalized" result in cache
    if (cachenow) {
        if (cache_put3(CACHE_QMDD_CGATE_RANGE, sylvan_false, AADD_TARGET(q),
                       GATE_OPID_64(gate, c_first, c_last, k, t, 0),
                       res)) {
            sylvan_stats_count(QMDD_CGATE_CACHEDPUT);
        }
    }
    // Multiply root amp of result with the input root amp
    AMP new_root_amp = wgt_mul(AADD_WEIGHT(q), AADD_WEIGHT(res));
    res = aadd_bundle(AADD_TARGET(res), new_root_amp);
    return res;
}

/******************************</Applying gates>*******************************/





/*********************<Applying (controlled) sub-circuits>*********************/

QMDD
qmdd_circuit_swap(QMDD qmdd, BDDVAR qubit1, BDDVAR qubit2)
{
    if (qubit1 > qubit2) {
        BDDVAR tmp = qubit2;
        qubit2 = qubit1;
        qubit1 = tmp;
    }

    QMDD res;

    // CNOT
    res = qmdd_cgate(qmdd, GATEID_X, qubit1, qubit2);
    // upside down CNOT (equivalent)
    res = qmdd_gate(res, GATEID_H, qubit1);
    res = qmdd_cgate(res, GATEID_Z, qubit1, qubit2);
    res = qmdd_gate(res, GATEID_H, qubit1);
    // CNOT
    res = qmdd_cgate(res, GATEID_X, qubit1, qubit2);

    return res;
}

QMDD
qmdd_circuit_reverse_range(QMDD qmdd, BDDVAR first, BDDVAR last)
{
    QMDD res = qmdd;
    BDDVAR a, b;
    int num_qubits = (last - first) + 1;
    for (int j = 0; j < (int)(num_qubits/2); j++) {
        a = first + j;
        b = last  - j;
        res = qmdd_circuit_swap(res, a, b);
    }
    return res;
}

QMDD
qmdd_circuit_QFT(QMDD qmdd, BDDVAR first, BDDVAR last)
{
    int k;
    QMDD res = qmdd;
    BDDVAR a, b;
    for (a = first; a <= last; a++) {
        
        // H gate on current qubit
        res = qmdd_gate(res, GATEID_H, a);

        // Controlled phase gates on all qubits below
        for (b = a+1; b <= last; b++) {
            k = (b - a) + 1;
            res = qmdd_cgate(res, GATEID_Rk(k), a, b);
        }
    }

    // Note that we're not swapping the qubit order in this function

    return res;
}

QMDD
qmdd_circuit_QFT_inv(QMDD qmdd, BDDVAR first, BDDVAR last)
{
    int k;
    QMDD res = qmdd;
    BDDVAR a, b;

    // Note that we're not swapping the qubit order in this function
    
    // H gates and phase gates (but now backwards)
    for (a = last + 1; a-- > first; ) { // weird for-loop because BDDVARs are unsigned

        // Controlled phase gates (negative angles this time)
        for (b = last; b >= (a+1); b--){
            k = (b - a) + 1;
            res = qmdd_cgate(res, GATEID_Rk_dag(k), a, b);
        }

        // H on current qubit
        res = qmdd_gate(res, GATEID_H, a);
    }

    return res;
}

QMDD
qmdd_circuit(QMDD qmdd, circuit_id_t circ_id, BDDVAR t1, BDDVAR t2)
{
    switch (circ_id) {  // don't judge me please
        case CIRCID_swap          : return qmdd_circuit_swap(qmdd, t1, t2);
        case CIRCID_reverse_range : return qmdd_circuit_reverse_range(qmdd, t1, t2);
        case CIRCID_QFT           : return qmdd_circuit_QFT(qmdd, t1, t2);
        case CIRCID_QFT_inv       : return qmdd_circuit_QFT_inv(qmdd, t1, t2);
        default :
            assert ("Invalid circuit ID" && false);
            return AADD_TERMINAL;
    }
}

TASK_IMPL_6(QMDD, qmdd_ccircuit, QMDD, qmdd, circuit_id_t, circ_id, BDDVAR*, cs, uint32_t, ci, BDDVAR, t1, BDDVAR, t2)
{
    // Cache lookup
    QMDD res;
    bool cachenow = 1;
    if (cachenow) {
        if (cache_get3(CACHE_QMDD_SUBCIRC, GATE_OPID_40(0, ci, 0), qmdd, 
                       GATE_OPID_64(circ_id, cs[0], cs[1], cs[2], t1, t2),
                       &res)) {
            return res;
        }
    }

    // Get current control qubit
    BDDVAR c = cs[ci];

    // If no more control qubits, apply sub-circ here
    if (c == AADD_INVALID_VAR || ci > MAX_CONTROLS) {
        res = qmdd_circuit(qmdd, circ_id, t1, t2);
        // the gates in qmdd_circuit already took care of multiplying the input 
        // root amp with normalization, so no need to do that here again
    }
    else {
        BDDVAR var;
        QMDD low, high;
        aadd_get_topvar(qmdd, c, &var, &low, &high);
        assert(var <= c);

        if (var == c) {
            ci++; // next control qubit
            high = CALL(qmdd_ccircuit, high, circ_id, cs, ci, t1, t2);
            ci--;
        }
        else {
            // recursive call to both children
            aadd_refs_spawn(SPAWN(qmdd_ccircuit, high, circ_id, cs, ci, t1, t2));
            low = CALL(qmdd_ccircuit, low, circ_id, cs, ci, t1, t2);
            aadd_refs_push(low);
            high = aadd_refs_sync(SYNC(qmdd_ccircuit));
            aadd_refs_pop(1);
        }
        res = aadd_makenode(var, low, high);
        // Multiply root amp of sum with input root amp 
        AMP new_root_amp = wgt_mul(AADD_WEIGHT(qmdd), AADD_WEIGHT(res));
        res = aadd_bundle(AADD_TARGET(res), new_root_amp);
    }
    
    // Add to cache, return
    if (cachenow) {
        cache_put3(CACHE_QMDD_SUBCIRC, GATE_OPID_40(0, ci, 0), qmdd, 
                   GATE_OPID_64(circ_id, cs[0], cs[1], cs[2], t1, t2), 
                   res);
    }
    return res;
}

QMDD
qmdd_all_control_phase_rec(QMDD qmdd, BDDVAR k, BDDVAR n, bool *x)
{
    assert(k < n);
    
    bool skipped_k = false;
    aaddnode_t node;
    if (AADD_TARGET(qmdd) == AADD_TERMINAL) {
        skipped_k = true;
    }
    else {
        node = AADD_GETNODE(AADD_TARGET(qmdd));
        BDDVAR var = aaddnode_getvar(node);
        if(var > k) {
            skipped_k = true;
        }
    }

    QMDD low, high;
    if (skipped_k) {
        // insert skipped node
        low  = aadd_bundle(AADD_TARGET(qmdd), AADD_ONE);
        high = aadd_bundle(AADD_TARGET(qmdd), AADD_ONE);
    }
    else {
        // case var == k (var < k shouldn't happen)
        aaddnode_getchilderen(node, &low, &high);
    }

    // terminal case, apply phase depending on x[k] (control k on 0 or 1)
    if (k == (n-1)) {
        if (x[k] == 1) {
            AMP new_amp = wgt_mul(AADD_WEIGHT(high), AADD_MIN_ONE);
            high = aadd_bundle(AADD_TARGET(high), new_amp);
        }
        else {
            AMP new_amp = wgt_mul(AADD_WEIGHT(low), AADD_MIN_ONE);
            low = aadd_bundle(AADD_TARGET(low), new_amp);
        }
    }
    // non terminal case, choose low/high depending on x[k] (control k on 0 or 1)
    else {
        if (x[k] == 1) {
            k++; // next level
            high = qmdd_all_control_phase_rec(high, k, n, x);
            k--;
        }
        else {
            k++;
            low = qmdd_all_control_phase_rec(low, k, n, x);
            k--;
        }
    }

    QMDD res = aadd_makenode(k, low, high);

    // multiply by existing edge weight on qmdd
    AMP new_root_amp = wgt_mul(AADD_WEIGHT(qmdd), AADD_WEIGHT(res));
    res = aadd_bundle(AADD_TARGET(res), new_root_amp);
    return res;
}

QMDD
qmdd_all_control_phase(QMDD qmdd, BDDVAR n, bool *x)
{
    return qmdd_all_control_phase_rec(qmdd, 0, n, x);
}

/********************</Applying (controlled) sub-circuits>*********************/





/***********************<Measurements and probabilities>***********************/

QMDD
qmdd_measure_qubit(QMDD qmdd, BDDVAR k, BDDVAR nvars, int *m, double *p)
{
    if (k == 0) return qmdd_measure_q0(qmdd, nvars, m, p);
    qmdd = qmdd_circuit_swap(qmdd, 0, k);
    qmdd = qmdd_measure_q0(qmdd, nvars, m, p);
    qmdd = qmdd_circuit_swap(qmdd, 0, k);
    return qmdd;
}

QMDD
qmdd_measure_q0(QMDD qmdd, BDDVAR nvars, int *m, double *p)
{  
    // get probabilities for q0 = |0> and q0 = |1>
    double prob_low, prob_high, prob_root;

    QMDD low, high;
    BDDVAR var;
    aadd_get_topvar(qmdd, 0, &var, &low, &high);

    if (testing_mode) assert(qmdd_is_unitvector(qmdd, nvars));

    // TODO: don't use doubles here but allow for mpreal ?
    // (e.g. by using AMPs)
    prob_low  = qmdd_unnormed_prob(low,  1, nvars);
    prob_high = qmdd_unnormed_prob(high, 1, nvars);
    prob_root = qmdd_amp_to_prob(AADD_WEIGHT(qmdd));
    prob_low  *= prob_root;
    prob_high *= prob_root;
    if (fabs(prob_low + prob_high - 1.0) > 1e-6) {
        fprintf(stderr, "WARNING: prob sum = %.14lf\n", prob_low + prob_high);
    }

    // flip a coin
    float rnd = ((float)rand())/((float)RAND_MAX);
    *m = (rnd < prob_low) ? 0 : 1;
    *p = prob_low;

    // produce post-measurement state
    AMP norm;
    if (*m == 0) {
        high = aadd_bundle(AADD_TERMINAL, AADD_ZERO);
        norm = qmdd_amp_from_prob(prob_low);
    }
    else {
        low  = aadd_bundle(AADD_TERMINAL, AADD_ZERO);
        norm = qmdd_amp_from_prob(prob_high);
    }

    QMDD res = aadd_makenode(0, low, high);

    AMP new_root_amp = wgt_mul(AADD_WEIGHT(qmdd), AADD_WEIGHT(res));
    new_root_amp     = wgt_div(new_root_amp, norm);

    res = aadd_bundle(AADD_TARGET(res), new_root_amp);
    res = qmdd_remove_global_phase(res);
    return res;
}

QMDD
qmdd_measure_all(QMDD qmdd, BDDVAR n, bool* ms, double *p)
{
    aaddnode_t node;
    bool skipped;
    BDDVAR var;
    double prob_low, prob_high, prob_path = 1.0, prob_roots = 1.0;

    for (BDDVAR k=0; k < n; k++) {
        // find relevant node (assuming it should be the next one)
        skipped = false;
        if (AADD_TARGET(qmdd) == AADD_TERMINAL) {
            skipped = true;
        }
        else {
            node = AADD_GETNODE(AADD_TARGET(qmdd));
            var = aaddnode_getvar(node);
            assert(var >= k);
            if (var > k) skipped = true;
        }
        QMDD low, high;
        if (skipped) {
            // if skipped q0 is a don't care, treat separately?
            low  = aadd_bundle(AADD_TARGET(qmdd), AADD_ONE);
            high = aadd_bundle(AADD_TARGET(qmdd), AADD_ONE);
        }
        else {
            aaddnode_getchilderen(node, &low, &high);
        }

        prob_low  = qmdd_unnormed_prob(low,  k+1, n);
        prob_high = qmdd_unnormed_prob(high, k+1, n);
        prob_roots *= qmdd_amp_to_prob(AADD_WEIGHT(qmdd));
        prob_high = prob_high * prob_roots / prob_path;
        prob_low  = prob_low  * prob_roots / prob_path;

        if (fabs(prob_low + prob_high - 1.0) > sylvan_edge_weights_tolerance()) {
            fprintf(stderr, "WARNING: prob sum = %.14lf\n", prob_low + prob_high);
        }

        // flip a coin
        float rnd = ((float)rand())/((float)RAND_MAX);
        ms[k] = (rnd < prob_low) ? 0 : 1;

        // Get next edge
        qmdd        = (ms[k] == 0) ? low : high;
        prob_path *= (ms[k] == 0) ? prob_low : prob_high;
    }

    *p = prob_path;

    //TODO: replace below ith "return qmdd_create_basis_state(n, ms);""
    QMDD low, high, prev = AADD_TERMINAL;

    for (int k = n-1; k >= 0; k--) {
        if (ms[k] == 0) {
            low = aadd_bundle(AADD_TARGET(prev), AADD_ONE);
            high = aadd_bundle(AADD_TERMINAL, AADD_ZERO);
        }
        else if (ms[k] == 1) {
            low = aadd_bundle(AADD_TERMINAL, AADD_ZERO);
            high = aadd_bundle(AADD_TARGET(prev), AADD_ONE);
        }
        // add node to unique table
        prev = aadd_makenode(k, low, high);
    }
    return prev;
}

// Container for disguising doubles as ints so they can go in Sylvan's cache
// (see also union "hack" in mtbdd_satcount)
typedef union {
    double   as_double;
    uint64_t as_int;
} double_hack_t;

TASK_IMPL_3(double, qmdd_unnormed_prob, QMDD, qmdd, BDDVAR, topvar, BDDVAR, nvars)
{
    assert(topvar <= nvars);

    if (topvar == nvars) {
        assert(AADD_TARGET(qmdd) == AADD_TERMINAL);
        return qmdd_amp_to_prob(AADD_WEIGHT(qmdd));
    }

    // Look in cache
    bool cachenow = 1;
    if (cachenow) {
        uint64_t prob_bits;
        if (cache_get3(CACHE_QMDD_PROB, sylvan_false, qmdd, QMDD_PARAM_PACK_16(topvar, nvars), &prob_bits)) {
            sylvan_stats_count(QMDD_PROB_CACHED);
            double_hack_t container = (double_hack_t) prob_bits;
            return container.as_double;
        }
    }

    // Check if the node we want is being skipped
    BDDVAR var;
    QMDD low, high;
    aadd_get_topvar(qmdd, topvar, &var, &low, &high);

    double prob_low, prob_high, prob_root, prob_res; // "prob" = absolute amps squared
    BDDVAR nextvar = topvar + 1;

    SPAWN(qmdd_unnormed_prob, high, nextvar, nvars);
    prob_low  = CALL(qmdd_unnormed_prob, low, nextvar, nvars);
    prob_high = SYNC(qmdd_unnormed_prob);
    prob_root = qmdd_amp_to_prob(AADD_WEIGHT(qmdd));
    prob_res  = prob_root * (prob_low + prob_high);

    // Put in cache and return
    if (cachenow) {
        double_hack_t container = (double_hack_t) prob_res;
        if (cache_put3(CACHE_QMDD_PROB, sylvan_false, qmdd, QMDD_PARAM_PACK_16(topvar, nvars), container.as_int))
            sylvan_stats_count(QMDD_PROB_CACHEDPUT);
    }
    return prob_res;
}

complex_t
qmdd_get_amplitude(QMDD q, bool *x, BDDVAR nqubits)
{
    // QMDD is indexed q_0, ..., q_{n-1} but |x> is assumed q_{n-1}, ..., q_0,
    // so we temporarily reverse x.
    reverse_bit_array(x, nqubits);
    complex_t res;
    weight_value(aadd_getvalue(q, x), &res);
    reverse_bit_array(x, nqubits);
    return res;
}

double
qmdd_amp_to_prob(AMP a)
{
    complex_t c;
    weight_value(a, &c);
    double abs = flt_sqrt( c.r*c.r + c.i*c.i );
    return (abs*abs);
}

AMP
qmdd_amp_from_prob(double a)
{
    complex_t c;
    c.r = flt_sqrt(a);
    c.i = 0;
    return weight_lookup(&c);
}

/**********************</Measurements and probabilities>***********************/





/*******************************<Miscellaneous>********************************/

QMDD
qmdd_remove_global_phase(QMDD qmdd)
{
    // Remove global phase by replacing amp of qmdd with absolute value of amp
    AMP abs = wgt_abs(AADD_WEIGHT(qmdd));
    QMDD res = aadd_bundle(AADD_TARGET(qmdd), abs);
    return res;
}

void
qmdd_set_periodic_gc_nodetable(int every_n_gates)
{
    periodic_gc_nodetable = every_n_gates;
}

/******************************</Miscellaneous>********************************/





/*******************************<Logging stats>********************************/

bool qmdd_stats_logging = false;
uint32_t statslog_granularity = 1;
uint64_t statslog_buffer = 10; // TODO: remove manual buffer flushing
FILE *qmdd_logfile;
uint64_t nodes_peak = 0;
double nodes_avg = 0;
uint64_t logcounter = 0;
uint64_t logtrycounter = 0;

void
qmdd_stats_start(FILE *out)
{
    if (out == NULL) return;
    qmdd_stats_logging = true;
    qmdd_logfile = out;
    fprintf(qmdd_logfile, "nodes, amps\n");
    nodes_peak = 0;
    logcounter = 0;
    logtrycounter = 0;
}

void
qmdd_stats_set_granularity(uint32_t g)
{
    if (g == 0) statslog_granularity = 1;
    else statslog_granularity = g;
}

void
qmdd_stats_log(QMDD qmdd)
{
    if (!qmdd_stats_logging) return;

    // only log every 'statslog_granularity' calls of this function
    if (logtrycounter++ % statslog_granularity != 0) return;

    // Insert info
    uint64_t num_nodes = aadd_countnodes(qmdd);
    uint64_t num_amps  = sylvan_edge_weights_count_entries();
    fprintf(qmdd_logfile, "%" PRIu64 ",%" PRIu64 "\n", num_nodes, num_amps);
    logcounter++;

    // manually flush every 'statslog_buffer' entries
    if (logcounter % statslog_buffer == 0)
        fflush(qmdd_logfile);

    // peak nodes
    if (num_nodes > nodes_peak)
        nodes_peak = num_nodes;
    
    // (online) avg nodes
    double a = 1.0/(double)logcounter;
    double b = 1.0 - a;
    nodes_avg = a * (double)num_nodes + b * nodes_avg;
}

uint64_t
qmdd_stats_get_nodes_peak()
{
    return nodes_peak;
}

double
qmdd_stats_get_nodes_avg()
{
    return nodes_avg;
}

uint64_t
qmdd_stats_get_logcounter()
{
    return logtrycounter;
}

void
qmdd_stats_finish()
{
    if (!qmdd_stats_logging) return;
    fflush(qmdd_logfile);
    qmdd_stats_logging = false;
    nodes_peak = 0;
    logcounter = 0;
    logtrycounter = 0;
}

/******************************</Logging stats>********************************/





/********************************<Debug stuff>*********************************/

void
qmdd_set_testing_mode(bool on)
{
    testing_mode = on;
}

bool
qmdd_is_close_to_unitvector(QMDD qmdd, BDDVAR n, double tol)
{
    bool WRITE_TO_FILE = false;
    bool has_next = true;
    AMP a;
    bool x[n];
    for(BDDVAR k=0; k<n; k++) x[k] = 0;

    double sum_abs_squares = 0.0;
    while(has_next){
        a = aadd_getvalue(qmdd, x);
        sum_abs_squares += qmdd_amp_to_prob(a);
        has_next = _next_bitstring(x, n);
    }

    if (fabs(sum_abs_squares - 1.0) < tol) {
        if (WRITE_TO_FILE) {
            FILE *fp;
            fp = fopen("is_unitvector_true.dot", "w");
            aadd_fprintdot(fp, qmdd, false);
            fclose(fp);
        }
        return true;
    }
    else {
        //printf("probs sum to %.30lf\n", sum_abs_squares);
        if (WRITE_TO_FILE) {
            FILE *fp;
            fp = fopen("is_unitvector_false.dot", "w");
            aadd_fprintdot(fp, qmdd, false);
            fclose(fp);
        }
        return false;
    }
}

bool
qmdd_is_unitvector(QMDD qmdd, BDDVAR n)
{
    return qmdd_is_close_to_unitvector(qmdd, n, sylvan_edge_weights_tolerance()*10);
}

double
qmdd_get_magnitude(QMDD qmdd, BDDVAR n)
{
    return qmdd_unnormed_prob(qmdd, 0, n);
}

/*******************************</Debug stuff>*********************************/
