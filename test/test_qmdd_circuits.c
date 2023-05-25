#include <inttypes.h>
#include <stdio.h>
#include <time.h>

#include "qsylvan.h"
#include "test_assert.h"

bool VERBOSE = true;

int test_swap_circuit()
{
    QMDD q;
    bool x3[] = {0,0,0};
    AMP a;

    q = qmdd_create_basis_state(3, x3);
    q = qmdd_gate(q, GATEID_X, 2);
    x3[2] = 0; x3[1] = 0; x3[0] = 0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 0; x3[1] = 0; x3[0] = 1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 1; x3[1] = 0; x3[0] = 0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ONE);
    x3[2] = 1; x3[1] = 0; x3[0] = 1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    test_assert(aadd_is_ordered(q, 3));

    q = qmdd_circuit_swap(q, 0, 2);
    x3[2] = 0; x3[1] = 0; x3[0] = 0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 0; x3[1] = 0; x3[0] = 1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ONE);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 1; x3[1] = 0; x3[0] = 0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 1; x3[1] = 0; x3[0] = 1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    test_assert(aadd_is_ordered(q, 3));

    q = qmdd_circuit_swap(q, 0, 1);
    x3[2] = 0; x3[1] = 0; x3[0] = 0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 0; x3[1] = 0; x3[0] = 1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ONE);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 1; x3[1] = 0; x3[0] = 0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 1; x3[1] = 0; x3[0] = 1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    test_assert(aadd_is_ordered(q, 3));

    // TODO: more tests

    if(VERBOSE) printf("qmdd swap gates:           ok\n");
    return 0;
}

int test_cswap_circuit()
{
    QMDD q;
    bool x3[] = {0,0,0};
    AMP a;

    uint32_t cs[3];
    cs[0] = 0; cs[1] = AADD_INVALID_VAR; cs[2] = AADD_INVALID_VAR; // control only on q0

    x3[2]=1; x3[1]=0; x3[0]=0; 
    q = qmdd_create_basis_state(3, x3);
    q = qmdd_ccircuit(q, CIRCID_swap, cs, 1, 2); // control is |0>, nothing should happen
    x3[2] = 0; x3[1] = 0; x3[0] = 0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 0; x3[1] = 0; x3[0] = 1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 1; x3[1] = 0; x3[0] = 0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ONE);
    x3[2] = 1; x3[1] = 0; x3[0] = 1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);

    x3[2]=1; x3[1]=0; x3[0]=1; 
    q = qmdd_create_basis_state(3, x3);
    q = qmdd_ccircuit(q, CIRCID_swap, cs, 1, 2); // control is |1>, should swap q1 and q2
    x3[2] = 0; x3[1] = 0; x3[0] = 0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 0; x3[1] = 0; x3[0] = 1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ONE);
    x3[2] = 1; x3[1] = 0; x3[0] = 0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 1; x3[1] = 0; x3[0] = 1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);

    x3[2]=1; x3[1]=0; x3[0]=0; 
    q = qmdd_create_basis_state(3, x3);
    q = qmdd_gate(q, GATEID_H, 0); // input state: |10+> = 1/sqrt(2)(|100> + |101>)
    x3[2] = 0; x3[1] = 0; x3[0] = 0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 0; x3[1] = 0; x3[0] = 1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 1; x3[1] = 0; x3[0] = 0; a = aadd_getvalue(q, x3); test_assert(a == weight_lookup(cmake(1.0/flt_sqrt(2.0),0)));
    x3[2] = 1; x3[1] = 0; x3[0] = 1; a = aadd_getvalue(q, x3); test_assert(a == weight_lookup(cmake(1.0/flt_sqrt(2.0),0)));
    x3[2] = 1; x3[1] = 1; x3[0] = 0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    q = qmdd_ccircuit(q, CIRCID_swap, cs, 1, 2); // control is |+>, expected output: 1/sqrt(2)(|100> + |011>)
    x3[2] = 0; x3[1] = 0; x3[0] = 0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 0; x3[1] = 0; x3[0] = 1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; a = aadd_getvalue(q, x3); test_assert(a == weight_lookup(cmake(1.0/flt_sqrt(2.0),0)));
    x3[2] = 1; x3[1] = 0; x3[0] = 0; a = aadd_getvalue(q, x3); test_assert(a == weight_lookup(cmake(1.0/flt_sqrt(2.0),0)));
    x3[2] = 1; x3[1] = 0; x3[0] = 1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);

    x3[2]=1; x3[1]=0; x3[0]=0; 
    q = qmdd_create_basis_state(3, x3);
    q = qmdd_gate(q, GATEID_H, 0); 
    q = qmdd_gate(q, GATEID_Z, 0); // input state: |10-> = 1/sqrt(2)(|100> - |101>)
    x3[2] = 0; x3[1] = 0; x3[0] = 0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 0; x3[1] = 0; x3[0] = 1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 1; x3[1] = 0; x3[0] = 0; a = aadd_getvalue(q, x3); test_assert(a == weight_lookup(cmake(1.0/flt_sqrt(2.0),0)));
    x3[2] = 1; x3[1] = 0; x3[0] = 1; a = aadd_getvalue(q, x3); test_assert(a == weight_lookup(cmake(-1.0/flt_sqrt(2.0),0)));
    x3[2] = 1; x3[1] = 1; x3[0] = 0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    q = qmdd_ccircuit(q, CIRCID_swap, cs, 1, 2); // control is |->, expected output: 1/sqrt(2)(|100> - |011>)
    x3[2] = 0; x3[1] = 0; x3[0] = 0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 0; x3[1] = 0; x3[0] = 1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; a = aadd_getvalue(q, x3); test_assert(a == weight_lookup(cmake(-1.0/flt_sqrt(2.0),0)));
    x3[2] = 1; x3[1] = 0; x3[0] = 0; a = aadd_getvalue(q, x3); test_assert(a == weight_lookup(cmake(1.0/flt_sqrt(2.0),0)));
    x3[2] = 1; x3[1] = 0; x3[0] = 1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);

    x3[2]=1; x3[1]=0; x3[0]=0; 
    q = qmdd_create_basis_state(3, x3);
    q = qmdd_gate(q, GATEID_H, 0); 
    q = qmdd_gate(q, GATEID_S, 0); // input state: |10>|+i> = 1/sqrt(2)(|100> + i|101>)
    x3[2] = 0; x3[1] = 0; x3[0] = 0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 0; x3[1] = 0; x3[0] = 1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 1; x3[1] = 0; x3[0] = 0; a = aadd_getvalue(q, x3); test_assert(a == weight_lookup(cmake(1.0/flt_sqrt(2.0),0)));
    x3[2] = 1; x3[1] = 0; x3[0] = 1; a = aadd_getvalue(q, x3); test_assert(a == weight_lookup(cmake(0,1.0/flt_sqrt(2.0))));
    x3[2] = 1; x3[1] = 1; x3[0] = 0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    q = qmdd_ccircuit(q, CIRCID_swap, cs, 1, 2); // control is |+i>, expected output: 1/sqrt(2)(|100> + i|011>)
    x3[2] = 0; x3[1] = 0; x3[0] = 0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 0; x3[1] = 0; x3[0] = 1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; a = aadd_getvalue(q, x3); test_assert(a == weight_lookup(cmake(0,1.0/flt_sqrt(2.0))));
    x3[2] = 1; x3[1] = 0; x3[0] = 0; a = aadd_getvalue(q, x3); test_assert(a == weight_lookup(cmake(1.0/flt_sqrt(2.0),0)));
    x3[2] = 1; x3[1] = 0; x3[0] = 1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 0; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 1; a = aadd_getvalue(q, x3); test_assert(a == AADD_ZERO);

    // TODO: more tests ?

    if(VERBOSE) printf("qmdd c-swap gates:         ok\n");
    return 0;
}

int test_tensor_product()
{
    QMDD q0, q1, qTest, qRef;
    bool x3[] = {0,1,0};
    bool x2[] = {1,0};
    bool x3_2[] = {0,1,0,1,0};
    bool x2_3[] = {1,0,0,1,0};

    q0 = qmdd_create_basis_state(3, x3);
    q1 = qmdd_create_basis_state(2, x2);
    test_assert(aadd_countnodes(q0) == 4);
    test_assert(aadd_countnodes(q1) == 3);

    // q0 (tensor) q1
    qTest = aadd_vec_tensor_prod(q0, q1, 3);
    qRef  = qmdd_create_basis_state(5, x3_2);
    test_assert(aadd_is_ordered(qTest, 5));
    test_assert(aadd_countnodes(qTest) == 6);
    test_assert(aadd_equivalent(qTest, qRef, 5, false, true));
    test_assert(aadd_equivalent(qTest, qRef, 5, true, true));
    test_assert(qTest == qRef);

    // q1 (tensor) q0
    qTest = aadd_vec_tensor_prod(q1, q0, 2);
    qRef  = qmdd_create_basis_state(5, x2_3);
    test_assert(aadd_is_ordered(qTest, 5));
    test_assert(aadd_countnodes(qTest) == 6);
    test_assert(aadd_equivalent(qTest, qRef, 5, false, true));
    test_assert(aadd_equivalent(qTest, qRef, 5, true, true));
    test_assert(qTest == qRef);


    // TODO: test on other states than basis states


    if(VERBOSE) printf("qmdd tensor (on vecs):     ok\n");
    return 0;
}

// measurement of q0 + sanity checks
int test_measure_random_state(QMDD qmdd, BDDVAR nvars)
{
    QMDD qm;
    int m; double p;

    test_assert(aadd_is_ordered(qmdd, nvars));
    test_assert(qmdd_is_unitvector(qmdd, nvars));
    qm = qmdd_measure_q0(qmdd, nvars, &m, &p);
    test_assert(aadd_is_ordered(qmdd, nvars));
    test_assert(qmdd_is_unitvector(qm, nvars));
    if ((flt_abs(p - 1.0) < cmap_get_tolerance()) || (flt_abs(p - 0.0) < cmap_get_tolerance())){
        qmdd = qmdd_remove_global_phase(qmdd); // measurement removes global phase
        test_assert(aadd_equivalent(qmdd, qm, nvars, false, false));
        test_assert(aadd_equivalent(qmdd, qm, nvars, true,  false));
        test_assert(qm == qmdd);
    } else {
        test_assert(aadd_countnodes(qm) < aadd_countnodes(qmdd));
    }

    return 0;
}

int test_measurements()
{
    QMDD q, qPM;
    AMP a;
    bool x2[] = {0,0};
    bool x3[] = {0,0,0};
    int m;
    int repeat = 10;
    double prob;
    srand(time(NULL));

    // kets labeld as |q2, q1, q0>
    for(int i=0; i < repeat; i++) {
        // |000> 
        x3[2]=0; x3[1]=0; x3[0]=0;
        q   = qmdd_create_basis_state(3, x3);
        qPM = qmdd_measure_qubit(q, 0, 3, &m, &prob);
        test_assert(m == 0);
        test_assert(prob == 1.0);
        test_assert(aadd_equivalent(q, qPM, 3, false, false));
        test_assert(aadd_equivalent(q, qPM, 3, true, false));
        test_assert(q == qPM);

        // |010>
        x3[2]=0; x3[1]=1; x3[0]=0;
        q   = qmdd_create_basis_state(3, x3);
        qPM = qmdd_measure_qubit(q, 0, 3, &m, &prob);
        test_assert(m == 0);
        test_assert(prob == 1.0);
        test_assert(aadd_equivalent(q, qPM, 3, false, false));
        test_assert(aadd_equivalent(q, qPM, 3, true, false));
        test_assert(q == qPM);

        // |011>
        x3[2]=0; x3[1]=1; x3[0]=1;
        q   = qmdd_create_basis_state(3, x3);
        qPM = qmdd_measure_qubit(q, 0, 3, &m, &prob);
        test_assert(m == 1);
        test_assert(prob == 0.0);
        test_assert(aadd_equivalent(q, qPM, 3, false, false));
        test_assert(aadd_equivalent(q, qPM, 3, true, false));
        test_assert(q == qPM);

        // |00+>
        x3[2]=0; x3[1]=0; x3[0]=0;
        q   = qmdd_create_basis_state(3, x3);
        q   = qmdd_gate(q, GATEID_H, 0);
        qPM = qmdd_measure_qubit(q, 0, 3, &m, &prob);
        test_assert(flt_abs(prob - 0.5) < cmap_get_tolerance());
        x3[2]=0; x3[1]=0; x3[0]=m; // either |000> or |001> depending on m
        q = qmdd_create_basis_state(3, x3); 
        test_assert(aadd_equivalent(q, qPM, 3, false, true));
        test_assert(aadd_equivalent(q, qPM, 3, true, false));
        test_assert(q == qPM);

        // |0+0>
        x3[2]=0; x3[1]=0; x3[0]=0;
        q   = qmdd_create_basis_state(3, x3);
        q   = qmdd_gate(q, GATEID_H, 1);
        qPM = qmdd_measure_qubit(q, 1, 3, &m, &prob);
        test_assert(flt_abs(prob - 0.5) < cmap_get_tolerance());
        x3[2]=0; x3[1]=m; x3[0]=0; // either |000> or |010> depending on m
        q = qmdd_create_basis_state(3, x3);
        test_assert(aadd_equivalent(q, qPM, 3, false, true));
        test_assert(aadd_equivalent(q, qPM, 3, true, false));
        test_assert(q == qPM);

        // |+00>
        x3[2]=0; x3[1]=0; x3[0]=0;
        q   = qmdd_create_basis_state(3, x3);
        q   = qmdd_gate(q, GATEID_H, 2);
        qPM = qmdd_measure_qubit(q, 2, 3, &m, &prob);
        test_assert(flt_abs(prob - 0.5) < cmap_get_tolerance());
        x3[2]=m; x3[1]=0; x3[0]=0; // either |000> or |100> depending on m
        q = qmdd_create_basis_state(3, x3);
        test_assert(aadd_equivalent(q, qPM, 3, false, true));
        test_assert(aadd_equivalent(q, qPM, 3, true, false));
        test_assert(q == qPM);

        // |00->
        x3[2]=0; x3[1]=0; x3[0]=1;
        q   = qmdd_create_basis_state(3, x3);
        q   = qmdd_gate(q, GATEID_H, 0);
        qPM = qmdd_measure_qubit(q, 0, 3, &m, &prob);
        test_assert(flt_abs(prob - 0.5) < cmap_get_tolerance());
        x3[2]=0; x3[1]=0; x3[0]=m; // either |000> or |001> depending on m
        q = qmdd_create_basis_state(3, x3); 
        test_assert(aadd_equivalent(q, qPM, 3, false, true));
        test_assert(aadd_equivalent(q, qPM, 3, true, false));
        test_assert(q == qPM);
        
        // |+++>, measure q0
        x3[2]=0; x3[1]=0; x3[0]=0;
        q   = qmdd_create_basis_state(3, x3);
        q   = qmdd_gate(q, GATEID_H, 0);
        q   = qmdd_gate(q, GATEID_H, 1);
        q   = qmdd_gate(q, GATEID_H, 2);
        qPM = qmdd_measure_qubit(q, 0, 3, &m, &prob); 
        test_assert(flt_abs(prob - 0.5) < cmap_get_tolerance());
        x3[2]=0; x3[1]=0; x3[0]=m; // either |++0> or |++1> depending on m
        q = qmdd_create_basis_state(3, x3); 
        q = qmdd_gate(q, GATEID_H, 1);
        q = qmdd_gate(q, GATEID_H, 2);
        test_assert(aadd_equivalent(q, qPM, 3, false, true));
        test_assert(aadd_equivalent(q, qPM, 3, true, false));
        test_assert(q == qPM); 

        // |+++>, measure q1
        x3[2]=0; x3[1]=0; x3[0]=0;
        q   = qmdd_create_basis_state(3, x3);
        q   = qmdd_gate(q, GATEID_H, 0);
        q   = qmdd_gate(q, GATEID_H, 1);
        q   = qmdd_gate(q, GATEID_H, 2);
        qPM = qmdd_measure_qubit(q, 1, 3, &m, &prob);
        test_assert(flt_abs(prob - 0.5) < cmap_get_tolerance());
        x3[2]=0; x3[1]=m; x3[0]=0; // either |+0+> or |+1+> depending on m
        q = qmdd_create_basis_state(3, x3); 
        q = qmdd_gate(q, GATEID_H, 0);
        q = qmdd_gate(q, GATEID_H, 2);
        test_assert(aadd_equivalent(q, qPM, 3, false, true));
        test_assert(aadd_equivalent(q, qPM, 3, true, false));
        test_assert(q == qPM);

        // [1/2, 1/2, 1/2, -1/2] = 1/2(|00> + |01> + |10> - |11>)_(q1, q0)
        x2[1]=0; x2[0]=0;
        q   = qmdd_create_basis_state(2, x2);
        q   = qmdd_gate(q, GATEID_H, 0);
        q   = qmdd_gate(q, GATEID_H, 1);
        q   = qmdd_cgate(q,GATEID_Z, 0, 1);
        x2[1]=0; x2[0]=0; a = aadd_getvalue(q, x2); test_assert(a == weight_lookup(cmake(0.5,0)));
        x2[1]=0; x2[0]=1; a = aadd_getvalue(q, x2); test_assert(a == weight_lookup(cmake(0.5,0)));
        x2[1]=1; x2[0]=0; a = aadd_getvalue(q, x2); test_assert(a == weight_lookup(cmake(0.5,0)));
        x2[1]=1; x2[0]=1; a = aadd_getvalue(q, x2); test_assert(a == weight_lookup(cmake(-0.5,0)));
        qPM = qmdd_measure_qubit(q, 0, 2, &m, &prob);
        test_assert(flt_abs(prob - 0.5) < cmap_get_tolerance());
        if (m == 0) { // expect 1/sqrt(2)(|00> + |10>)
            x2[1]=0; x2[0]=0; a = aadd_getvalue(qPM, x2); test_assert(a == weight_lookup(cmake(1.0/flt_sqrt(2.0),0)));
            x2[1]=0; x2[0]=1; a = aadd_getvalue(qPM, x2); test_assert(a == AADD_ZERO);
            x2[1]=1; x2[0]=0; a = aadd_getvalue(qPM, x2); test_assert(a == weight_lookup(cmake(1.0/flt_sqrt(2.0),0)));
            x2[1]=1; x2[0]=1; a = aadd_getvalue(qPM, x2); test_assert(a == AADD_ZERO);
        }
        if (m == 1) { // expect 1/sqrt(2)(|01> - |11>)
            x2[1]=0; x2[0]=0; a = aadd_getvalue(qPM, x2); test_assert(a == AADD_ZERO);
            x2[1]=0; x2[0]=1; a = aadd_getvalue(qPM, x2); test_assert(a == weight_lookup(cmake(1.0/flt_sqrt(2.0),0)));
            x2[1]=1; x2[0]=0; a = aadd_getvalue(qPM, x2); test_assert(a == AADD_ZERO);
            x2[1]=1; x2[0]=1; a = aadd_getvalue(qPM, x2); test_assert(a == weight_lookup(cmake(-1.0/flt_sqrt(2.0),0)));
        }
    }

    // Test measure all
    bool ms[3] = {0};
    int m_zer[3] = {0};
    for(int i=0; i < repeat; i++) {
        
        // |000>  
        x3[2]=0; x3[1]=0; x3[0]=0;
        q   = qmdd_create_basis_state(3, x3);
        qPM = qmdd_measure_all(q, 3, ms, &prob);
        test_assert(ms[0] == 0);
        test_assert(ms[1] == 0);
        test_assert(ms[2] == 0);
        test_assert(prob == 1.0);
        test_assert(aadd_equivalent(q, qPM, 3, false, false));
        test_assert(aadd_equivalent(q, qPM, 3, true, false));
        test_assert(q == qPM);

        // |010>
         x3[2]=0; x3[1]=1; x3[0]=0;
        q   = qmdd_create_basis_state(3, x3);
        qPM = qmdd_measure_all(q, 3, ms, &prob);
        test_assert(ms[0] == 0);
        test_assert(ms[1] == 1);
        test_assert(ms[2] == 0);
        test_assert(prob == 1.0);
        test_assert(aadd_equivalent(q, qPM, 3, false, false));
        test_assert(aadd_equivalent(q, qPM, 3, true, false));
        test_assert(q == qPM);

        // |011>
        x3[2]=0; x3[1]=1; x3[0]=1;
        q   = qmdd_create_basis_state(3, x3);
        qPM = qmdd_measure_all(q, 3, ms, &prob);
        test_assert(ms[0] == 1);
        test_assert(ms[1] == 1);
        test_assert(ms[2] == 0);
        test_assert(prob == 1.0);
        test_assert(aadd_equivalent(q, qPM, 3, false, false));
        test_assert(aadd_equivalent(q, qPM, 3, true, false));
        test_assert(q == qPM);

       // |00+>
        x3[2]=0; x3[1]=0; x3[0]=0;
        q   = qmdd_create_basis_state(3, x3);
        q   = qmdd_gate(q, GATEID_H, 0);
        qPM = qmdd_measure_all(q, 3, ms, &prob);
        test_assert(flt_abs(prob - 0.5) < cmap_get_tolerance());
        x3[2]=0; x3[1]=0; x3[0]=ms[0]; // either |000> or |001> depending on m
        q = qmdd_create_basis_state(3, x3); 
        test_assert(aadd_equivalent(q, qPM, 3, false, true));
        test_assert(aadd_equivalent(q, qPM, 3, true, false));
        test_assert(q == qPM);
        if (ms[0] == 0) m_zer[0] += 1;

        // |0+0>
        x3[2]=0; x3[1]=0; x3[0]=0;
        q   = qmdd_create_basis_state(3, x3);
        q   = qmdd_gate(q, GATEID_H, 1);
        qPM = qmdd_measure_all(q, 3, ms, &prob);
        test_assert(flt_abs(prob - 0.5) < cmap_get_tolerance());
        x3[2]=0; x3[1]=ms[1]; x3[0]=0; // either |000> or |010> depending on m
        q = qmdd_create_basis_state(3, x3);
        test_assert(aadd_equivalent(q, qPM, 3, false, true));
        test_assert(aadd_equivalent(q, qPM, 3, true, false));
        test_assert(q == qPM);
        if (ms[1] == 0) m_zer[1] += 1;

        // |+00>
        x3[2]=0; x3[1]=0; x3[0]=0;
        q   = qmdd_create_basis_state(3, x3);
        q   = qmdd_gate(q, GATEID_H, 2);
        qPM = qmdd_measure_all(q, 3, ms, &prob);
        test_assert(flt_abs(prob - 0.5) < cmap_get_tolerance());
        x3[2]=ms[2]; x3[1]=0; x3[0]=0; // either |000> or |100> depending on m
        q = qmdd_create_basis_state(3, x3);
        test_assert(aadd_equivalent(q, qPM, 3, false, true));
        test_assert(aadd_equivalent(q, qPM, 3, true, false));
        test_assert(q == qPM);
        if (ms[2] == 0) m_zer[2] += 1;
    }
    
    // TODO: more tests

    if(VERBOSE) printf("qmdd measurements:         ok\n");
    return 0;
}

int test_QFT()
{
    QMDD q3, q5, qref3, qref5;
    AMP a;
    complex_t c, cref;

    // 3 qubit QFT
    bool x3[] = {0,1,1}; // little endian (q0, q1, q2)
    q3 = qmdd_create_basis_state(3, x3);
    qref3 = qmdd_create_basis_state(3, x3);
    q3 = qmdd_circuit(q3, CIRCID_QFT, 0, 2);
    q3 = qmdd_circuit(q3, CIRCID_reverse_range, 0, 2);

    // check approx equal against output from qiskit
    x3[2]=0; x3[1]=0; x3[0]=0; a = aadd_getvalue(q3, x3); weight_value(a, &c); cref = cmake(3.5355339059327384e-01,-8.6595605623549353e-17);  test_assert(weight_approx_eq(&c, &cref));
    x3[2]=0; x3[1]=0; x3[0]=1; a = aadd_getvalue(q3, x3); weight_value(a, &c); cref = cmake(-3.5355339059327384e-01,8.6595605623549353e-17);  test_assert(weight_approx_eq(&c, &cref));
    x3[2]=0; x3[1]=1; x3[0]=0; a = aadd_getvalue(q3, x3); weight_value(a, &c); cref = cmake(-1.0824450702943669e-16,-3.5355339059327384e-01); test_assert(weight_approx_eq(&c, &cref));
    x3[2]=0; x3[1]=1; x3[0]=1; a = aadd_getvalue(q3, x3); weight_value(a, &c); cref = cmake(1.0824450702943669e-16,3.5355339059327384e-01);   test_assert(weight_approx_eq(&c, &cref));
    x3[2]=1; x3[1]=0; x3[0]=0; a = aadd_getvalue(q3, x3); weight_value(a, &c); cref = cmake(-2.5000000000000000e-01,2.5000000000000017e-01);  test_assert(weight_approx_eq(&c, &cref));
    x3[2]=1; x3[1]=0; x3[0]=1; a = aadd_getvalue(q3, x3); weight_value(a, &c); cref = cmake(2.5000000000000000e-01,-2.5000000000000017e-01);  test_assert(weight_approx_eq(&c, &cref));
    x3[2]=1; x3[1]=1; x3[0]=0; a = aadd_getvalue(q3, x3); weight_value(a, &c); cref = cmake(2.5000000000000017e-01,2.5000000000000000e-01);   test_assert(weight_approx_eq(&c, &cref));
    x3[2]=1; x3[1]=1; x3[0]=1; a = aadd_getvalue(q3, x3); weight_value(a, &c); cref = cmake(-2.5000000000000017e-01,-2.5000000000000000e-01); test_assert(weight_approx_eq(&c, &cref));
    test_assert(qmdd_is_unitvector(q3, 3));
    test_assert(aadd_is_ordered(q3, 3));

    // inverse QFT
    q3 = qmdd_circuit(q3, CIRCID_reverse_range, 0, 2);
    q3 = qmdd_circuit(q3, CIRCID_QFT_inv, 0, 2);
    test_assert(aadd_equivalent(q3, qref3, 3, false, false));
    test_assert(aadd_equivalent(q3, qref3, 3, true, false));
    test_assert(q3 == qref3);

    // 5 qubit QFT
    bool x5[] = {0,1,1,0,1};
    q5 = qmdd_create_basis_state(5, x5);
    qref5 = qmdd_create_basis_state(5, x5);
    q5 = qmdd_circuit(q5, CIRCID_QFT, 0, 4);
    q5 = qmdd_circuit(q5, CIRCID_reverse_range, 0, 4);
    test_assert(aadd_is_ordered(q5, 5));

    // check approx equal against output from qiskit
    x5[4]=0; x5[3]=0; x5[2]=0; x5[1]=0; x5[0]=0; a = aadd_getvalue(q5, x5); weight_value(a, &c); cref = cmake(1.7677669529663692e-01,-6.4946704217662027e-17);  test_assert(weight_approx_eq(&c, &cref));
    x5[4]=0; x5[3]=0; x5[2]=0; x5[1]=0; x5[0]=1; a = aadd_getvalue(q5, x5); weight_value(a, &c); cref = cmake(-1.7677669529663692e-01,6.4946704217662027e-17);  test_assert(weight_approx_eq(&c, &cref));
    x5[4]=0; x5[3]=0; x5[2]=0; x5[1]=1; x5[0]=0; a = aadd_getvalue(q5, x5); weight_value(a, &c); cref = cmake(7.5771154920605696e-17,1.7677669529663692e-01);   test_assert(weight_approx_eq(&c, &cref));
    x5[4]=0; x5[3]=0; x5[2]=0; x5[1]=1; x5[0]=1; a = aadd_getvalue(q5, x5); weight_value(a, &c); cref = cmake(-7.5771154920605696e-17,-1.7677669529663692e-01); test_assert(weight_approx_eq(&c, &cref));
    x5[4]=0; x5[3]=0; x5[2]=1; x5[1]=0; x5[0]=0; a = aadd_getvalue(q5, x5); weight_value(a, &c); cref = cmake(-1.2500000000000011e-01,-1.2499999999999999e-01); test_assert(weight_approx_eq(&c, &cref));
    x5[4]=0; x5[3]=0; x5[2]=1; x5[1]=0; x5[0]=1; a = aadd_getvalue(q5, x5); weight_value(a, &c); cref = cmake(1.2500000000000011e-01,1.2499999999999999e-01);   test_assert(weight_approx_eq(&c, &cref));
    x5[4]=0; x5[3]=0; x5[2]=1; x5[1]=1; x5[0]=0; a = aadd_getvalue(q5, x5); weight_value(a, &c); cref = cmake(1.2499999999999999e-01,-1.2500000000000011e-01);  test_assert(weight_approx_eq(&c, &cref));
    x5[4]=0; x5[3]=0; x5[2]=1; x5[1]=1; x5[0]=1; a = aadd_getvalue(q5, x5); weight_value(a, &c); cref = cmake(-1.2499999999999999e-01,1.2500000000000011e-01);  test_assert(weight_approx_eq(&c, &cref));
    x5[4]=0; x5[3]=1; x5[2]=0; x5[1]=0; x5[0]=0; a = aadd_getvalue(q5, x5); weight_value(a, &c); cref = cmake(6.7649512518274585e-02,-1.6332037060954713e-01);  test_assert(weight_approx_eq(&c, &cref));
    x5[4]=0; x5[3]=1; x5[2]=0; x5[1]=0; x5[0]=1; a = aadd_getvalue(q5, x5); weight_value(a, &c); cref = cmake(-6.7649512518274585e-02,1.6332037060954713e-01);  test_assert(weight_approx_eq(&c, &cref));
    x5[4]=0; x5[3]=1; x5[2]=0; x5[1]=1; x5[0]=0; a = aadd_getvalue(q5, x5); weight_value(a, &c); cref = cmake(1.6332037060954713e-01,6.7649512518274571e-02);   test_assert(weight_approx_eq(&c, &cref));
    x5[4]=0; x5[3]=1; x5[2]=0; x5[1]=1; x5[0]=1; a = aadd_getvalue(q5, x5); weight_value(a, &c); cref = cmake(-1.6332037060954713e-01,-6.7649512518274571e-02); test_assert(weight_approx_eq(&c, &cref));
    x5[4]=0; x5[3]=1; x5[2]=1; x5[1]=0; x5[0]=0; a = aadd_getvalue(q5, x5); weight_value(a, &c); cref = cmake(-1.6332037060954710e-01,6.7649512518274696e-02);  test_assert(weight_approx_eq(&c, &cref));
    x5[4]=0; x5[3]=1; x5[2]=1; x5[1]=0; x5[0]=1; a = aadd_getvalue(q5, x5); weight_value(a, &c); cref = cmake(1.6332037060954710e-01,-6.7649512518274696e-02);  test_assert(weight_approx_eq(&c, &cref));
    x5[4]=0; x5[3]=1; x5[2]=1; x5[1]=1; x5[0]=0; a = aadd_getvalue(q5, x5); weight_value(a, &c); cref = cmake(-6.7649512518274710e-02,-1.6332037060954710e-01); test_assert(weight_approx_eq(&c, &cref));
    x5[4]=0; x5[3]=1; x5[2]=1; x5[1]=1; x5[0]=1; a = aadd_getvalue(q5, x5); weight_value(a, &c); cref = cmake(6.7649512518274710e-02,1.6332037060954710e-01);   test_assert(weight_approx_eq(&c, &cref));
    x5[4]=1; x5[3]=0; x5[2]=0; x5[1]=0; x5[0]=0; a = aadd_getvalue(q5, x5); weight_value(a, &c); cref = cmake(-1.4698445030241986e-01,9.8211869798387877e-02);  test_assert(weight_approx_eq(&c, &cref));
    x5[4]=1; x5[3]=0; x5[2]=0; x5[1]=0; x5[0]=1; a = aadd_getvalue(q5, x5); weight_value(a, &c); cref = cmake(1.4698445030241986e-01,-9.8211869798387877e-02);  test_assert(weight_approx_eq(&c, &cref));
    x5[4]=1; x5[3]=0; x5[2]=0; x5[1]=1; x5[0]=0; a = aadd_getvalue(q5, x5); weight_value(a, &c); cref = cmake(-9.8211869798387877e-02,-1.4698445030241986e-01); test_assert(weight_approx_eq(&c, &cref));
    x5[4]=1; x5[3]=0; x5[2]=0; x5[1]=1; x5[0]=1; a = aadd_getvalue(q5, x5); weight_value(a, &c); cref = cmake(9.8211869798387877e-02,1.4698445030241986e-01);   test_assert(weight_approx_eq(&c, &cref));
    x5[4]=1; x5[3]=0; x5[2]=1; x5[1]=0; x5[0]=0; a = aadd_getvalue(q5, x5); weight_value(a, &c); cref = cmake(1.7337998066526852e-01,3.4487422410367806e-02);   test_assert(weight_approx_eq(&c, &cref));
    x5[4]=1; x5[3]=0; x5[2]=1; x5[1]=0; x5[0]=1; a = aadd_getvalue(q5, x5); weight_value(a, &c); cref = cmake(-1.7337998066526852e-01,-3.4487422410367806e-02); test_assert(weight_approx_eq(&c, &cref));
    x5[4]=1; x5[3]=0; x5[2]=1; x5[1]=1; x5[0]=0; a = aadd_getvalue(q5, x5); weight_value(a, &c); cref = cmake(-3.4487422410367799e-02,1.7337998066526852e-01);  test_assert(weight_approx_eq(&c, &cref));
    x5[4]=1; x5[3]=0; x5[2]=1; x5[1]=1; x5[0]=1; a = aadd_getvalue(q5, x5); weight_value(a, &c); cref = cmake(3.4487422410367799e-02,-1.7337998066526852e-01);  test_assert(weight_approx_eq(&c, &cref));
    x5[4]=1; x5[3]=1; x5[2]=0; x5[1]=0; x5[0]=0; a = aadd_getvalue(q5, x5); weight_value(a, &c); cref = cmake(3.4487422410367972e-02,1.7337998066526850e-01);   test_assert(weight_approx_eq(&c, &cref));
    x5[4]=1; x5[3]=1; x5[2]=0; x5[1]=0; x5[0]=1; a = aadd_getvalue(q5, x5); weight_value(a, &c); cref = cmake(-3.4487422410367972e-02,-1.7337998066526850e-01); test_assert(weight_approx_eq(&c, &cref));
    x5[4]=1; x5[3]=1; x5[2]=0; x5[1]=1; x5[0]=0; a = aadd_getvalue(q5, x5); weight_value(a, &c); cref = cmake(-1.7337998066526850e-01,3.4487422410367986e-02);  test_assert(weight_approx_eq(&c, &cref));
    x5[4]=1; x5[3]=1; x5[2]=0; x5[1]=1; x5[0]=1; a = aadd_getvalue(q5, x5); weight_value(a, &c); cref = cmake(1.7337998066526850e-01,-3.4487422410367986e-02);  test_assert(weight_approx_eq(&c, &cref));
    x5[4]=1; x5[3]=1; x5[2]=1; x5[1]=0; x5[0]=0; a = aadd_getvalue(q5, x5); weight_value(a, &c); cref = cmake(9.8211869798387752e-02,-1.4698445030241994e-01);  test_assert(weight_approx_eq(&c, &cref));
    x5[4]=1; x5[3]=1; x5[2]=1; x5[1]=0; x5[0]=1; a = aadd_getvalue(q5, x5); weight_value(a, &c); cref = cmake(-9.8211869798387752e-02,1.4698445030241994e-01);  test_assert(weight_approx_eq(&c, &cref));
    x5[4]=1; x5[3]=1; x5[2]=1; x5[1]=1; x5[0]=0; a = aadd_getvalue(q5, x5); weight_value(a, &c); cref = cmake(1.4698445030241994e-01,9.8211869798387752e-02);   test_assert(weight_approx_eq(&c, &cref));
    x5[4]=1; x5[3]=1; x5[2]=1; x5[1]=1; x5[0]=1; a = aadd_getvalue(q5, x5); weight_value(a, &c); cref = cmake(-1.4698445030241994e-01,-9.8211869798387752e-02); test_assert(weight_approx_eq(&c, &cref));
    test_assert(qmdd_is_unitvector(q5, 5));

    // inverse QFT
    q5 = qmdd_circuit(q5, CIRCID_reverse_range, 0, 4);
    q5 = qmdd_circuit(q5, CIRCID_QFT_inv, 0, 4);
    test_assert(aadd_equivalent(q5, qref5, 5, false, false));
    test_assert(aadd_equivalent(q5, qref5, 5, true, false));
    test_assert(q5 == qref5);
    
    
    if(VERBOSE) printf("qmdd QFT:                  ok\n");
    return 0;
}

int test_5qubit_circuit()
{
    QMDD q, qref;
    uint64_t node_count;
    int n_qubits = 5;
    bool x5[] = {0,0,0,0,0};
    FILE *fp;
    
    // 5 qubit state
    qref = qmdd_create_basis_state(n_qubits, x5);
    q    = qmdd_create_basis_state(n_qubits, x5);

    // 32 gates 
    q = qmdd_cgate(q, GATEID_Z, 1, 2);       q = qmdd_gate(q, GATEID_X, 2);       test_assert(qmdd_is_unitvector(q, 5));
    q = qmdd_cgate(q, GATEID_Z, 3, 4);       q = qmdd_cgate(q, GATEID_X, 1, 3);   test_assert(qmdd_is_unitvector(q, 5));
    q = qmdd_gate(q, GATEID_Z, 1);           q = qmdd_gate(q, GATEID_H, 4);       test_assert(qmdd_is_unitvector(q, 5));
    q = qmdd_gate(q, GATEID_H, 0);           q = qmdd_cgate(q, GATEID_X, 1, 3);   test_assert(qmdd_is_unitvector(q, 5));
    q = qmdd_gate(q, GATEID_H, 3);           q = qmdd_gate(q, GATEID_H, 0);       test_assert(qmdd_is_unitvector(q, 5));
    q = qmdd_gate(q, GATEID_Z, 1);           q = qmdd_cgate(q, GATEID_X, 1, 2);   test_assert(qmdd_is_unitvector(q, 5));
    q = qmdd_gate(q, GATEID_X, 1);           q = qmdd_cgate(q, GATEID_X, 0, 4);   test_assert(qmdd_is_unitvector(q, 5));
    q = qmdd_gate(q, GATEID_H, 4);           q = qmdd_cgate(q, GATEID_X, 0, 1);   test_assert(qmdd_is_unitvector(q, 5));
    q = qmdd_cgate(q, GATEID_Z, 0, 4);       q = qmdd_cgate(q, GATEID_Z, 0, 4);   test_assert(qmdd_is_unitvector(q, 5));
    q = qmdd_cgate(q, GATEID_Z, 2, 3);       q = qmdd_gate(q, GATEID_Z, 0);       test_assert(qmdd_is_unitvector(q, 5));
    q = qmdd_cgate(q, GATEID_X, 3, 4);       q = qmdd_cgate(q, GATEID_Z, 0, 2);   test_assert(qmdd_is_unitvector(q, 5));
    q = qmdd_gate(q, GATEID_Z, 3);           q = qmdd_cgate(q, GATEID_Z, 1, 3);   test_assert(qmdd_is_unitvector(q, 5));
    q = qmdd_cgate(q, GATEID_X, 0, 3);       q = qmdd_cgate(q, GATEID_Z, 0, 4);   test_assert(qmdd_is_unitvector(q, 5));
    q = qmdd_cgate(q, GATEID_Z, 0, 4);       q = qmdd_gate(q, GATEID_H, 1);       test_assert(qmdd_is_unitvector(q, 5));
    q = qmdd_gate(q, GATEID_H, 3);           q = qmdd_gate(q, GATEID_X, 1);       test_assert(qmdd_is_unitvector(q, 5));
    q = qmdd_cgate(q, GATEID_Z, 1, 2);       q = qmdd_gate(q, GATEID_Z, 1);       test_assert(qmdd_is_unitvector(q, 5));
    node_count = aadd_countnodes(q);
    if (test_measure_random_state(q, n_qubits)) return 1;

    // inverse
    q = qmdd_gate(q, GATEID_Z, 1);           q = qmdd_cgate(q, GATEID_Z, 1, 2);   test_assert(qmdd_is_unitvector(q, 5));
    q = qmdd_gate(q, GATEID_X, 1);           q = qmdd_gate(q, GATEID_H, 3);       test_assert(qmdd_is_unitvector(q, 5));
    q = qmdd_gate(q, GATEID_H, 1);           q = qmdd_cgate(q, GATEID_Z, 0, 4);   test_assert(qmdd_is_unitvector(q, 5));
    q = qmdd_cgate(q, GATEID_Z, 0, 4);       q = qmdd_cgate(q, GATEID_X, 0, 3);   test_assert(qmdd_is_unitvector(q, 5));
    q = qmdd_cgate(q, GATEID_Z, 1, 3);       q = qmdd_gate(q, GATEID_Z, 3);       test_assert(qmdd_is_unitvector(q, 5));
    q = qmdd_cgate(q, GATEID_Z, 0, 2);       q = qmdd_cgate(q, GATEID_X, 3, 4);   test_assert(qmdd_is_unitvector(q, 5));
    q = qmdd_gate(q, GATEID_Z, 0);           q = qmdd_cgate(q, GATEID_Z, 2, 3);   test_assert(qmdd_is_unitvector(q, 5));
    q = qmdd_cgate(q, GATEID_Z, 0, 4);       q = qmdd_cgate(q, GATEID_Z, 0, 4);   test_assert(qmdd_is_unitvector(q, 5));
    q = qmdd_cgate(q, GATEID_X, 0, 1);       q = qmdd_gate(q, GATEID_H, 4);       test_assert(qmdd_is_unitvector(q, 5));
    q = qmdd_cgate(q, GATEID_X, 0, 4);       q = qmdd_gate(q, GATEID_X, 1);       test_assert(qmdd_is_unitvector(q, 5));
    q = qmdd_cgate(q, GATEID_X, 1, 2);       q = qmdd_gate(q, GATEID_Z, 1);       test_assert(qmdd_is_unitvector(q, 5));
    q = qmdd_gate(q, GATEID_H, 0);           q = qmdd_gate(q, GATEID_H, 3);       test_assert(qmdd_is_unitvector(q, 5));
    q = qmdd_cgate(q, GATEID_X, 1, 3);       q = qmdd_gate(q, GATEID_H, 0);       test_assert(qmdd_is_unitvector(q, 5));
    q = qmdd_gate(q, GATEID_H, 4);           q = qmdd_gate(q, GATEID_Z, 1);       test_assert(qmdd_is_unitvector(q, 5));
    q = qmdd_cgate(q, GATEID_X, 1, 3);       q = qmdd_cgate(q, GATEID_Z, 3, 4);   test_assert(qmdd_is_unitvector(q, 5));
    q = qmdd_gate(q, GATEID_X, 2);           q = qmdd_cgate(q, GATEID_Z, 1, 2);   test_assert(qmdd_is_unitvector(q, 5));

    test_assert(aadd_equivalent(q, qref, n_qubits, false, VERBOSE)); // check approx equiv
    test_assert(aadd_equivalent(q, qref, n_qubits, true,  VERBOSE)); // check exact equiv
    test_assert(q == qref);

    fp = fopen("5_qubit_res.dot", "w");
    aadd_fprintdot(fp, q, true);
    fclose(fp);

    if(VERBOSE) printf("qmdd 5 qubit circuit:      ok (%" PRIu64 " nodes)\n", node_count);
    return 0;
}

int test_10qubit_circuit()
{
    QMDD q, qref;
    uint64_t node_count;
    bool x10[] = {0,0,0,0,0,0,0,0,0,0};
    
    // 10 qubit state
    qref = qmdd_create_basis_state(10, x10);
    q    = qmdd_create_basis_state(10, x10);

    // 30 random* Clifford gates            *chosen by a fair dice roll
    q = qmdd_cgate(q, GATEID_X, 1, 3);       q = qmdd_gate(q, GATEID_H, 0);       test_assert(qmdd_is_unitvector(q, 10));
    q = qmdd_gate(q, GATEID_X, 6);           q = qmdd_cgate(q, GATEID_X, 6, 9);   test_assert(qmdd_is_unitvector(q, 10));
    q = qmdd_gate(q, GATEID_H, 4);           q = qmdd_cgate(q, GATEID_X, 3, 5);   test_assert(qmdd_is_unitvector(q, 10));
    q = qmdd_gate(q, GATEID_H, 1);           q = qmdd_gate(q, GATEID_X, 1);       test_assert(qmdd_is_unitvector(q, 10));
    q = qmdd_cgate(q, GATEID_X, 3, 8);       q = qmdd_cgate(q, GATEID_Z, 3, 6);   test_assert(qmdd_is_unitvector(q, 10));
    q = qmdd_gate(q, GATEID_Z, 3);           q = qmdd_cgate(q, GATEID_X, 0, 7);   test_assert(qmdd_is_unitvector(q, 10));
    q = qmdd_cgate(q, GATEID_X, 1, 9);       q = qmdd_gate(q, GATEID_H, 4);       test_assert(qmdd_is_unitvector(q, 10));
    q = qmdd_cgate(q, GATEID_X, 0, 2);       q = qmdd_gate(q, GATEID_X, 2);       test_assert(qmdd_is_unitvector(q, 10));
    q = qmdd_cgate(q, GATEID_X, 5, 8);       q = qmdd_cgate(q, GATEID_X, 0, 4);   test_assert(qmdd_is_unitvector(q, 10));
    q = qmdd_cgate(q, GATEID_X, 0, 8);       q = qmdd_cgate(q, GATEID_X, 6, 9);   test_assert(qmdd_is_unitvector(q, 10));
    q = qmdd_cgate(q, GATEID_X, 0, 9);       q = qmdd_gate(q, GATEID_X, 9);       test_assert(qmdd_is_unitvector(q, 10));
    q = qmdd_cgate(q, GATEID_X, 4, 9);       q = qmdd_cgate(q, GATEID_Z, 2, 7);   test_assert(qmdd_is_unitvector(q, 10));
    q = qmdd_cgate(q, GATEID_Z, 7, 8);       q = qmdd_gate(q, GATEID_X, 7);       test_assert(qmdd_is_unitvector(q, 10));
    q = qmdd_gate(q, GATEID_Z, 2);           q = qmdd_gate(q, GATEID_Z, 7);       test_assert(qmdd_is_unitvector(q, 10));
    q = qmdd_gate(q, GATEID_X, 6);           q = qmdd_gate(q, GATEID_X, 1);       test_assert(qmdd_is_unitvector(q, 10));
    node_count = aadd_countnodes(q);
    if (test_measure_random_state(q, 10)) return 1;

    // inverse
    q = qmdd_gate(q, GATEID_X, 1);           q = qmdd_gate(q, GATEID_X, 6);       test_assert(qmdd_is_unitvector(q, 10));
    q = qmdd_gate(q, GATEID_Z, 7);           q = qmdd_gate(q, GATEID_Z, 2);       test_assert(qmdd_is_unitvector(q, 10));
    q = qmdd_gate(q, GATEID_X, 7);           q = qmdd_cgate(q, GATEID_Z, 7, 8);   test_assert(qmdd_is_unitvector(q, 10));
    q = qmdd_cgate(q, GATEID_Z, 2, 7);       q = qmdd_cgate(q, GATEID_X, 4, 9);   test_assert(qmdd_is_unitvector(q, 10));
    q = qmdd_gate(q, GATEID_X, 9);           q = qmdd_cgate(q, GATEID_X, 0, 9);   test_assert(qmdd_is_unitvector(q, 10));
    q = qmdd_cgate(q, GATEID_X, 6, 9);       q = qmdd_cgate(q, GATEID_X, 0, 8);   test_assert(qmdd_is_unitvector(q, 10));
    q = qmdd_cgate(q, GATEID_X, 0, 4);       q = qmdd_cgate(q, GATEID_X, 5, 8);   test_assert(qmdd_is_unitvector(q, 10));
    q = qmdd_gate(q, GATEID_X, 2);           q = qmdd_cgate(q, GATEID_X, 0, 2);   test_assert(qmdd_is_unitvector(q, 10));
    q = qmdd_gate(q, GATEID_H, 4);           q = qmdd_cgate(q, GATEID_X, 1, 9);   test_assert(qmdd_is_unitvector(q, 10));
    q = qmdd_cgate(q, GATEID_X, 0, 7);       q = qmdd_gate(q, GATEID_Z, 3);       test_assert(qmdd_is_unitvector(q, 10));
    q = qmdd_cgate(q, GATEID_Z, 3, 6);       q = qmdd_cgate(q, GATEID_X, 3, 8);   test_assert(qmdd_is_unitvector(q, 10));
    q = qmdd_gate(q, GATEID_X, 1);           q = qmdd_gate(q, GATEID_H, 1);       test_assert(qmdd_is_unitvector(q, 10));
    q = qmdd_cgate(q, GATEID_X, 3, 5);       q = qmdd_gate(q, GATEID_H, 4);       test_assert(qmdd_is_unitvector(q, 10));
    q = qmdd_cgate(q, GATEID_X, 6, 9);       q = qmdd_gate(q, GATEID_X, 6);       test_assert(qmdd_is_unitvector(q, 10));
    q = qmdd_gate(q, GATEID_H, 0);           q = qmdd_cgate(q, GATEID_X, 1, 3);   test_assert(qmdd_is_unitvector(q, 10));

    test_assert(aadd_equivalent(q, qref, 10, false, VERBOSE)); // check approx equiv
    test_assert(aadd_equivalent(q, qref, 10, true,  VERBOSE)); // check exact equiv
    test_assert(q == qref);

    if(VERBOSE) printf("qmdd 10 qubit circuit:     ok (%" PRIu64 " nodes)\n", node_count);
    return 0;
}

int test_20qubit_circuit()
{
    QMDD q, qref;
    uint64_t node_count;
    bool x20[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

    // 20 qubit state
    qref = qmdd_create_basis_state(20, x20);
    q    = qmdd_create_basis_state(20, x20);

    // 100 gates
    q = qmdd_cgate(q, GATEID_Z, 4, 18);      q = qmdd_gate(q, GATEID_H, 16);          q = qmdd_cgate(q, GATEID_X, 1, 12);      q = qmdd_gate(q, GATEID_Z, 4);
    q = qmdd_cgate(q, GATEID_Z, 9, 10);      q = qmdd_cgate(q, GATEID_Z, 0, 4);       q = qmdd_cgate(q, GATEID_Z, 1, 16);      q = qmdd_cgate(q, GATEID_Z, 13, 16);
    q = qmdd_cgate(q, GATEID_Z, 7, 11);      q = qmdd_cgate(q, GATEID_X, 3, 5);       q = qmdd_cgate(q, GATEID_Z, 1, 4);       q = qmdd_cgate(q, GATEID_X, 6, 16);
    q = qmdd_cgate(q, GATEID_X, 3, 18);      q = qmdd_cgate(q, GATEID_X, 2, 15);      q = qmdd_cgate(q, GATEID_X, 7, 10);      q = qmdd_gate(q, GATEID_Z, 6);
    q = qmdd_cgate(q, GATEID_X, 3, 6);       q = qmdd_cgate(q, GATEID_Z, 11, 16);     q = qmdd_cgate(q, GATEID_X, 5, 19);      q = qmdd_gate(q, GATEID_Z, 18);
    q = qmdd_cgate(q, GATEID_Z, 14, 15);     q = qmdd_cgate(q, GATEID_Z, 10, 12);     q = qmdd_gate(q, GATEID_H, 8);           q = qmdd_gate(q, GATEID_X, 9);
    q = qmdd_gate(q, GATEID_X, 8);           q = qmdd_cgate(q, GATEID_X, 7, 18);      q = qmdd_gate(q, GATEID_X, 17);          q = qmdd_gate(q, GATEID_Z, 11);
    q = qmdd_cgate(q, GATEID_X, 12, 16);     q = qmdd_gate(q, GATEID_X, 18);          q = qmdd_gate(q, GATEID_Z, 4);           q = qmdd_gate(q, GATEID_X, 18);
    q = qmdd_cgate(q, GATEID_X, 4, 10);      q = qmdd_gate(q, GATEID_X, 15);          q = qmdd_cgate(q, GATEID_Z, 16, 18);     q = qmdd_cgate(q, GATEID_Z, 0, 15);
    q = qmdd_cgate(q, GATEID_X, 7, 10);      q = qmdd_gate(q, GATEID_X, 18);          q = qmdd_gate(q, GATEID_Z, 16);          q = qmdd_cgate(q, GATEID_X, 7, 18);
    q = qmdd_gate(q, GATEID_X, 16);          q = qmdd_gate(q, GATEID_X, 2);           q = qmdd_cgate(q, GATEID_X, 9, 10);      q = qmdd_gate(q, GATEID_X, 6);
    q = qmdd_gate(q, GATEID_X, 18);          q = qmdd_gate(q, GATEID_Z, 11);          q = qmdd_cgate(q, GATEID_Z, 4, 5);       q = qmdd_gate(q, GATEID_X, 1);
    q = qmdd_cgate(q, GATEID_Z, 2, 19);      q = qmdd_cgate(q, GATEID_X, 8, 9);       q = qmdd_cgate(q, GATEID_Z, 10, 12);     q = qmdd_cgate(q, GATEID_Z, 11, 16);
    q = qmdd_cgate(q, GATEID_X, 13, 19);     q = qmdd_cgate(q, GATEID_Z, 1, 3);       q = qmdd_gate(q, GATEID_X, 6);           q = qmdd_gate(q, GATEID_X, 15);
    q = qmdd_gate(q, GATEID_Z, 0);           q = qmdd_cgate(q, GATEID_X, 0, 15);      q = qmdd_gate(q, GATEID_H, 16);          q = qmdd_gate(q, GATEID_Z, 8);
    q = qmdd_cgate(q, GATEID_X, 12, 14);     q = qmdd_cgate(q, GATEID_Z, 2, 18);      q = qmdd_cgate(q, GATEID_X, 12, 15);     q = qmdd_gate(q, GATEID_X, 9);
    q = qmdd_gate(q, GATEID_Z, 12);          q = qmdd_gate(q, GATEID_X, 3);           q = qmdd_gate(q, GATEID_X, 0);           q = qmdd_cgate(q, GATEID_X, 1, 4);
    q = qmdd_gate(q, GATEID_H, 1);           q = qmdd_gate(q, GATEID_X, 19);          q = qmdd_gate(q, GATEID_X, 5);           q = qmdd_cgate(q, GATEID_Z, 2, 16);
    q = qmdd_gate(q, GATEID_X, 4);           q = qmdd_cgate(q, GATEID_X, 9, 11);      q = qmdd_cgate(q, GATEID_X, 0, 7);       q = qmdd_gate(q, GATEID_Z, 12);
    q = qmdd_cgate(q, GATEID_X, 9, 11);      q = qmdd_gate(q, GATEID_Z, 13);          q = qmdd_cgate(q, GATEID_X, 12, 16);     q = qmdd_gate(q, GATEID_Z, 10);
    q = qmdd_gate(q, GATEID_X, 4);           q = qmdd_gate(q, GATEID_Z, 16);          q = qmdd_cgate(q, GATEID_Z, 4, 17);      q = qmdd_gate(q, GATEID_Z, 7);
    q = qmdd_gate(q, GATEID_H, 4);           q = qmdd_cgate(q, GATEID_Z, 6, 7);       q = qmdd_cgate(q, GATEID_X, 12, 19);     q = qmdd_gate(q, GATEID_Z, 15);
    q = qmdd_cgate(q, GATEID_X, 5, 11);      q = qmdd_cgate(q, GATEID_X, 9, 17);      q = qmdd_gate(q, GATEID_Z, 3);           q = qmdd_cgate(q, GATEID_X, 11, 18);
    q = qmdd_cgate(q, GATEID_Z, 5, 15);      q = qmdd_cgate(q, GATEID_X, 0, 15);      q = qmdd_cgate(q, GATEID_X, 1, 6);       q = qmdd_cgate(q, GATEID_X, 8, 16);
    q = qmdd_cgate(q, GATEID_X, 5, 19);      q = qmdd_cgate(q, GATEID_Z, 3, 18);      q = qmdd_cgate(q, GATEID_X, 5, 8);       q = qmdd_cgate(q, GATEID_Z, 14, 18);
    node_count = aadd_countnodes(q);
    if (test_measure_random_state(q, 20)) return 1;

    // inverse
    q = qmdd_cgate(q, GATEID_Z, 14, 18);     q = qmdd_cgate(q, GATEID_X, 5, 8);       q = qmdd_cgate(q, GATEID_Z, 3, 18);      q = qmdd_cgate(q, GATEID_X, 5, 19);
    q = qmdd_cgate(q, GATEID_X, 8, 16);      q = qmdd_cgate(q, GATEID_X, 1, 6);       q = qmdd_cgate(q, GATEID_X, 0, 15);      q = qmdd_cgate(q, GATEID_Z, 5, 15);
    q = qmdd_cgate(q, GATEID_X, 11, 18);     q = qmdd_gate(q, GATEID_Z, 3);           q = qmdd_cgate(q, GATEID_X, 9, 17);      q = qmdd_cgate(q, GATEID_X, 5, 11);
    q = qmdd_gate(q, GATEID_Z, 15);          q = qmdd_cgate(q, GATEID_X, 12, 19);     q = qmdd_cgate(q, GATEID_Z, 6, 7);       q = qmdd_gate(q, GATEID_H, 4);
    q = qmdd_gate(q, GATEID_Z, 7);           q = qmdd_cgate(q, GATEID_Z, 4, 17);      q = qmdd_gate(q, GATEID_Z, 16);          q = qmdd_gate(q, GATEID_X, 4);
    q = qmdd_gate(q, GATEID_Z, 10);          q = qmdd_cgate(q, GATEID_X, 12, 16);     q = qmdd_gate(q, GATEID_Z, 13);          q = qmdd_cgate(q, GATEID_X, 9, 11);
    q = qmdd_gate(q, GATEID_Z, 12);          q = qmdd_cgate(q, GATEID_X, 0, 7);       q = qmdd_cgate(q, GATEID_X, 9, 11);      q = qmdd_gate(q, GATEID_X, 4);
    q = qmdd_cgate(q, GATEID_Z, 2, 16);      q = qmdd_gate(q, GATEID_X, 5);           q = qmdd_gate(q, GATEID_X, 19);          q = qmdd_gate(q, GATEID_H, 1);
    q = qmdd_cgate(q, GATEID_X, 1, 4);       q = qmdd_gate(q, GATEID_X, 0);           q = qmdd_gate(q, GATEID_X, 3);           q = qmdd_gate(q, GATEID_Z, 12);
    q = qmdd_gate(q, GATEID_X, 9);           q = qmdd_cgate(q, GATEID_X, 12, 15);     q = qmdd_cgate(q, GATEID_Z, 2, 18);      q = qmdd_cgate(q, GATEID_X, 12, 14);
    q = qmdd_gate(q, GATEID_Z, 8);           q = qmdd_gate(q, GATEID_H, 16);          q = qmdd_cgate(q, GATEID_X, 0, 15);      q = qmdd_gate(q, GATEID_Z, 0);
    q = qmdd_gate(q, GATEID_X, 15);          q = qmdd_gate(q, GATEID_X, 6);           q = qmdd_cgate(q, GATEID_Z, 1, 3);       q = qmdd_cgate(q, GATEID_X, 13, 19);
    q = qmdd_cgate(q, GATEID_Z, 11, 16);     q = qmdd_cgate(q, GATEID_Z, 10, 12);     q = qmdd_cgate(q, GATEID_X, 8, 9);       q = qmdd_cgate(q, GATEID_Z, 2, 19);
    q = qmdd_gate(q, GATEID_X, 1);           q = qmdd_cgate(q, GATEID_Z, 4, 5);       q = qmdd_gate(q, GATEID_Z, 11);          q = qmdd_gate(q, GATEID_X, 18);
    q = qmdd_gate(q, GATEID_X, 6);           q = qmdd_cgate(q, GATEID_X, 9, 10);      q = qmdd_gate(q, GATEID_X, 2);           q = qmdd_gate(q, GATEID_X, 16);
    q = qmdd_cgate(q, GATEID_X, 7, 18);      q = qmdd_gate(q, GATEID_Z, 16);          q = qmdd_gate(q, GATEID_X, 18);          q = qmdd_cgate(q, GATEID_X, 7, 10);
    q = qmdd_cgate(q, GATEID_Z, 0, 15);      q = qmdd_cgate(q, GATEID_Z, 16, 18);     q = qmdd_gate(q, GATEID_X, 15);          q = qmdd_cgate(q, GATEID_X, 4, 10);
    q = qmdd_gate(q, GATEID_X, 18);          q = qmdd_gate(q, GATEID_Z, 4);           q = qmdd_gate(q, GATEID_X, 18);          q = qmdd_cgate(q, GATEID_X, 12, 16);
    q = qmdd_gate(q, GATEID_Z, 11);          q = qmdd_gate(q, GATEID_X, 17);          q = qmdd_cgate(q, GATEID_X, 7, 18);      q = qmdd_gate(q, GATEID_X, 8);
    q = qmdd_gate(q, GATEID_X, 9);           q = qmdd_gate(q, GATEID_H, 8);           q = qmdd_cgate(q, GATEID_Z, 10, 12);     q = qmdd_cgate(q, GATEID_Z, 14, 15);
    q = qmdd_gate(q, GATEID_Z, 18);          q = qmdd_cgate(q, GATEID_X, 5, 19);      q = qmdd_cgate(q, GATEID_Z, 11, 16);     q = qmdd_cgate(q, GATEID_X, 3, 6);
    q = qmdd_gate(q, GATEID_Z, 6);           q = qmdd_cgate(q, GATEID_X, 7, 10);      q = qmdd_cgate(q, GATEID_X, 2, 15);      q = qmdd_cgate(q, GATEID_X, 3, 18);
    q = qmdd_cgate(q, GATEID_X, 6, 16);      q = qmdd_cgate(q, GATEID_Z, 1, 4);       q = qmdd_cgate(q, GATEID_X, 3, 5);       q = qmdd_cgate(q, GATEID_Z, 7, 11);
    q = qmdd_cgate(q, GATEID_Z, 13, 16);     q = qmdd_cgate(q, GATEID_Z, 1, 16);      q = qmdd_cgate(q, GATEID_Z, 0, 4);       q = qmdd_cgate(q, GATEID_Z, 9, 10);
    q = qmdd_gate(q, GATEID_Z, 4);           q = qmdd_cgate(q, GATEID_X, 1, 12);      q = qmdd_gate(q, GATEID_H, 16);          q = qmdd_cgate(q, GATEID_Z, 4, 18);

    test_assert(aadd_equivalent(q, qref, 20, false, VERBOSE)); // check approx equiv
    test_assert(aadd_equivalent(q, qref, 20, true,  VERBOSE)); // check exact equiv
    test_assert(q == qref);

    if(VERBOSE) printf("qmdd 20 qubit circuit:     ok (%" PRIu64 " nodes)\n", node_count);
    return 0;
}

int run_qmdd_tests()
{
    // we are not testing garbage collection
    sylvan_gc_disable();

    // circuits
    if (test_swap_circuit()) return 1;
    if (test_cswap_circuit()) return 1;
    if (test_tensor_product()) return 1;
    if (test_measurements()) return 1;
    if (test_5qubit_circuit()) return 1;
    if (test_10qubit_circuit()) return 1;
    //if (test_20qubit_circuit()) return 1;
    if (test_QFT()) return 1;

    return 0;
}

int test_with(bool it_ref, int amps_backend, int norm_strat) 
{
    // Standard Lace initialization
    int workers = 1;
    lace_start(workers, 0);
    printf("%d worker(s), ", workers);
    
    // Simple Sylvan initialization
    sylvan_set_sizes(1LL<<25, 1LL<<25, 1LL<<16, 1LL<<16);
    sylvan_init_package();
    qsylvan_init_simulator(1LL<<11, -1, amps_backend, norm_strat);
    qmdd_set_testing_mode(true); // turn on internal sanity tests
    qsylvan_set_iterative_refinement(it_ref);

    printf("it_ref = %d, amps backend = %d, norm strategy = %d:\n", it_ref, amps_backend, norm_strat);
    int res = run_qmdd_tests();

    sylvan_quit();
    lace_stop();
    return res;
}

int runtests()
{
    for (int it_ref = 0; it_ref <= 1; it_ref++) {
        for (int backend = 0; backend < n_backends; backend++) {
            for (int norm_strat = 0; norm_strat < n_norm_strategies; norm_strat++) {
                if (test_with((bool)it_ref, backend, norm_strat)) return 1;
            }
        }
    }
    return 0;
}

int main()
{
    return runtests();
}
