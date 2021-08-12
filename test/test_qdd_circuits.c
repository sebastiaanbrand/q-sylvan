#include <stdio.h>
#include <time.h>

#include "qsylvan.h"
#include "test_assert.h"

bool VERBOSE = true;

int test_swap_circuit()
{
    QDD q;
    bool x3[] = {0,0,0};
    AMP a;

    LACE_ME;

    q = qdd_create_basis_state(3, x3);
    q = qdd_gate(q, GATEID_X, 2);
    x3[2] = 0; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ONE);
    x3[2] = 1; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    test_assert(qdd_is_ordered(q, 3));

    q = qdd_circuit_swap(q, 0, 2);
    x3[2] = 0; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ONE);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    test_assert(qdd_is_ordered(q, 3));

    q = qdd_circuit_swap(q, 0, 1);
    x3[2] = 0; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ONE);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    test_assert(qdd_is_ordered(q, 3));

    // TODO: more tests

    if(VERBOSE) printf("qdd swap gates:           ok\n");
    return 0;
}

int test_cswap_circuit()
{
    QDD q;
    bool x3[] = {0,0,0};
    AMP a;

    uint32_t cs[3];
    cs[0] = 0; cs[1] = QDD_INVALID_VAR; cs[2] = QDD_INVALID_VAR; // control only on q0

    LACE_ME;

    x3[2]=1; x3[1]=0; x3[0]=0; 
    q = qdd_create_basis_state(3, x3);
    q = qdd_ccircuit(q, CIRCID_swap, cs, 1, 2); // control is |0>, nothing should happen
    x3[2] = 0; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ONE);
    x3[2] = 1; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);

    x3[2]=1; x3[1]=0; x3[0]=1; 
    q = qdd_create_basis_state(3, x3);
    q = qdd_ccircuit(q, CIRCID_swap, cs, 1, 2); // control is |1>, should swap q1 and q2
    x3[2] = 0; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ONE);
    x3[2] = 1; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);

    x3[2]=1; x3[1]=0; x3[0]=0; 
    q = qdd_create_basis_state(3, x3);
    q = qdd_gate(q, GATEID_H, 0); // input state: |10+> = 1/sqrt(2)(|100> + |101>)
    x3[2] = 0; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == comp_lookup(comp_make(1.0/flt_sqrt(2.0),0)));
    x3[2] = 1; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == comp_lookup(comp_make(1.0/flt_sqrt(2.0),0)));
    x3[2] = 1; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    q = qdd_ccircuit(q, CIRCID_swap, cs, 1, 2); // control is |+>, expected output: 1/sqrt(2)(|100> + |011>)
    x3[2] = 0; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == comp_lookup(comp_make(1.0/flt_sqrt(2.0),0)));
    x3[2] = 1; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == comp_lookup(comp_make(1.0/flt_sqrt(2.0),0)));
    x3[2] = 1; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);

    x3[2]=1; x3[1]=0; x3[0]=0; 
    q = qdd_create_basis_state(3, x3);
    q = qdd_gate(q, GATEID_H, 0); 
    q = qdd_gate(q, GATEID_Z, 0); // input state: |10-> = 1/sqrt(2)(|100> - |101>)
    x3[2] = 0; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == comp_lookup(comp_make(1.0/flt_sqrt(2.0),0)));
    x3[2] = 1; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == comp_lookup(comp_make(-1.0/flt_sqrt(2.0),0)));
    x3[2] = 1; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    q = qdd_ccircuit(q, CIRCID_swap, cs, 1, 2); // control is |->, expected output: 1/sqrt(2)(|100> - |011>)
    x3[2] = 0; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == comp_lookup(comp_make(-1.0/flt_sqrt(2.0),0)));
    x3[2] = 1; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == comp_lookup(comp_make(1.0/flt_sqrt(2.0),0)));
    x3[2] = 1; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);

    x3[2]=1; x3[1]=0; x3[0]=0; 
    q = qdd_create_basis_state(3, x3);
    q = qdd_gate(q, GATEID_H, 0); 
    q = qdd_gate(q, GATEID_S, 0); // input state: |10>|+i> = 1/sqrt(2)(|100> + i|101>)
    x3[2] = 0; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == comp_lookup(comp_make(1.0/flt_sqrt(2.0),0)));
    x3[2] = 1; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == comp_lookup(comp_make(0,1.0/flt_sqrt(2.0))));
    x3[2] = 1; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    q = qdd_ccircuit(q, CIRCID_swap, cs, 1, 2); // control is |+i>, expected output: 1/sqrt(2)(|100> + i|011>)
    x3[2] = 0; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 0; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == comp_lookup(comp_make(0,1.0/flt_sqrt(2.0))));
    x3[2] = 1; x3[1] = 0; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == comp_lookup(comp_make(1.0/flt_sqrt(2.0),0)));
    x3[2] = 1; x3[1] = 0; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 0; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);
    x3[2] = 1; x3[1] = 1; x3[0] = 1; a = qdd_get_amplitude(q, x3); test_assert(a == C_ZERO);

    // TODO: more tests ?

    if(VERBOSE) printf("qdd c-swap gates:         ok\n");
    return 0;
}

int test_tensor_product()
{
    LACE_ME;

    QDD q0, q1, qTest, qRef;
    bool x3[] = {0,1,0};
    bool x2[] = {1,0};
    bool x3_2[] = {0,1,0,1,0};
    bool x2_3[] = {1,0,0,1,0};

    q0 = qdd_create_basis_state(3, x3);
    q1 = qdd_create_basis_state(2, x2);
    test_assert(qdd_countnodes(q0) == 4);
    test_assert(qdd_countnodes(q1) == 3);

    // q0 (tensor) q1
    qTest = qdd_vec_tensor_prod(q0, q1, 3);
    qRef  = qdd_create_basis_state(5, x3_2);
    test_assert(qdd_is_ordered(qTest, 5));
    test_assert(qdd_countnodes(qTest) == 6);
    test_assert(qdd_equivalent(qTest, qRef, 5, false, true));
    test_assert(qdd_equivalent(qTest, qRef, 5, true, true));
    test_assert(qTest == qRef);

    // q1 (tensor) q0
    qTest = qdd_vec_tensor_prod(q1, q0, 2);
    qRef  = qdd_create_basis_state(5, x2_3);
    test_assert(qdd_is_ordered(qTest, 5));
    test_assert(qdd_countnodes(qTest) == 6);
    test_assert(qdd_equivalent(qTest, qRef, 5, false, true));
    test_assert(qdd_equivalent(qTest, qRef, 5, true, true));
    test_assert(qTest == qRef);


    // TODO: test on other states than basis states


    if(VERBOSE) printf("qdd tensor (on vecs):     ok\n");
    return 0;
}

// measurement of q0 + sanity checks
int test_measure_random_state(QDD qdd, BDDVAR nvars)
{
    QDD qm;
    int m; double p;

    LACE_ME;

    test_assert(qdd_is_ordered(qdd, nvars));
    test_assert(qdd_is_unitvector(qdd, nvars));
    qm = qdd_measure_q0(qdd, nvars, &m, &p);
    test_assert(qdd_is_ordered(qdd, nvars));
    test_assert(qdd_is_unitvector(qm, nvars));
    if ((flt_abs(p - 1.0) < cmap_get_tolerance()) || (flt_abs(p - 0.0) < cmap_get_tolerance())){
        qdd = qdd_remove_global_phase(qdd); // measurement removes global phase
        test_assert(qdd_equivalent(qdd, qm, nvars, false, false));
        test_assert(qdd_equivalent(qdd, qm, nvars, true,  false));
        test_assert(qm == qdd);
    } else {
        test_assert(qdd_countnodes(qm) < qdd_countnodes(qdd));
    }

    return 0;
}

int test_measurements()
{
    QDD q, qPM;
    AMP a;
    bool x2[] = {0,0};
    bool x3[] = {0,0,0};
    int m;
    int repeat = 10;
    double prob;

    LACE_ME;
    srand(time(NULL));

    // kets labeld as |q2, q1, q0>
    for(int i=0; i < repeat; i++) {
        // |000> 
        x3[2]=0; x3[1]=0; x3[0]=0;
        q   = qdd_create_basis_state(3, x3);
        qPM = qdd_measure_qubit(q, 0, 3, &m, &prob);
        test_assert(m == 0);
        test_assert(prob == 1.0);
        test_assert(qdd_equivalent(q, qPM, 3, false, false));
        test_assert(qdd_equivalent(q, qPM, 3, true, false));
        test_assert(q == qPM);

        // |010>
        x3[2]=0; x3[1]=1; x3[0]=0;
        q   = qdd_create_basis_state(3, x3);
        qPM = qdd_measure_qubit(q, 0, 3, &m, &prob);
        test_assert(m == 0);
        test_assert(prob == 1.0);
        test_assert(qdd_equivalent(q, qPM, 3, false, false));
        test_assert(qdd_equivalent(q, qPM, 3, true, false));
        test_assert(q == qPM);

        // |011>
        x3[2]=0; x3[1]=1; x3[0]=1;
        q   = qdd_create_basis_state(3, x3);
        qPM = qdd_measure_qubit(q, 0, 3, &m, &prob);
        test_assert(m == 1);
        test_assert(prob == 0.0);
        test_assert(qdd_equivalent(q, qPM, 3, false, false));
        test_assert(qdd_equivalent(q, qPM, 3, true, false));
        test_assert(q == qPM);

        // |00+>
        x3[2]=0; x3[1]=0; x3[0]=0;
        q   = qdd_create_basis_state(3, x3);
        q   = qdd_gate(q, GATEID_H, 0);
        qPM = qdd_measure_qubit(q, 0, 3, &m, &prob);
        test_assert(flt_abs(prob - 0.5) < cmap_get_tolerance());
        x3[2]=0; x3[1]=0; x3[0]=m; // either |000> or |001> depending on m
        q = qdd_create_basis_state(3, x3); 
        test_assert(qdd_equivalent(q, qPM, 3, false, true));
        test_assert(qdd_equivalent(q, qPM, 3, true, false));
        test_assert(q == qPM);

        // |0+0>
        x3[2]=0; x3[1]=0; x3[0]=0;
        q   = qdd_create_basis_state(3, x3);
        q   = qdd_gate(q, GATEID_H, 1);
        qPM = qdd_measure_qubit(q, 1, 3, &m, &prob);
        test_assert(flt_abs(prob - 0.5) < cmap_get_tolerance());
        x3[2]=0; x3[1]=m; x3[0]=0; // either |000> or |010> depending on m
        q = qdd_create_basis_state(3, x3);
        test_assert(qdd_equivalent(q, qPM, 3, false, true));
        test_assert(qdd_equivalent(q, qPM, 3, true, false));
        test_assert(q == qPM);

        // |+00>
        x3[2]=0; x3[1]=0; x3[0]=0;
        q   = qdd_create_basis_state(3, x3);
        q   = qdd_gate(q, GATEID_H, 2);
        qPM = qdd_measure_qubit(q, 2, 3, &m, &prob);
        test_assert(flt_abs(prob - 0.5) < cmap_get_tolerance());
        x3[2]=m; x3[1]=0; x3[0]=0; // either |000> or |100> depending on m
        q = qdd_create_basis_state(3, x3);
        test_assert(qdd_equivalent(q, qPM, 3, false, true));
        test_assert(qdd_equivalent(q, qPM, 3, true, false));
        test_assert(q == qPM);

        // |00->
        x3[2]=0; x3[1]=0; x3[0]=1;
        q   = qdd_create_basis_state(3, x3);
        q   = qdd_gate(q, GATEID_H, 0);
        qPM = qdd_measure_qubit(q, 0, 3, &m, &prob);
        test_assert(flt_abs(prob - 0.5) < cmap_get_tolerance());
        x3[2]=0; x3[1]=0; x3[0]=m; // either |000> or |001> depending on m
        q = qdd_create_basis_state(3, x3); 
        test_assert(qdd_equivalent(q, qPM, 3, false, true));
        test_assert(qdd_equivalent(q, qPM, 3, true, false));
        test_assert(q == qPM);
        
        // |+++>, measure q0
        x3[2]=0; x3[1]=0; x3[0]=0;
        q   = qdd_create_basis_state(3, x3);
        q   = qdd_gate(q, GATEID_H, 0);
        q   = qdd_gate(q, GATEID_H, 1);
        q   = qdd_gate(q, GATEID_H, 2);
        qPM = qdd_measure_qubit(q, 0, 3, &m, &prob); 
        test_assert(flt_abs(prob - 0.5) < cmap_get_tolerance());
        x3[2]=0; x3[1]=0; x3[0]=m; // either |++0> or |++1> depending on m
        q = qdd_create_basis_state(3, x3); 
        q = qdd_gate(q, GATEID_H, 1);
        q = qdd_gate(q, GATEID_H, 2);
        test_assert(qdd_equivalent(q, qPM, 3, false, true));
        test_assert(qdd_equivalent(q, qPM, 3, true, false));
        test_assert(q == qPM); 

        // |+++>, measure q1
        x3[2]=0; x3[1]=0; x3[0]=0;
        q   = qdd_create_basis_state(3, x3);
        q   = qdd_gate(q, GATEID_H, 0);
        q   = qdd_gate(q, GATEID_H, 1);
        q   = qdd_gate(q, GATEID_H, 2);
        qPM = qdd_measure_qubit(q, 1, 3, &m, &prob);
        test_assert(flt_abs(prob - 0.5) < cmap_get_tolerance());
        x3[2]=0; x3[1]=m; x3[0]=0; // either |+0+> or |+1+> depending on m
        q = qdd_create_basis_state(3, x3); 
        q = qdd_gate(q, GATEID_H, 0);
        q = qdd_gate(q, GATEID_H, 2);
        test_assert(qdd_equivalent(q, qPM, 3, false, true));
        test_assert(qdd_equivalent(q, qPM, 3, true, false));
        test_assert(q == qPM);

        // [1/2, 1/2, 1/2, -1/2] = 1/2(|00> + |01> + |10> - |11>)_(q1, q0)
        x2[1]=0; x2[0]=0;
        q   = qdd_create_basis_state(2, x2);
        q   = qdd_gate(q, GATEID_H, 0);
        q   = qdd_gate(q, GATEID_H, 1);
        q   = qdd_cgate(q,GATEID_Z, 0, 1);
        x2[1]=0; x2[0]=0; a = qdd_get_amplitude(q, x2); test_assert(a == comp_lookup(comp_make(0.5,0)));
        x2[1]=0; x2[0]=1; a = qdd_get_amplitude(q, x2); test_assert(a == comp_lookup(comp_make(0.5,0)));
        x2[1]=1; x2[0]=0; a = qdd_get_amplitude(q, x2); test_assert(a == comp_lookup(comp_make(0.5,0)));
        x2[1]=1; x2[0]=1; a = qdd_get_amplitude(q, x2); test_assert(a == comp_lookup(comp_make(-0.5,0)));
        qPM = qdd_measure_qubit(q, 0, 2, &m, &prob);
        test_assert(flt_abs(prob - 0.5) < cmap_get_tolerance());
        if (m == 0) { // expect 1/sqrt(2)(|00> + |10>)
            x2[1]=0; x2[0]=0; a = qdd_get_amplitude(qPM, x2); test_assert(a == comp_lookup(comp_make(1.0/flt_sqrt(2.0),0)));
            x2[1]=0; x2[0]=1; a = qdd_get_amplitude(qPM, x2); test_assert(a == C_ZERO);
            x2[1]=1; x2[0]=0; a = qdd_get_amplitude(qPM, x2); test_assert(a == comp_lookup(comp_make(1.0/flt_sqrt(2.0),0)));
            x2[1]=1; x2[0]=1; a = qdd_get_amplitude(qPM, x2); test_assert(a == C_ZERO);
        }
        if (m == 1) { // expect 1/sqrt(2)(|01> - |11>)
            x2[1]=0; x2[0]=0; a = qdd_get_amplitude(qPM, x2); test_assert(a == C_ZERO);
            x2[1]=0; x2[0]=1; a = qdd_get_amplitude(qPM, x2); test_assert(a == comp_lookup(comp_make(1.0/flt_sqrt(2.0),0)));
            x2[1]=1; x2[0]=0; a = qdd_get_amplitude(qPM, x2); test_assert(a == C_ZERO);
            x2[1]=1; x2[0]=1; a = qdd_get_amplitude(qPM, x2); test_assert(a == comp_lookup(comp_make(-1.0/flt_sqrt(2.0),0)));
        }
    }

    // Test measure all
    bool ms[3] = {0};
    int m_zer[3] = {0};
    for(int i=0; i < repeat; i++) {
        
        // |000>  
        x3[2]=0; x3[1]=0; x3[0]=0;
        q   = qdd_create_basis_state(3, x3);
        qPM = qdd_measure_all(q, 3, ms, &prob);
        test_assert(ms[0] == 0);
        test_assert(ms[1] == 0);
        test_assert(ms[2] == 0);
        test_assert(prob == 1.0);
        test_assert(qdd_equivalent(q, qPM, 3, false, false));
        test_assert(qdd_equivalent(q, qPM, 3, true, false));
        test_assert(q == qPM);

        // |010>
         x3[2]=0; x3[1]=1; x3[0]=0;
        q   = qdd_create_basis_state(3, x3);
        qPM = qdd_measure_all(q, 3, ms, &prob);
        test_assert(ms[0] == 0);
        test_assert(ms[1] == 1);
        test_assert(ms[2] == 0);
        test_assert(prob == 1.0);
        test_assert(qdd_equivalent(q, qPM, 3, false, false));
        test_assert(qdd_equivalent(q, qPM, 3, true, false));
        test_assert(q == qPM);

        // |011>
        x3[2]=0; x3[1]=1; x3[0]=1;
        q   = qdd_create_basis_state(3, x3);
        qPM = qdd_measure_all(q, 3, ms, &prob);
        test_assert(ms[0] == 1);
        test_assert(ms[1] == 1);
        test_assert(ms[2] == 0);
        test_assert(prob == 1.0);
        test_assert(qdd_equivalent(q, qPM, 3, false, false));
        test_assert(qdd_equivalent(q, qPM, 3, true, false));
        test_assert(q == qPM);

       // |00+>
        x3[2]=0; x3[1]=0; x3[0]=0;
        q   = qdd_create_basis_state(3, x3);
        q   = qdd_gate(q, GATEID_H, 0);
        qPM = qdd_measure_all(q, 3, ms, &prob);
        test_assert(flt_abs(prob - 0.5) < cmap_get_tolerance());
        x3[2]=0; x3[1]=0; x3[0]=ms[0]; // either |000> or |001> depending on m
        q = qdd_create_basis_state(3, x3); 
        test_assert(qdd_equivalent(q, qPM, 3, false, true));
        test_assert(qdd_equivalent(q, qPM, 3, true, false));
        test_assert(q == qPM);
        if (ms[0] == 0) m_zer[0] += 1;

        // |0+0>
        x3[2]=0; x3[1]=0; x3[0]=0;
        q   = qdd_create_basis_state(3, x3);
        q   = qdd_gate(q, GATEID_H, 1);
        qPM = qdd_measure_all(q, 3, ms, &prob);
        test_assert(flt_abs(prob - 0.5) < cmap_get_tolerance());
        x3[2]=0; x3[1]=ms[1]; x3[0]=0; // either |000> or |010> depending on m
        q = qdd_create_basis_state(3, x3);
        test_assert(qdd_equivalent(q, qPM, 3, false, true));
        test_assert(qdd_equivalent(q, qPM, 3, true, false));
        test_assert(q == qPM);
        if (ms[1] == 0) m_zer[1] += 1;

        // |+00>
        x3[2]=0; x3[1]=0; x3[0]=0;
        q   = qdd_create_basis_state(3, x3);
        q   = qdd_gate(q, GATEID_H, 2);
        qPM = qdd_measure_all(q, 3, ms, &prob);
        test_assert(flt_abs(prob - 0.5) < cmap_get_tolerance());
        x3[2]=ms[2]; x3[1]=0; x3[0]=0; // either |000> or |100> depending on m
        q = qdd_create_basis_state(3, x3);
        test_assert(qdd_equivalent(q, qPM, 3, false, true));
        test_assert(qdd_equivalent(q, qPM, 3, true, false));
        test_assert(q == qPM);
        if (ms[2] == 0) m_zer[2] += 1;
    }
    
    // TODO: more tests

    if(VERBOSE) printf("qdd measurements:         ok\n");
    return 0;
}

int test_QFT()
{
    QDD q3, q5, qref3, qref5;
    AMP a;

    LACE_ME;

    // 3 qubit QFT
    bool x3[] = {0,1,1}; // little endian (q0, q1, q2)
    q3 = qdd_create_basis_state(3, x3);
    qref3 = qdd_create_basis_state(3, x3);
    q3 = qdd_circuit(q3, CIRCID_QFT, 0, 2);
    q3 = qdd_circuit(q3, CIRCID_reverse_range, 0, 2);

    // check approx equal against output from qiskit
    x3[2]=0; x3[1]=0; x3[0]=0; a = qdd_get_amplitude(q3, x3); test_assert(comp_approx_equal(comp_value(a),comp_make(3.5355339059327384e-01,-8.6595605623549353e-17)));
    x3[2]=0; x3[1]=0; x3[0]=1; a = qdd_get_amplitude(q3, x3); test_assert(comp_approx_equal(comp_value(a),comp_make(-3.5355339059327384e-01,8.6595605623549353e-17)));
    x3[2]=0; x3[1]=1; x3[0]=0; a = qdd_get_amplitude(q3, x3); test_assert(comp_approx_equal(comp_value(a),comp_make(-1.0824450702943669e-16,-3.5355339059327384e-01)));
    x3[2]=0; x3[1]=1; x3[0]=1; a = qdd_get_amplitude(q3, x3); test_assert(comp_approx_equal(comp_value(a),comp_make(1.0824450702943669e-16,3.5355339059327384e-01)));
    x3[2]=1; x3[1]=0; x3[0]=0; a = qdd_get_amplitude(q3, x3); test_assert(comp_approx_equal(comp_value(a),comp_make(-2.5000000000000000e-01,2.5000000000000017e-01)));
    x3[2]=1; x3[1]=0; x3[0]=1; a = qdd_get_amplitude(q3, x3); test_assert(comp_approx_equal(comp_value(a),comp_make(2.5000000000000000e-01,-2.5000000000000017e-01)));
    x3[2]=1; x3[1]=1; x3[0]=0; a = qdd_get_amplitude(q3, x3); test_assert(comp_approx_equal(comp_value(a),comp_make(2.5000000000000017e-01,2.5000000000000000e-01)));
    x3[2]=1; x3[1]=1; x3[0]=1; a = qdd_get_amplitude(q3, x3); test_assert(comp_approx_equal(comp_value(a),comp_make(-2.5000000000000017e-01,-2.5000000000000000e-01)));
    test_assert(qdd_is_unitvector(q3, 3));
    test_assert(qdd_is_ordered(q3, 3));

    // inverse QFT
    q3 = qdd_circuit(q3, CIRCID_reverse_range, 0, 2);
    q3 = qdd_circuit(q3, CIRCID_QFT_inv, 0, 2);
    test_assert(qdd_equivalent(q3, qref3, 3, false, false));
    test_assert(qdd_equivalent(q3, qref3, 3, true, false));
    test_assert(q3 == qref3);

    // 5 qubit QFT
    bool x5[] = {0,1,1,0,1};
    q5 = qdd_create_basis_state(5, x5);
    qref5 = qdd_create_basis_state(5, x5);
    q5 = qdd_circuit(q5, CIRCID_QFT, 0, 4);
    q5 = qdd_circuit(q5, CIRCID_reverse_range, 0, 4);
    test_assert(qdd_is_ordered(q5, 5));

    // check approx equal against output from qiskit
    x5[4]=0; x5[3]=0; x5[2]=0; x5[1]=0; x5[0]=0; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(1.7677669529663692e-01,-6.4946704217662027e-17)));
    x5[4]=0; x5[3]=0; x5[2]=0; x5[1]=0; x5[0]=1; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(-1.7677669529663692e-01,6.4946704217662027e-17)));
    x5[4]=0; x5[3]=0; x5[2]=0; x5[1]=1; x5[0]=0; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(7.5771154920605696e-17,1.7677669529663692e-01)));
    x5[4]=0; x5[3]=0; x5[2]=0; x5[1]=1; x5[0]=1; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(-7.5771154920605696e-17,-1.7677669529663692e-01)));
    x5[4]=0; x5[3]=0; x5[2]=1; x5[1]=0; x5[0]=0; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(-1.2500000000000011e-01,-1.2499999999999999e-01)));
    x5[4]=0; x5[3]=0; x5[2]=1; x5[1]=0; x5[0]=1; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(1.2500000000000011e-01,1.2499999999999999e-01)));
    x5[4]=0; x5[3]=0; x5[2]=1; x5[1]=1; x5[0]=0; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(1.2499999999999999e-01,-1.2500000000000011e-01)));
    x5[4]=0; x5[3]=0; x5[2]=1; x5[1]=1; x5[0]=1; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(-1.2499999999999999e-01,1.2500000000000011e-01)));
    x5[4]=0; x5[3]=1; x5[2]=0; x5[1]=0; x5[0]=0; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(6.7649512518274585e-02,-1.6332037060954713e-01)));
    x5[4]=0; x5[3]=1; x5[2]=0; x5[1]=0; x5[0]=1; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(-6.7649512518274585e-02,1.6332037060954713e-01)));
    x5[4]=0; x5[3]=1; x5[2]=0; x5[1]=1; x5[0]=0; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(1.6332037060954713e-01,6.7649512518274571e-02)));
    x5[4]=0; x5[3]=1; x5[2]=0; x5[1]=1; x5[0]=1; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(-1.6332037060954713e-01,-6.7649512518274571e-02)));
    x5[4]=0; x5[3]=1; x5[2]=1; x5[1]=0; x5[0]=0; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(-1.6332037060954710e-01,6.7649512518274696e-02)));
    x5[4]=0; x5[3]=1; x5[2]=1; x5[1]=0; x5[0]=1; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(1.6332037060954710e-01,-6.7649512518274696e-02)));
    x5[4]=0; x5[3]=1; x5[2]=1; x5[1]=1; x5[0]=0; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(-6.7649512518274710e-02,-1.6332037060954710e-01)));
    x5[4]=0; x5[3]=1; x5[2]=1; x5[1]=1; x5[0]=1; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(6.7649512518274710e-02,1.6332037060954710e-01)));
    x5[4]=1; x5[3]=0; x5[2]=0; x5[1]=0; x5[0]=0; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(-1.4698445030241986e-01,9.8211869798387877e-02)));
    x5[4]=1; x5[3]=0; x5[2]=0; x5[1]=0; x5[0]=1; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(1.4698445030241986e-01,-9.8211869798387877e-02)));
    x5[4]=1; x5[3]=0; x5[2]=0; x5[1]=1; x5[0]=0; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(-9.8211869798387877e-02,-1.4698445030241986e-01)));
    x5[4]=1; x5[3]=0; x5[2]=0; x5[1]=1; x5[0]=1; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(9.8211869798387877e-02,1.4698445030241986e-01)));
    x5[4]=1; x5[3]=0; x5[2]=1; x5[1]=0; x5[0]=0; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(1.7337998066526852e-01,3.4487422410367806e-02)));
    x5[4]=1; x5[3]=0; x5[2]=1; x5[1]=0; x5[0]=1; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(-1.7337998066526852e-01,-3.4487422410367806e-02)));
    x5[4]=1; x5[3]=0; x5[2]=1; x5[1]=1; x5[0]=0; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(-3.4487422410367799e-02,1.7337998066526852e-01)));
    x5[4]=1; x5[3]=0; x5[2]=1; x5[1]=1; x5[0]=1; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(3.4487422410367799e-02,-1.7337998066526852e-01)));
    x5[4]=1; x5[3]=1; x5[2]=0; x5[1]=0; x5[0]=0; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(3.4487422410367972e-02,1.7337998066526850e-01)));
    x5[4]=1; x5[3]=1; x5[2]=0; x5[1]=0; x5[0]=1; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(-3.4487422410367972e-02,-1.7337998066526850e-01)));
    x5[4]=1; x5[3]=1; x5[2]=0; x5[1]=1; x5[0]=0; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(-1.7337998066526850e-01,3.4487422410367986e-02)));
    x5[4]=1; x5[3]=1; x5[2]=0; x5[1]=1; x5[0]=1; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(1.7337998066526850e-01,-3.4487422410367986e-02)));
    x5[4]=1; x5[3]=1; x5[2]=1; x5[1]=0; x5[0]=0; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(9.8211869798387752e-02,-1.4698445030241994e-01)));
    x5[4]=1; x5[3]=1; x5[2]=1; x5[1]=0; x5[0]=1; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(-9.8211869798387752e-02,1.4698445030241994e-01)));
    x5[4]=1; x5[3]=1; x5[2]=1; x5[1]=1; x5[0]=0; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(1.4698445030241994e-01,9.8211869798387752e-02)));
    x5[4]=1; x5[3]=1; x5[2]=1; x5[1]=1; x5[0]=1; a = qdd_get_amplitude(q5, x5); test_assert(comp_approx_equal(comp_value(a), comp_make(-1.4698445030241994e-01,-9.8211869798387752e-02)));
    test_assert(qdd_is_unitvector(q5, 5));

    // inverse QFT
    q5 = qdd_circuit(q5, CIRCID_reverse_range, 0, 4);
    q5 = qdd_circuit(q5, CIRCID_QFT_inv, 0, 4);
    test_assert(qdd_equivalent(q5, qref5, 5, false, false));
    test_assert(qdd_equivalent(q5, qref5, 5, true, false));
    test_assert(q5 == qref5);
    
    
    if(VERBOSE) printf("qdd QFT:                  ok\n");
    return 0;
}

int test_5qubit_circuit()
{
    QDD q, qref;
    uint64_t node_count;
    int n_qubits = 5;
    bool x5[] = {0,0,0,0,0};
    FILE *fp;
    
    LACE_ME;

    // 5 qubit state
    qref = qdd_create_basis_state(n_qubits, x5);
    q    = qdd_create_basis_state(n_qubits, x5);

    // 32 gates 
    q = qdd_cgate(q, GATEID_Z, 1, 2);       q = qdd_gate(q, GATEID_X, 2);       test_assert(qdd_is_unitvector(q, 5));
    q = qdd_cgate(q, GATEID_Z, 3, 4);       q = qdd_cgate(q, GATEID_X, 1, 3);   test_assert(qdd_is_unitvector(q, 5));
    q = qdd_gate(q, GATEID_Z, 1);           q = qdd_gate(q, GATEID_H, 4);       test_assert(qdd_is_unitvector(q, 5));
    q = qdd_gate(q, GATEID_H, 0);           q = qdd_cgate(q, GATEID_X, 1, 3);   test_assert(qdd_is_unitvector(q, 5));
    q = qdd_gate(q, GATEID_H, 3);           q = qdd_gate(q, GATEID_H, 0);       test_assert(qdd_is_unitvector(q, 5));
    q = qdd_gate(q, GATEID_Z, 1);           q = qdd_cgate(q, GATEID_X, 1, 2);   test_assert(qdd_is_unitvector(q, 5));
    q = qdd_gate(q, GATEID_X, 1);           q = qdd_cgate(q, GATEID_X, 0, 4);   test_assert(qdd_is_unitvector(q, 5));
    q = qdd_gate(q, GATEID_H, 4);           q = qdd_cgate(q, GATEID_X, 0, 1);   test_assert(qdd_is_unitvector(q, 5));
    q = qdd_cgate(q, GATEID_Z, 0, 4);       q = qdd_cgate(q, GATEID_Z, 0, 4);   test_assert(qdd_is_unitvector(q, 5));
    q = qdd_cgate(q, GATEID_Z, 2, 3);       q = qdd_gate(q, GATEID_Z, 0);       test_assert(qdd_is_unitvector(q, 5));
    q = qdd_cgate(q, GATEID_X, 3, 4);       q = qdd_cgate(q, GATEID_Z, 0, 2);   test_assert(qdd_is_unitvector(q, 5));
    q = qdd_gate(q, GATEID_Z, 3);           q = qdd_cgate(q, GATEID_Z, 1, 3);   test_assert(qdd_is_unitvector(q, 5));
    q = qdd_cgate(q, GATEID_X, 0, 3);       q = qdd_cgate(q, GATEID_Z, 0, 4);   test_assert(qdd_is_unitvector(q, 5));
    q = qdd_cgate(q, GATEID_Z, 0, 4);       q = qdd_gate(q, GATEID_H, 1);       test_assert(qdd_is_unitvector(q, 5));
    q = qdd_gate(q, GATEID_H, 3);           q = qdd_gate(q, GATEID_X, 1);       test_assert(qdd_is_unitvector(q, 5));
    q = qdd_cgate(q, GATEID_Z, 1, 2);       q = qdd_gate(q, GATEID_Z, 1);       test_assert(qdd_is_unitvector(q, 5));
    node_count = qdd_countnodes(q);
    if (test_measure_random_state(q, n_qubits)) return 1;

    // inverse
    q = qdd_gate(q, GATEID_Z, 1);           q = qdd_cgate(q, GATEID_Z, 1, 2);   test_assert(qdd_is_unitvector(q, 5));
    q = qdd_gate(q, GATEID_X, 1);           q = qdd_gate(q, GATEID_H, 3);       test_assert(qdd_is_unitvector(q, 5));
    q = qdd_gate(q, GATEID_H, 1);           q = qdd_cgate(q, GATEID_Z, 0, 4);   test_assert(qdd_is_unitvector(q, 5));
    q = qdd_cgate(q, GATEID_Z, 0, 4);       q = qdd_cgate(q, GATEID_X, 0, 3);   test_assert(qdd_is_unitvector(q, 5));
    q = qdd_cgate(q, GATEID_Z, 1, 3);       q = qdd_gate(q, GATEID_Z, 3);       test_assert(qdd_is_unitvector(q, 5));
    q = qdd_cgate(q, GATEID_Z, 0, 2);       q = qdd_cgate(q, GATEID_X, 3, 4);   test_assert(qdd_is_unitvector(q, 5));
    q = qdd_gate(q, GATEID_Z, 0);           q = qdd_cgate(q, GATEID_Z, 2, 3);   test_assert(qdd_is_unitvector(q, 5));
    q = qdd_cgate(q, GATEID_Z, 0, 4);       q = qdd_cgate(q, GATEID_Z, 0, 4);   test_assert(qdd_is_unitvector(q, 5));
    q = qdd_cgate(q, GATEID_X, 0, 1);       q = qdd_gate(q, GATEID_H, 4);       test_assert(qdd_is_unitvector(q, 5));
    q = qdd_cgate(q, GATEID_X, 0, 4);       q = qdd_gate(q, GATEID_X, 1);       test_assert(qdd_is_unitvector(q, 5));
    q = qdd_cgate(q, GATEID_X, 1, 2);       q = qdd_gate(q, GATEID_Z, 1);       test_assert(qdd_is_unitvector(q, 5));
    q = qdd_gate(q, GATEID_H, 0);           q = qdd_gate(q, GATEID_H, 3);       test_assert(qdd_is_unitvector(q, 5));
    q = qdd_cgate(q, GATEID_X, 1, 3);       q = qdd_gate(q, GATEID_H, 0);       test_assert(qdd_is_unitvector(q, 5));
    q = qdd_gate(q, GATEID_H, 4);           q = qdd_gate(q, GATEID_Z, 1);       test_assert(qdd_is_unitvector(q, 5));
    q = qdd_cgate(q, GATEID_X, 1, 3);       q = qdd_cgate(q, GATEID_Z, 3, 4);   test_assert(qdd_is_unitvector(q, 5));
    q = qdd_gate(q, GATEID_X, 2);           q = qdd_cgate(q, GATEID_Z, 1, 2);   test_assert(qdd_is_unitvector(q, 5));

    test_assert(qdd_equivalent(q, qref, n_qubits, false, VERBOSE)); // check approx equiv
    test_assert(qdd_equivalent(q, qref, n_qubits, true,  VERBOSE)); // check exact equiv
    test_assert(q == qref);

    fp = fopen("5_qubit_res.dot", "w");
    qdd_fprintdot(fp, q, true);
    fclose(fp);

    if(VERBOSE) printf("qdd 5 qubit circuit:      ok (%ld nodes)\n", node_count);
    return 0;
}

int test_10qubit_circuit()
{
    QDD q, qref;
    uint64_t node_count;
    bool x10[] = {0,0,0,0,0,0,0,0,0,0};
    
    LACE_ME;

    // 10 qubit state
    qref = qdd_create_basis_state(10, x10);
    q    = qdd_create_basis_state(10, x10);

    // 30 random* Clifford gates            *chosen by a fair dice roll
    q = qdd_cgate(q, GATEID_X, 1, 3);       q = qdd_gate(q, GATEID_H, 0);       test_assert(qdd_is_unitvector(q, 10));
    q = qdd_gate(q, GATEID_X, 6);           q = qdd_cgate(q, GATEID_X, 6, 9);   test_assert(qdd_is_unitvector(q, 10));
    q = qdd_gate(q, GATEID_H, 4);           q = qdd_cgate(q, GATEID_X, 3, 5);   test_assert(qdd_is_unitvector(q, 10));
    q = qdd_gate(q, GATEID_H, 1);           q = qdd_gate(q, GATEID_X, 1);       test_assert(qdd_is_unitvector(q, 10));
    q = qdd_cgate(q, GATEID_X, 3, 8);       q = qdd_cgate(q, GATEID_Z, 3, 6);   test_assert(qdd_is_unitvector(q, 10));
    q = qdd_gate(q, GATEID_Z, 3);           q = qdd_cgate(q, GATEID_X, 0, 7);   test_assert(qdd_is_unitvector(q, 10));
    q = qdd_cgate(q, GATEID_X, 1, 9);       q = qdd_gate(q, GATEID_H, 4);       test_assert(qdd_is_unitvector(q, 10));
    q = qdd_cgate(q, GATEID_X, 0, 2);       q = qdd_gate(q, GATEID_X, 2);       test_assert(qdd_is_unitvector(q, 10));
    q = qdd_cgate(q, GATEID_X, 5, 8);       q = qdd_cgate(q, GATEID_X, 0, 4);   test_assert(qdd_is_unitvector(q, 10));
    q = qdd_cgate(q, GATEID_X, 0, 8);       q = qdd_cgate(q, GATEID_X, 6, 9);   test_assert(qdd_is_unitvector(q, 10));
    q = qdd_cgate(q, GATEID_X, 0, 9);       q = qdd_gate(q, GATEID_X, 9);       test_assert(qdd_is_unitvector(q, 10));
    q = qdd_cgate(q, GATEID_X, 4, 9);       q = qdd_cgate(q, GATEID_Z, 2, 7);   test_assert(qdd_is_unitvector(q, 10));
    q = qdd_cgate(q, GATEID_Z, 7, 8);       q = qdd_gate(q, GATEID_X, 7);       test_assert(qdd_is_unitvector(q, 10));
    q = qdd_gate(q, GATEID_Z, 2);           q = qdd_gate(q, GATEID_Z, 7);       test_assert(qdd_is_unitvector(q, 10));
    q = qdd_gate(q, GATEID_X, 6);           q = qdd_gate(q, GATEID_X, 1);       test_assert(qdd_is_unitvector(q, 10));
    node_count = qdd_countnodes(q);
    if (test_measure_random_state(q, 10)) return 1;

    // inverse
    q = qdd_gate(q, GATEID_X, 1);           q = qdd_gate(q, GATEID_X, 6);       test_assert(qdd_is_unitvector(q, 10));
    q = qdd_gate(q, GATEID_Z, 7);           q = qdd_gate(q, GATEID_Z, 2);       test_assert(qdd_is_unitvector(q, 10));
    q = qdd_gate(q, GATEID_X, 7);           q = qdd_cgate(q, GATEID_Z, 7, 8);   test_assert(qdd_is_unitvector(q, 10));
    q = qdd_cgate(q, GATEID_Z, 2, 7);       q = qdd_cgate(q, GATEID_X, 4, 9);   test_assert(qdd_is_unitvector(q, 10));
    q = qdd_gate(q, GATEID_X, 9);           q = qdd_cgate(q, GATEID_X, 0, 9);   test_assert(qdd_is_unitvector(q, 10));
    q = qdd_cgate(q, GATEID_X, 6, 9);       q = qdd_cgate(q, GATEID_X, 0, 8);   test_assert(qdd_is_unitvector(q, 10));
    q = qdd_cgate(q, GATEID_X, 0, 4);       q = qdd_cgate(q, GATEID_X, 5, 8);   test_assert(qdd_is_unitvector(q, 10));
    q = qdd_gate(q, GATEID_X, 2);           q = qdd_cgate(q, GATEID_X, 0, 2);   test_assert(qdd_is_unitvector(q, 10));
    q = qdd_gate(q, GATEID_H, 4);           q = qdd_cgate(q, GATEID_X, 1, 9);   test_assert(qdd_is_unitvector(q, 10));
    q = qdd_cgate(q, GATEID_X, 0, 7);       q = qdd_gate(q, GATEID_Z, 3);       test_assert(qdd_is_unitvector(q, 10));
    q = qdd_cgate(q, GATEID_Z, 3, 6);       q = qdd_cgate(q, GATEID_X, 3, 8);   test_assert(qdd_is_unitvector(q, 10));
    q = qdd_gate(q, GATEID_X, 1);           q = qdd_gate(q, GATEID_H, 1);       test_assert(qdd_is_unitvector(q, 10));
    q = qdd_cgate(q, GATEID_X, 3, 5);       q = qdd_gate(q, GATEID_H, 4);       test_assert(qdd_is_unitvector(q, 10));
    q = qdd_cgate(q, GATEID_X, 6, 9);       q = qdd_gate(q, GATEID_X, 6);       test_assert(qdd_is_unitvector(q, 10));
    q = qdd_gate(q, GATEID_H, 0);           q = qdd_cgate(q, GATEID_X, 1, 3);   test_assert(qdd_is_unitvector(q, 10));

    test_assert(qdd_equivalent(q, qref, 10, false, VERBOSE)); // check approx equiv
    test_assert(qdd_equivalent(q, qref, 10, true,  VERBOSE)); // check exact equiv
    test_assert(q == qref);

    if(VERBOSE) printf("qdd 10 qubit circuit:     ok (%ld nodes)\n", node_count);
    return 0;
}

int test_20qubit_circuit()
{
    QDD q, qref;
    uint64_t node_count;
    bool x20[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

    LACE_ME;

    // 20 qubit state
    qref = qdd_create_basis_state(20, x20);
    q    = qdd_create_basis_state(20, x20);

    // 100 gates
    q = qdd_cgate(q, GATEID_Z, 4, 18);      q = qdd_gate(q, GATEID_H, 16);          q = qdd_cgate(q, GATEID_X, 1, 12);      q = qdd_gate(q, GATEID_Z, 4);
    q = qdd_cgate(q, GATEID_Z, 9, 10);      q = qdd_cgate(q, GATEID_Z, 0, 4);       q = qdd_cgate(q, GATEID_Z, 1, 16);      q = qdd_cgate(q, GATEID_Z, 13, 16);
    q = qdd_cgate(q, GATEID_Z, 7, 11);      q = qdd_cgate(q, GATEID_X, 3, 5);       q = qdd_cgate(q, GATEID_Z, 1, 4);       q = qdd_cgate(q, GATEID_X, 6, 16);
    q = qdd_cgate(q, GATEID_X, 3, 18);      q = qdd_cgate(q, GATEID_X, 2, 15);      q = qdd_cgate(q, GATEID_X, 7, 10);      q = qdd_gate(q, GATEID_Z, 6);
    q = qdd_cgate(q, GATEID_X, 3, 6);       q = qdd_cgate(q, GATEID_Z, 11, 16);     q = qdd_cgate(q, GATEID_X, 5, 19);      q = qdd_gate(q, GATEID_Z, 18);
    q = qdd_cgate(q, GATEID_Z, 14, 15);     q = qdd_cgate(q, GATEID_Z, 10, 12);     q = qdd_gate(q, GATEID_H, 8);           q = qdd_gate(q, GATEID_X, 9);
    q = qdd_gate(q, GATEID_X, 8);           q = qdd_cgate(q, GATEID_X, 7, 18);      q = qdd_gate(q, GATEID_X, 17);          q = qdd_gate(q, GATEID_Z, 11);
    q = qdd_cgate(q, GATEID_X, 12, 16);     q = qdd_gate(q, GATEID_X, 18);          q = qdd_gate(q, GATEID_Z, 4);           q = qdd_gate(q, GATEID_X, 18);
    q = qdd_cgate(q, GATEID_X, 4, 10);      q = qdd_gate(q, GATEID_X, 15);          q = qdd_cgate(q, GATEID_Z, 16, 18);     q = qdd_cgate(q, GATEID_Z, 0, 15);
    q = qdd_cgate(q, GATEID_X, 7, 10);      q = qdd_gate(q, GATEID_X, 18);          q = qdd_gate(q, GATEID_Z, 16);          q = qdd_cgate(q, GATEID_X, 7, 18);
    q = qdd_gate(q, GATEID_X, 16);          q = qdd_gate(q, GATEID_X, 2);           q = qdd_cgate(q, GATEID_X, 9, 10);      q = qdd_gate(q, GATEID_X, 6);
    q = qdd_gate(q, GATEID_X, 18);          q = qdd_gate(q, GATEID_Z, 11);          q = qdd_cgate(q, GATEID_Z, 4, 5);       q = qdd_gate(q, GATEID_X, 1);
    q = qdd_cgate(q, GATEID_Z, 2, 19);      q = qdd_cgate(q, GATEID_X, 8, 9);       q = qdd_cgate(q, GATEID_Z, 10, 12);     q = qdd_cgate(q, GATEID_Z, 11, 16);
    q = qdd_cgate(q, GATEID_X, 13, 19);     q = qdd_cgate(q, GATEID_Z, 1, 3);       q = qdd_gate(q, GATEID_X, 6);           q = qdd_gate(q, GATEID_X, 15);
    q = qdd_gate(q, GATEID_Z, 0);           q = qdd_cgate(q, GATEID_X, 0, 15);      q = qdd_gate(q, GATEID_H, 16);          q = qdd_gate(q, GATEID_Z, 8);
    q = qdd_cgate(q, GATEID_X, 12, 14);     q = qdd_cgate(q, GATEID_Z, 2, 18);      q = qdd_cgate(q, GATEID_X, 12, 15);     q = qdd_gate(q, GATEID_X, 9);
    q = qdd_gate(q, GATEID_Z, 12);          q = qdd_gate(q, GATEID_X, 3);           q = qdd_gate(q, GATEID_X, 0);           q = qdd_cgate(q, GATEID_X, 1, 4);
    q = qdd_gate(q, GATEID_H, 1);           q = qdd_gate(q, GATEID_X, 19);          q = qdd_gate(q, GATEID_X, 5);           q = qdd_cgate(q, GATEID_Z, 2, 16);
    q = qdd_gate(q, GATEID_X, 4);           q = qdd_cgate(q, GATEID_X, 9, 11);      q = qdd_cgate(q, GATEID_X, 0, 7);       q = qdd_gate(q, GATEID_Z, 12);
    q = qdd_cgate(q, GATEID_X, 9, 11);      q = qdd_gate(q, GATEID_Z, 13);          q = qdd_cgate(q, GATEID_X, 12, 16);     q = qdd_gate(q, GATEID_Z, 10);
    q = qdd_gate(q, GATEID_X, 4);           q = qdd_gate(q, GATEID_Z, 16);          q = qdd_cgate(q, GATEID_Z, 4, 17);      q = qdd_gate(q, GATEID_Z, 7);
    q = qdd_gate(q, GATEID_H, 4);           q = qdd_cgate(q, GATEID_Z, 6, 7);       q = qdd_cgate(q, GATEID_X, 12, 19);     q = qdd_gate(q, GATEID_Z, 15);
    q = qdd_cgate(q, GATEID_X, 5, 11);      q = qdd_cgate(q, GATEID_X, 9, 17);      q = qdd_gate(q, GATEID_Z, 3);           q = qdd_cgate(q, GATEID_X, 11, 18);
    q = qdd_cgate(q, GATEID_Z, 5, 15);      q = qdd_cgate(q, GATEID_X, 0, 15);      q = qdd_cgate(q, GATEID_X, 1, 6);       q = qdd_cgate(q, GATEID_X, 8, 16);
    q = qdd_cgate(q, GATEID_X, 5, 19);      q = qdd_cgate(q, GATEID_Z, 3, 18);      q = qdd_cgate(q, GATEID_X, 5, 8);       q = qdd_cgate(q, GATEID_Z, 14, 18);
    node_count = qdd_countnodes(q);
    if (test_measure_random_state(q, 20)) return 1;

    // inverse
    q = qdd_cgate(q, GATEID_Z, 14, 18);     q = qdd_cgate(q, GATEID_X, 5, 8);       q = qdd_cgate(q, GATEID_Z, 3, 18);      q = qdd_cgate(q, GATEID_X, 5, 19);
    q = qdd_cgate(q, GATEID_X, 8, 16);      q = qdd_cgate(q, GATEID_X, 1, 6);       q = qdd_cgate(q, GATEID_X, 0, 15);      q = qdd_cgate(q, GATEID_Z, 5, 15);
    q = qdd_cgate(q, GATEID_X, 11, 18);     q = qdd_gate(q, GATEID_Z, 3);           q = qdd_cgate(q, GATEID_X, 9, 17);      q = qdd_cgate(q, GATEID_X, 5, 11);
    q = qdd_gate(q, GATEID_Z, 15);          q = qdd_cgate(q, GATEID_X, 12, 19);     q = qdd_cgate(q, GATEID_Z, 6, 7);       q = qdd_gate(q, GATEID_H, 4);
    q = qdd_gate(q, GATEID_Z, 7);           q = qdd_cgate(q, GATEID_Z, 4, 17);      q = qdd_gate(q, GATEID_Z, 16);          q = qdd_gate(q, GATEID_X, 4);
    q = qdd_gate(q, GATEID_Z, 10);          q = qdd_cgate(q, GATEID_X, 12, 16);     q = qdd_gate(q, GATEID_Z, 13);          q = qdd_cgate(q, GATEID_X, 9, 11);
    q = qdd_gate(q, GATEID_Z, 12);          q = qdd_cgate(q, GATEID_X, 0, 7);       q = qdd_cgate(q, GATEID_X, 9, 11);      q = qdd_gate(q, GATEID_X, 4);
    q = qdd_cgate(q, GATEID_Z, 2, 16);      q = qdd_gate(q, GATEID_X, 5);           q = qdd_gate(q, GATEID_X, 19);          q = qdd_gate(q, GATEID_H, 1);
    q = qdd_cgate(q, GATEID_X, 1, 4);       q = qdd_gate(q, GATEID_X, 0);           q = qdd_gate(q, GATEID_X, 3);           q = qdd_gate(q, GATEID_Z, 12);
    q = qdd_gate(q, GATEID_X, 9);           q = qdd_cgate(q, GATEID_X, 12, 15);     q = qdd_cgate(q, GATEID_Z, 2, 18);      q = qdd_cgate(q, GATEID_X, 12, 14);
    q = qdd_gate(q, GATEID_Z, 8);           q = qdd_gate(q, GATEID_H, 16);          q = qdd_cgate(q, GATEID_X, 0, 15);      q = qdd_gate(q, GATEID_Z, 0);
    q = qdd_gate(q, GATEID_X, 15);          q = qdd_gate(q, GATEID_X, 6);           q = qdd_cgate(q, GATEID_Z, 1, 3);       q = qdd_cgate(q, GATEID_X, 13, 19);
    q = qdd_cgate(q, GATEID_Z, 11, 16);     q = qdd_cgate(q, GATEID_Z, 10, 12);     q = qdd_cgate(q, GATEID_X, 8, 9);       q = qdd_cgate(q, GATEID_Z, 2, 19);
    q = qdd_gate(q, GATEID_X, 1);           q = qdd_cgate(q, GATEID_Z, 4, 5);       q = qdd_gate(q, GATEID_Z, 11);          q = qdd_gate(q, GATEID_X, 18);
    q = qdd_gate(q, GATEID_X, 6);           q = qdd_cgate(q, GATEID_X, 9, 10);      q = qdd_gate(q, GATEID_X, 2);           q = qdd_gate(q, GATEID_X, 16);
    q = qdd_cgate(q, GATEID_X, 7, 18);      q = qdd_gate(q, GATEID_Z, 16);          q = qdd_gate(q, GATEID_X, 18);          q = qdd_cgate(q, GATEID_X, 7, 10);
    q = qdd_cgate(q, GATEID_Z, 0, 15);      q = qdd_cgate(q, GATEID_Z, 16, 18);     q = qdd_gate(q, GATEID_X, 15);          q = qdd_cgate(q, GATEID_X, 4, 10);
    q = qdd_gate(q, GATEID_X, 18);          q = qdd_gate(q, GATEID_Z, 4);           q = qdd_gate(q, GATEID_X, 18);          q = qdd_cgate(q, GATEID_X, 12, 16);
    q = qdd_gate(q, GATEID_Z, 11);          q = qdd_gate(q, GATEID_X, 17);          q = qdd_cgate(q, GATEID_X, 7, 18);      q = qdd_gate(q, GATEID_X, 8);
    q = qdd_gate(q, GATEID_X, 9);           q = qdd_gate(q, GATEID_H, 8);           q = qdd_cgate(q, GATEID_Z, 10, 12);     q = qdd_cgate(q, GATEID_Z, 14, 15);
    q = qdd_gate(q, GATEID_Z, 18);          q = qdd_cgate(q, GATEID_X, 5, 19);      q = qdd_cgate(q, GATEID_Z, 11, 16);     q = qdd_cgate(q, GATEID_X, 3, 6);
    q = qdd_gate(q, GATEID_Z, 6);           q = qdd_cgate(q, GATEID_X, 7, 10);      q = qdd_cgate(q, GATEID_X, 2, 15);      q = qdd_cgate(q, GATEID_X, 3, 18);
    q = qdd_cgate(q, GATEID_X, 6, 16);      q = qdd_cgate(q, GATEID_Z, 1, 4);       q = qdd_cgate(q, GATEID_X, 3, 5);       q = qdd_cgate(q, GATEID_Z, 7, 11);
    q = qdd_cgate(q, GATEID_Z, 13, 16);     q = qdd_cgate(q, GATEID_Z, 1, 16);      q = qdd_cgate(q, GATEID_Z, 0, 4);       q = qdd_cgate(q, GATEID_Z, 9, 10);
    q = qdd_gate(q, GATEID_Z, 4);           q = qdd_cgate(q, GATEID_X, 1, 12);      q = qdd_gate(q, GATEID_H, 16);          q = qdd_cgate(q, GATEID_Z, 4, 18);

    test_assert(qdd_equivalent(q, qref, 20, false, VERBOSE)); // check approx equiv
    test_assert(qdd_equivalent(q, qref, 20, true,  VERBOSE)); // check exact equiv
    test_assert(q == qref);

    if(VERBOSE) printf("qdd 20 qubit circuit:     ok (%ld nodes)\n", node_count);
    return 0;
}

int run_qdd_tests()
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

int test_with(int amps_backend, int norm_strat) 
{
    // Standard Lace initialization
    int workers = 1;
    lace_init(workers, 0);
    printf("%d worker(s), ", workers);
    lace_startup(0, NULL, NULL);

    // Simple Sylvan initialization
    sylvan_set_sizes(1LL<<25, 1LL<<25, 1LL<<16, 1LL<<16);
    sylvan_init_package();
    sylvan_init_qdd(1LL<<11, -1, amps_backend, norm_strat);
    qdd_set_testing_mode(true); // turn on internal sanity tests

    printf("amps backend = %d, norm strategy = %d:\n", amps_backend, norm_strat);
    int res = run_qdd_tests();

    sylvan_quit();
    lace_exit();

    return res;
}

int runtests()
{
    for (int backend = 0; backend < n_backends; backend++) {
        for (int norm_strat = 0; norm_strat < n_norm_stragegies; norm_strat++) {
            if (test_with(backend, norm_strat)) return 1;
        }
    }
    return 0;
}

int main()
{
    return runtests();
}
