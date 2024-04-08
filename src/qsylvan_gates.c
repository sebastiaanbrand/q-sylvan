#include <qsylvan_gates.h>
#include <sylvan_int.h>
#include <sylvan_edge_weights_complex.h>


static long double Pi;    // set value of global Pi

uint64_t gates[n_predef_gates+256+256][4];

/********************** <dynamic custom rotation gates> ***********************/


// store complex values of dynamic gate to re-initialize gate after gc
complex_t dynamic_gate[4];

uint32_t
GATEID_Rz(fl_t theta)
{
    // clear cache to invalidate cached results for GATEID_dynamic
    sylvan_clear_cache();

    // initialize (and store for gc)
    dynamic_gate[0] = cmake_angle(-theta/2.0, 1);
    dynamic_gate[1] = czero();
    dynamic_gate[2] = czero();
    dynamic_gate[3] = cmake_angle(theta/2.0, 1);
    gates[GATEID_dynamic][0] = weight_lookup(&dynamic_gate[0]); // u00
    gates[GATEID_dynamic][1] = weight_lookup(&dynamic_gate[1]); // u01
    gates[GATEID_dynamic][2] = weight_lookup(&dynamic_gate[2]); // u10
    gates[GATEID_dynamic][3] = weight_lookup(&dynamic_gate[3]); // u11

    // return (temporary) gate_id for this gate
    return GATEID_dynamic;
}


uint32_t
GATEID_Rx(fl_t theta)
{
    // clear cache to invalidate cached results for GATEID_dynamic
    sylvan_clear_cache();

    // initialize (and store for gc)
    dynamic_gate[0] = cmake(flt_cos(theta/2.0), 0.0);
    dynamic_gate[1] = cmake(0.0, -flt_sin(theta/2.0));
    dynamic_gate[2] = cmake(0.0, -flt_sin(theta/2.0));
    dynamic_gate[3] = cmake(flt_cos(theta/2.0), 0.0);
    gates[GATEID_dynamic][0] = weight_lookup(&dynamic_gate[0]); // u00
    gates[GATEID_dynamic][1] = weight_lookup(&dynamic_gate[1]); // u01
    gates[GATEID_dynamic][2] = weight_lookup(&dynamic_gate[2]); // u10
    gates[GATEID_dynamic][3] = weight_lookup(&dynamic_gate[3]); // u11

    // return (temporary) gate_id for this gate
    return GATEID_dynamic;
}

uint32_t
GATEID_Ry(fl_t theta)
{
    // clear cache to invalidate cached results for GATEID_dynamic
    sylvan_clear_cache();

    // initialize (and store for gc)
    dynamic_gate[0] = cmake( flt_cos(theta/2.0), 0.0);
    dynamic_gate[1] = cmake(-flt_sin(theta/2.0), 0.0);
    dynamic_gate[2] = cmake( flt_sin(theta/2.0), 0.0);
    dynamic_gate[3] = cmake( flt_cos(theta/2.0), 0.0);
    gates[GATEID_dynamic][0] = weight_lookup(&dynamic_gate[0]); // u00
    gates[GATEID_dynamic][1] = weight_lookup(&dynamic_gate[1]); // u01
    gates[GATEID_dynamic][2] = weight_lookup(&dynamic_gate[2]); // u10
    gates[GATEID_dynamic][3] = weight_lookup(&dynamic_gate[3]); // u11

    // return (temporary) gate_id for this gate
    return GATEID_dynamic;
}

uint32_t
GATEID_Phase(fl_t theta)
{
    // clear cache to invalidate cached results for GATEID_dynamic
    sylvan_clear_cache();
    
    // initialize (and store for gc)
    dynamic_gate[0] = cmake(1.0, 0.0);
    dynamic_gate[1] = cmake(0.0, 0.0);
    dynamic_gate[2] = cmake(0.0, 0.0);
    dynamic_gate[3] = cmake_angle(theta, 1);
    gates[GATEID_dynamic][0] = weight_lookup(&dynamic_gate[0]); // u00
    gates[GATEID_dynamic][1] = weight_lookup(&dynamic_gate[1]); // u01
    gates[GATEID_dynamic][2] = weight_lookup(&dynamic_gate[2]); // u10
    gates[GATEID_dynamic][3] = weight_lookup(&dynamic_gate[3]); // u11

    // return (temporary) gate_id for this gate
    return GATEID_dynamic;
}

uint32_t
GATEID_U(fl_t theta, fl_t phi, fl_t lambda)
{
    // clear cache to invalidate cached results for GATEID_dynamic
    sylvan_clear_cache();

    // initialize (and store for gc)
    dynamic_gate[0] = cmake(flt_cos(theta/2.0), 0.0);
    dynamic_gate[1] = cmul(cmake_angle(lambda,1), cmake(-flt_sin(theta/2.0), 0));
    dynamic_gate[2] = cmul(cmake_angle(phi,1), cmake(flt_sin(theta/2.0), 0));
    dynamic_gate[3] = cmul(cmake_angle(phi+lambda,1), cmake(flt_cos(theta/2.0), 0));
    gates[GATEID_dynamic][0] = weight_lookup(&dynamic_gate[0]); // u00
    gates[GATEID_dynamic][1] = weight_lookup(&dynamic_gate[1]); // u01
    gates[GATEID_dynamic][2] = weight_lookup(&dynamic_gate[2]); // u10
    gates[GATEID_dynamic][3] = weight_lookup(&dynamic_gate[3]); // u11

    // return (temporary) gate_id for this gate
    return GATEID_dynamic;
}

/********************* </dynamic custom rotation gates> ***********************/


/*************************** <dynamic custom gates> ***************************/
void
qmdd_gates_init()
{
    Pi = 2.0 * flt_acos(0.0);

    // initialize 2x2 gates (complex values from gates currently stored in 
    // same table as complex amplitude values)
    uint32_t k;

    k = GATEID_I;
    gates[k][0] = AADD_ONE;  gates[k][1] = AADD_ZERO;
    gates[k][2] = AADD_ZERO; gates[k][3] = AADD_ONE;

    k = GATEID_X;
    gates[k][0] = AADD_ZERO; gates[k][1] = AADD_ONE;
    gates[k][2] = AADD_ONE;  gates[k][3] = AADD_ZERO;

    k = GATEID_Y;
    gates[k][0] = AADD_ZERO; gates[k][1] = complex_lookup(0.0, -1.0);
    gates[k][2] = complex_lookup(0.0, 1.0);  gates[k][3] = AADD_ZERO;

    k = GATEID_Z;
    gates[k][0] = AADD_ONE;  gates[k][1] = AADD_ZERO;
    gates[k][2] = AADD_ZERO; gates[k][3] = AADD_MIN_ONE;

    k = GATEID_H;
    gates[k][0] = gates[k][1] = gates[k][2] = complex_lookup(1.0/flt_sqrt(2.0),0);
    gates[k][3] = complex_lookup(-1.0/flt_sqrt(2.0),0);

    k = GATEID_S;
    gates[k][0] = AADD_ONE;  gates[k][1] = AADD_ZERO;
    gates[k][2] = AADD_ZERO; gates[k][3] = complex_lookup(0.0, 1.0);

    k = GATEID_Sdag;
    gates[k][0] = AADD_ONE;  gates[k][1] = AADD_ZERO;
    gates[k][2] = AADD_ZERO; gates[k][3] = complex_lookup(0.0, -1.0);

    k = GATEID_T;
    gates[k][0] = AADD_ONE;  gates[k][1] = AADD_ZERO;
    gates[k][2] = AADD_ZERO; gates[k][3] = complex_lookup(1.0/flt_sqrt(2.0), 1.0/flt_sqrt(2.0));

    k = GATEID_Tdag;
    gates[k][0] = AADD_ONE;  gates[k][1] = AADD_ZERO;
    gates[k][2] = AADD_ZERO; gates[k][3] = complex_lookup(1.0/flt_sqrt(2.0), -1.0/flt_sqrt(2.0));

    k = GATEID_sqrtX;
    gates[k][0] = complex_lookup(0.5, 0.5); gates[k][1] = complex_lookup(0.5,-0.5);
    gates[k][2] = complex_lookup(0.5,-0.5); gates[k][3] = complex_lookup(0.5, 0.5);

    k = GATEID_sqrtXdag;
    gates[k][0] = complex_lookup(0.5,-0.5); gates[k][1] = complex_lookup(0.5, 0.5);
    gates[k][2] = complex_lookup(0.5, 0.5); gates[k][3] = complex_lookup(0.5,-0.5);

    k = GATEID_sqrtY;
    gates[k][0] = complex_lookup(0.5, 0.5); gates[k][1] = complex_lookup(-0.5,-0.5);
    gates[k][2] = complex_lookup(0.5, 0.5); gates[k][3] = complex_lookup(0.5, 0.5);

    k = GATEID_sqrtYdag;
    gates[k][0] = complex_lookup(0.5,-0.5); gates[k][1] = complex_lookup(0.5,-0.5);
    gates[k][2] = complex_lookup(-0.5,0.5); gates[k][3] = complex_lookup(0.5,-0.5);

    qmdd_phase_gates_init(255);

    // init dynamic gate 
    // (necessary when qmdd_gates_init() is called after gc to re-init all gates)
    k = GATEID_dynamic;
    gates[k][0] = weight_lookup(&dynamic_gate[0]);
    gates[k][1] = weight_lookup(&dynamic_gate[1]);
    gates[k][2] = weight_lookup(&dynamic_gate[2]);
    gates[k][3] = weight_lookup(&dynamic_gate[3]);
}

void
qmdd_phase_gates_init(int n)
{
    // add gate R_k to gates table
    // (note that R_0 = I, R_1 = Z, R_2 = S, R_4 = T)
    uint32_t gate_id;
    fl_t angle;
    complex_t cartesian;
    for (int k=0; k<=n; k++) {
        // forward rotation
        angle = 2*Pi / (fl_t)(1<<k);
        cartesian = cmake_angle(angle, 1);
        gate_id = GATEID_Rk(k);
        gates[gate_id][0] = AADD_ONE;  gates[gate_id][1] = AADD_ZERO;
        gates[gate_id][2] = AADD_ZERO; gates[gate_id][3] = weight_lookup(&cartesian);

        // backward rotation
        angle = -2*Pi / (fl_t)(1<<k);
        cartesian = cmake_angle(angle, 1);
        gate_id = GATEID_Rk_dag(k);
        gates[gate_id][0] = AADD_ONE;  gates[gate_id][1] = AADD_ZERO;
        gates[gate_id][2] = AADD_ZERO; gates[gate_id][3] = weight_lookup(&cartesian);
    }
}
