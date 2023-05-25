#include <qsylvan_gates.h>
#include <sylvan_int.h>


static long double Pi;    // set value of global Pi

uint64_t gates[n_predef_gates+256+256+1000][4];
gate_id_t inv_gate_ids[n_predef_gates+256+256+1000];

/********************** <dynamic custom rotation gates> ***********************/

uint32_t next_custom_id; // set to 0 in init

uint32_t
get_custom_gate_id()
{
    if (next_custom_id + 2 >= num_dynamic_gates) {
        // max custom gates used, reset ID counter to 0 and clear opcache
        next_custom_id = 0;
        sylvan_clear_cache();
    }

    uint32_t res = next_custom_id + 1;
    next_custom_id += 2; // reserve next id for inverse of returned gateid

    return num_static_gates + res; // index offset by num_static_gates
}

uint32_t
GATEID_Rz(fl_t a)
{
    // get gate id for this gate
    uint32_t gate_id = get_custom_gate_id();

    // initialize gate
    double theta_over_2 = Pi * a;
    AMP u00, u11;
    u00 = weight_lookup(cmake_angle(-theta_over_2, 1));
    u11 = weight_lookup(cmake_angle(theta_over_2, 1));
    gates[gate_id][0] = u00;        gates[gate_id][1] = AADD_ZERO;
    gates[gate_id][2] = AADD_ZERO;  gates[gate_id][3] = u11;

    // also store the inverse of this gate (gate_id + 1 is reserved for this)
    gates[gate_id+1][0] = wgt_ccj(u00); gates[gate_id+1][1] = AADD_ZERO;
    gates[gate_id+1][2] = AADD_ZERO;    gates[gate_id+1][3] = wgt_ccj(u11);
    inv_gate_ids[gate_id] = gate_id+1;
    inv_gate_ids[gate_id+1] = gate_id;

    // return (temporary) gate_id for this gate
    return gate_id;
}

uint32_t
GATEID_Rx(fl_t a)
{
    // get gate id for this gate
    uint32_t gate_id = get_custom_gate_id();

    // initialize gate
    fl_t theta_over_2 = Pi * a;
    AMP u00, u01, u10, u11;
    u00 = weight_lookup(cmake(flt_cos(theta_over_2), 0.0));
    u01 = weight_lookup(cmake(0.0, -flt_sin(theta_over_2)));
    u10 = weight_lookup(cmake(0.0, -flt_sin(theta_over_2)));
    u11 = weight_lookup(cmake(flt_cos(theta_over_2), 0.0));
    gates[gate_id][0] = u00; gates[gate_id][1] = u01;
    gates[gate_id][2] = u10; gates[gate_id][3] = u11;

    // also store the inverse of this gate (gate_id + 1 is reserved for this)
    gates[gate_id+1][0] = wgt_ccj(u00); gates[gate_id+1][1] = wgt_ccj(u10);
    gates[gate_id+1][2] = wgt_ccj(u01); gates[gate_id+1][3] = wgt_ccj(u11);
    inv_gate_ids[gate_id] = gate_id+1;
    inv_gate_ids[gate_id+1] = gate_id;

    // return (temporary) gate_id for this gate
    return gate_id;
}

uint32_t
GATEID_Ry(fl_t a)
{
    // get gate id for this gate
    uint32_t gate_id = get_custom_gate_id();

    // initialize gate
    fl_t theta_over_2 = Pi * a;
    AMP u00, u01, u10, u11;
    u00 = weight_lookup(cmake(flt_cos(theta_over_2),  0.0));
    u01 = weight_lookup(cmake(-flt_sin(theta_over_2), 0.0));
    u10 = weight_lookup(cmake(flt_sin(theta_over_2),  0.0));
    u11 = weight_lookup(cmake(flt_cos(theta_over_2),  0.0));
    gates[gate_id][0] = u00; gates[gate_id][1] = u01;
    gates[gate_id][2] = u10; gates[gate_id][3] = u11;

    // also store the inverse of this gate (gate_id + 1 is reserved for this)
    gates[gate_id+1][0] = wgt_ccj(u00); gates[gate_id+1][1] = wgt_ccj(u10);
    gates[gate_id+1][2] = wgt_ccj(u01); gates[gate_id+1][3] = wgt_ccj(u11);
    inv_gate_ids[gate_id] = gate_id+1;
    inv_gate_ids[gate_id+1] = gate_id;

    // return (temporary) gate_id for this gate
    return gate_id;
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
    inv_gate_ids[GATEID_I] = GATEID_I;

    k = GATEID_X;
    gates[k][0] = AADD_ZERO; gates[k][1] = AADD_ONE;
    gates[k][2] = AADD_ONE;  gates[k][3] = AADD_ZERO;
    inv_gate_ids[GATEID_X] = GATEID_X;

    k = GATEID_Y;
    gates[k][0] = AADD_ZERO; gates[k][1] = weight_lookup(cmake(0.0, -1.0));
    gates[k][2] = weight_lookup(cmake(0.0, 1.0));  gates[k][3] = AADD_ZERO;
    inv_gate_ids[GATEID_Y] = GATEID_Y;

    k = GATEID_Z;
    gates[k][0] = AADD_ONE;  gates[k][1] = AADD_ZERO;
    gates[k][2] = AADD_ZERO; gates[k][3] = AADD_MIN_ONE;
    inv_gate_ids[GATEID_Z] = GATEID_Z;

    k = GATEID_H;
    gates[k][0] = gates[k][1] = gates[k][2] = weight_lookup(cmake(1.0/flt_sqrt(2.0),0));
    gates[k][3] = weight_lookup(cmake(-1.0/flt_sqrt(2.0),0));
    inv_gate_ids[GATEID_H] = GATEID_H;

    k = GATEID_S;
    gates[k][0] = AADD_ONE;  gates[k][1] = AADD_ZERO;
    gates[k][2] = AADD_ZERO; gates[k][3] = weight_lookup(cmake(0.0, 1.0));
    inv_gate_ids[GATEID_S] = GATEID_Sdag;

    k = GATEID_Sdag;
    gates[k][0] = AADD_ONE;  gates[k][1] = AADD_ZERO;
    gates[k][2] = AADD_ZERO; gates[k][3] = weight_lookup(cmake(0.0, -1.0));
    inv_gate_ids[GATEID_Sdag] = GATEID_S;

    k = GATEID_T;
    gates[k][0] = AADD_ONE;  gates[k][1] = AADD_ZERO;
    gates[k][2] = AADD_ZERO; gates[k][3] = weight_lookup(cmake(1.0/flt_sqrt(2.0), 1.0/flt_sqrt(2.0)));
    inv_gate_ids[GATEID_T] = GATEID_Tdag;

    k = GATEID_Tdag;
    gates[k][0] = AADD_ONE;  gates[k][1] = AADD_ZERO;
    gates[k][2] = AADD_ZERO; gates[k][3] = weight_lookup(cmake(1.0/flt_sqrt(2.0), -1.0/flt_sqrt(2.0)));
    inv_gate_ids[GATEID_Tdag] = GATEID_T;

    k = GATEID_sqrtX;
    gates[k][0] = weight_lookup(cmake(0.5, 0.5)); gates[k][1] = weight_lookup(cmake(0.5,-0.5));
    gates[k][2] = weight_lookup(cmake(0.5,-0.5)); gates[k][3] = weight_lookup(cmake(0.5, 0.5));
    inv_gate_ids[GATEID_sqrtX] = GATEID_sqrtXdag;

    k = GATEID_sqrtXdag;
    gates[k][0] = weight_lookup(cmake(0.5,-0.5)); gates[k][1] = weight_lookup(cmake(0.5, 0.5));
    gates[k][2] = weight_lookup(cmake(0.5, 0.5)); gates[k][3] = weight_lookup(cmake(0.5,-0.5));
    inv_gate_ids[GATEID_sqrtXdag] = GATEID_sqrtX;

    k = GATEID_sqrtY;
    gates[k][0] = weight_lookup(cmake(0.5, 0.5)); gates[k][1] = weight_lookup(cmake(-0.5,-0.5));
    gates[k][2] = weight_lookup(cmake(0.5, 0.5)); gates[k][3] = weight_lookup(cmake(0.5, 0.5));
    inv_gate_ids[GATEID_sqrtY] = GATEID_sqrtYdag;

    k = GATEID_sqrtYdag;
    gates[k][0] = weight_lookup(cmake(0.5,-0.5)); gates[k][1] = weight_lookup(cmake(0.5,-0.5));
    gates[k][2] = weight_lookup(cmake(-0.5,0.5)); gates[k][3] = weight_lookup(cmake(0.5,-0.5));
    inv_gate_ids[GATEID_sqrtYdag] = GATEID_sqrtY;

    qmdd_phase_gates_init(255);

    next_custom_id = 0;
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
        gates[gate_id][2] = AADD_ZERO; gates[gate_id][3] = weight_lookup(cartesian);

        // backward rotation
        angle = -2*Pi / (fl_t)(1<<k);
        cartesian = cmake_angle(angle, 1);
        gate_id = GATEID_Rk_dag(k);
        gates[gate_id][0] = AADD_ONE;  gates[gate_id][1] = AADD_ZERO;
        gates[gate_id][2] = AADD_ZERO; gates[gate_id][3] = weight_lookup(cartesian);

        inv_gate_ids[GATEID_Rk(k)] = GATEID_Rk_dag(k);
        inv_gate_ids[GATEID_Rk_dag(k)] = GATEID_Rk(k);
    }
}