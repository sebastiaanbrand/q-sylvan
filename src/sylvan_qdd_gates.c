#include "sylvan_qdd_gates.h"
#include "sylvan_qdd_complex.h"
#include "sylvan_int.h"


static long double Pi;    // set value of global Pi


/********************** <dynamic custom rotation gates> ***********************/

uint32_t next_custom_id; // set to 0 in init

uint32_t
get_custom_gate_id()
{
    next_custom_id++;
    if (next_custom_id >= num_dynamic_gates) {
        // max custom gates used, reset ID counter to 0 and clear opcache
        next_custom_id = 0;
        LACE_ME;
        sylvan_clear_cache();
    }
    return num_static_gates + next_custom_id; // index offset by num_static_gates
}

uint32_t
GATEID_Rz(fl_t a)
{
    // get gate id for this gate
    uint32_t gate_id = get_custom_gate_id();

    // initialize gate
    double theta_over_2 = Pi * a;
    AMP u00, u11;
    u00 = comp_lookup(comp_make_angle(-theta_over_2));
    u11 = comp_lookup(comp_make_angle(theta_over_2));
    gates[gate_id][0] = u00;    gates[gate_id][1] = C_ZERO;
    gates[gate_id][2] = C_ZERO; gates[gate_id][3] = u11;

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
    u00 = comp_lookup(comp_make(flt_cos(theta_over_2), 0.0));
    u01 = comp_lookup(comp_make(0.0, -flt_sin(theta_over_2)));
    u10 = comp_lookup(comp_make(0.0, -flt_sin(theta_over_2)));
    u11 = comp_lookup(comp_make(flt_cos(theta_over_2), 0.0));
    gates[gate_id][0] = u00; gates[gate_id][1] = u01;
    gates[gate_id][2] = u10; gates[gate_id][3] = u11;

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
    u00 = comp_lookup(comp_make(flt_cos(theta_over_2),  0.0));
    u01 = comp_lookup(comp_make(-flt_sin(theta_over_2), 0.0));
    u10 = comp_lookup(comp_make(flt_sin(theta_over_2),  0.0));
    u11 = comp_lookup(comp_make(flt_cos(theta_over_2),  0.0));
    gates[gate_id][0] = u00; gates[gate_id][1] = u01;
    gates[gate_id][2] = u10; gates[gate_id][3] = u11;

    // return (temporary) gate_id for this gate
    return gate_id;
}

/********************* </dynamic custom rotation gates> ***********************/


/*************************** <dynamic custom gates> ***************************/
void
qdd_gates_init()
{
    Pi = 2.0 * flt_acos(0.0);

    // initialize 2x2 gates (complex values from gates currently stored in 
    // same table as complex amplitude values)
    uint32_t k;

    k = GATEID_I;
    gates[k][0] = C_ONE;  gates[k][1] = C_ZERO;
    gates[k][2] = C_ZERO; gates[k][3] = C_ONE;

    k = GATEID_X;
    gates[k][0] = C_ZERO; gates[k][1] = C_ONE;
    gates[k][2] = C_ONE;  gates[k][3] = C_ZERO;

    k = GATEID_Y;
    gates[k][0] = C_ZERO; gates[k][1] = comp_lookup(comp_make(0.0, -1.0));
    gates[k][2] = comp_lookup(comp_make(0.0, 1.0));  gates[k][3] = C_ZERO;

    k = GATEID_Z;
    gates[k][0] = C_ONE;  gates[k][1] = C_ZERO;
    gates[k][2] = C_ZERO; gates[k][3] = C_MIN_ONE;

    k = GATEID_H;
    gates[k][0] = gates[k][1] = gates[k][2] = comp_lookup(comp_make(1.0/flt_sqrt(2.0),0));
    gates[k][3] = comp_lookup(comp_make(-1.0/flt_sqrt(2.0),0));

    k = GATEID_S;
    gates[k][0] = C_ONE;  gates[k][1] = C_ZERO;
    gates[k][2] = C_ZERO; gates[k][3] = comp_lookup(comp_make(0.0, 1.0));

    k = GATEID_Sdag;
    gates[k][0] = C_ONE;  gates[k][1] = C_ZERO;
    gates[k][2] = C_ZERO; gates[k][3] = comp_lookup(comp_make(0.0, -1.0));

    k = GATEID_T;
    gates[k][0] = C_ONE;  gates[k][1] = C_ZERO;
    gates[k][2] = C_ZERO; gates[k][3] = comp_lookup(comp_make(1.0/flt_sqrt(2.0), 1.0/flt_sqrt(2.0)));

    k = GATEID_Tdag;
    gates[k][0] = C_ONE;  gates[k][1] = C_ZERO;
    gates[k][2] = C_ZERO; gates[k][3] = comp_lookup(comp_make(1.0/flt_sqrt(2.0), -1.0/flt_sqrt(2.0)));

    k = GATEID_sqrtX;
    gates[k][0] = comp_lookup(comp_make(0.5, 0.5)); gates[k][1] = comp_lookup(comp_make(0.5,-0.5));
    gates[k][2] = comp_lookup(comp_make(0.5,-0.5)); gates[k][3] = comp_lookup(comp_make(0.5, 0.5));

    k = GATEID_sqrtXdag;
    gates[k][0] = comp_lookup(comp_make(0.5,-0.5)); gates[k][1] = comp_lookup(comp_make(0.5, 0.5));
    gates[k][2] = comp_lookup(comp_make(0.5, 0.5)); gates[k][3] = comp_lookup(comp_make(0.5,-0.5));

    k = GATEID_sqrtY;
    gates[k][0] = comp_lookup(comp_make(0.5, 0.5)); gates[k][1] = comp_lookup(comp_make(-0.5,-0.5));
    gates[k][2] = comp_lookup(comp_make(0.5, 0.5)); gates[k][3] = comp_lookup(comp_make(0.5, 0.5));

    k = GATEID_sqrtYdag;
    gates[k][0] = comp_lookup(comp_make(0.5,-0.5)); gates[k][1] = comp_lookup(comp_make(0.5,-0.5));
    gates[k][2] = comp_lookup(comp_make(-0.5,0.5)); gates[k][3] = comp_lookup(comp_make(0.5,-0.5));

    qdd_phase_gates_init(255);

    next_custom_id = 0;
}

void
qdd_phase_gates_init(int n)
{
    // add gate R_k to gates table
    // (note that R_0 = I, R_1 = Z, R_2 = S, R_4 = T)
    uint32_t gate_id;
    fl_t angle;
    complex_t cartesian;
    for (int k=0; k<=n; k++) {
        // forward rotation
        angle = 2*Pi / (fl_t)(1<<k);
        cartesian = comp_make_angle(angle);
        gate_id = GATEID_Rk(k);
        gates[gate_id][0] = C_ONE;  gates[gate_id][1] = C_ZERO;
        gates[gate_id][2] = C_ZERO; gates[gate_id][3] = comp_lookup(cartesian);

        // backward rotation
        angle = -2*Pi / (fl_t)(1<<k);
        cartesian = comp_make_angle(angle);
        gate_id = GATEID_Rk_dag(k);
        gates[gate_id][0] = C_ONE;  gates[gate_id][1] = C_ZERO;
        gates[gate_id][2] = C_ZERO; gates[gate_id][3] = comp_lookup(cartesian);
    }
}