#include <qsylvan_gates.h>
#include <sylvan_int.h>


static long double Pi;    // set value of global Pi

uint64_t gates[n_predef_gates+256+256+1000][4];

/********************** <dynamic custom rotation gates> ***********************/

uint32_t next_custom_id; // set to 0 in init

uint32_t
get_custom_gate_id()
{
    next_custom_id++;
    if (next_custom_id >= num_dynamic_gates) {
        // max custom gates used, reset ID counter to 0 and clear opcache
        next_custom_id = 0;
        sylvan_clear_cache();
    }
    return num_static_gates + next_custom_id; // index offset by num_static_gates
}

uint32_t
GATEID_Rz(fl_t theta)
{
    // get gate id for this gate
    uint32_t gate_id = get_custom_gate_id();

    // initialize gate
    AMP u00, u11;
    u00 = weight_lookup(cmake_angle(-theta/2.0, 1));
    u11 = weight_lookup(cmake_angle(theta/2.0, 1));
    gates[gate_id][0] = u00;    gates[gate_id][1] = AADD_ZERO;
    gates[gate_id][2] = AADD_ZERO; gates[gate_id][3] = u11;

    // return (temporary) gate_id for this gate
    return gate_id;
}

// TODO: add GATEID_Phase(fl_t a) (global phase diff with Rz)

uint32_t
GATEID_Rx(fl_t theta)
{
    // get gate id for this gate
    uint32_t gate_id = get_custom_gate_id();

    // initialize gate
    AMP u00, u01, u10, u11;
    u00 = weight_lookup(cmake(flt_cos(theta/2.0), 0.0));
    u01 = weight_lookup(cmake(0.0, -flt_sin(theta/2.0)));
    u10 = weight_lookup(cmake(0.0, -flt_sin(theta/2.0)));
    u11 = weight_lookup(cmake(flt_cos(theta/2.0), 0.0));
    gates[gate_id][0] = u00; gates[gate_id][1] = u01;
    gates[gate_id][2] = u10; gates[gate_id][3] = u11;

    // return (temporary) gate_id for this gate
    return gate_id;
}

uint32_t
GATEID_Ry(fl_t theta)
{
    // get gate id for this gate
    uint32_t gate_id = get_custom_gate_id();

    // initialize gate
    AMP u00, u01, u10, u11;
    u00 = weight_lookup(cmake(flt_cos(theta/2.0),  0.0));
    u01 = weight_lookup(cmake(-flt_sin(theta/2.0), 0.0));
    u10 = weight_lookup(cmake(flt_sin(theta/2.0),  0.0));
    u11 = weight_lookup(cmake(flt_cos(theta/2.0),  0.0));
    gates[gate_id][0] = u00; gates[gate_id][1] = u01;
    gates[gate_id][2] = u10; gates[gate_id][3] = u11;

    // return (temporary) gate_id for this gate
    return gate_id;
}

uint32_t
GATEID_Phase(fl_t theta)
{
    // get gate id for this gate
    uint32_t gate_id = get_custom_gate_id();

    // initialize gate
    AMP u11;
    u11 = weight_lookup(cmake_angle(theta, 1));
    gates[gate_id][0] = AADD_ONE;   gates[gate_id][1] = AADD_ZERO;
    gates[gate_id][2] = AADD_ZERO;  gates[gate_id][3] = u11;

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

    k = GATEID_X;
    gates[k][0] = AADD_ZERO; gates[k][1] = AADD_ONE;
    gates[k][2] = AADD_ONE;  gates[k][3] = AADD_ZERO;

    k = GATEID_Y;
    gates[k][0] = AADD_ZERO; gates[k][1] = weight_lookup(cmake(0.0, -1.0));
    gates[k][2] = weight_lookup(cmake(0.0, 1.0));  gates[k][3] = AADD_ZERO;

    k = GATEID_Z;
    gates[k][0] = AADD_ONE;  gates[k][1] = AADD_ZERO;
    gates[k][2] = AADD_ZERO; gates[k][3] = AADD_MIN_ONE;

    k = GATEID_H;
    gates[k][0] = gates[k][1] = gates[k][2] = weight_lookup(cmake(1.0/flt_sqrt(2.0),0));
    gates[k][3] = weight_lookup(cmake(-1.0/flt_sqrt(2.0),0));

    k = GATEID_S;
    gates[k][0] = AADD_ONE;  gates[k][1] = AADD_ZERO;
    gates[k][2] = AADD_ZERO; gates[k][3] = weight_lookup(cmake(0.0, 1.0));

    k = GATEID_Sdag;
    gates[k][0] = AADD_ONE;  gates[k][1] = AADD_ZERO;
    gates[k][2] = AADD_ZERO; gates[k][3] = weight_lookup(cmake(0.0, -1.0));

    k = GATEID_T;
    gates[k][0] = AADD_ONE;  gates[k][1] = AADD_ZERO;
    gates[k][2] = AADD_ZERO; gates[k][3] = weight_lookup(cmake(1.0/flt_sqrt(2.0), 1.0/flt_sqrt(2.0)));

    k = GATEID_Tdag;
    gates[k][0] = AADD_ONE;  gates[k][1] = AADD_ZERO;
    gates[k][2] = AADD_ZERO; gates[k][3] = weight_lookup(cmake(1.0/flt_sqrt(2.0), -1.0/flt_sqrt(2.0)));

    k = GATEID_sqrtX;
    gates[k][0] = weight_lookup(cmake(0.5, 0.5)); gates[k][1] = weight_lookup(cmake(0.5,-0.5));
    gates[k][2] = weight_lookup(cmake(0.5,-0.5)); gates[k][3] = weight_lookup(cmake(0.5, 0.5));

    k = GATEID_sqrtXdag;
    gates[k][0] = weight_lookup(cmake(0.5,-0.5)); gates[k][1] = weight_lookup(cmake(0.5, 0.5));
    gates[k][2] = weight_lookup(cmake(0.5, 0.5)); gates[k][3] = weight_lookup(cmake(0.5,-0.5));

    k = GATEID_sqrtY;
    gates[k][0] = weight_lookup(cmake(0.5, 0.5)); gates[k][1] = weight_lookup(cmake(-0.5,-0.5));
    gates[k][2] = weight_lookup(cmake(0.5, 0.5)); gates[k][3] = weight_lookup(cmake(0.5, 0.5));

    k = GATEID_sqrtYdag;
    gates[k][0] = weight_lookup(cmake(0.5,-0.5)); gates[k][1] = weight_lookup(cmake(0.5,-0.5));
    gates[k][2] = weight_lookup(cmake(-0.5,0.5)); gates[k][3] = weight_lookup(cmake(0.5,-0.5));

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
    }
}