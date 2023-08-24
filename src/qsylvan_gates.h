#ifndef SYLVAN_QMDD_GATES_H
#define SYLVAN_QMDD_GATES_H

#include <stdint.h>
#include <edge_weight_storage/flt.h>

typedef uint64_t AMP; // TODO: replace with AADD_WGT?

// GATE_ID's (gates are initialized in qmdd_gates_init)
// currently 24 bits available for this number (see GATE_OPID)
typedef enum predef_gates {
    GATEID_I,
    GATEID_X,
    GATEID_Y,
    GATEID_Z,
    GATEID_H,
    GATEID_S,
    GATEID_Sdag,
    GATEID_T,
    GATEID_Tdag,
    GATEID_sqrtX,
    GATEID_sqrtXdag,
    GATEID_sqrtY,
    GATEID_sqrtYdag,
    GATEID_dynamic, // used for Rz,Ry,Rz,Phase,U gates
    n_predef_gates
} gate_id_t;

static const uint64_t num_static_gates  = n_predef_gates+256+256; // predef gates + phase gates

// 2x2 gates, k := GATEID_U 
// gates[k][0] = u00 (top left)
// gates[k][1] = u01 (top right)
// gates[k][2] = u10 (bottom left)
// gates[k][3] = u11 (bottom right)
extern uint64_t gates[n_predef_gates+256+256][4]; // max 2^24 gates atm

void qmdd_gates_init();
// The next 255 gates are reserved for parameterized phase gates.
// The reason why these are initialized beforhand instead of on-demand is that 
// we would like a (for example) pi/16 gate to always have the same unique ID 
// throughout the entire run of the circuit.
// TODO: maybe just use GATEID_dynamic for this
void qmdd_phase_gates_init(int n);
static inline uint32_t GATEID_Rk(int k) { return k + n_predef_gates; };
// Another 255 parameterized phase gates, but this time with negative angles.
static inline uint32_t GATEID_Rk_dag(int k){ return k + (n_predef_gates+256); };

// Reserve 'num_dynamic_gates' gate IDs (comibned) for of the following: 
// Rx, Ry, Rz. These gate IDs are re-used when a user requires more than 
// 'num_dynamic_gates' custom gates.
/**
 * Rotation around x-axis with angle theta.
 * NOTE: This GATEID is re-used for parametrized gates. The returned ID only
 * corresponds to the Rx(theta) gate until the next parametrized gate is created.
 */
uint32_t GATEID_Rx(fl_t theta);
/**
 * Rotation around y-axis with angle theta.
 * NOTE: This GATEID is re-used for parametrized gates. The returned ID only
 * corresponds to the Ry(theta) gate until the next parametrized gate is created.
 * a custom gate id.
 */
uint32_t GATEID_Ry(fl_t theta);
/**
 * Rotation around z-axis with angle theta.
 * NOTE: This GATEID is re-used for parametrized gates. The returned ID only
 * corresponds to the Rz(theta) gate until the next parametrized gate is created.
 */
uint32_t GATEID_Rz(fl_t theta);
/**
 * Rotation around z-axis with angle theta (but different global phase than Rz)
 * NOTE: This GATEID is re-used for parametrized gates. The returned ID only
 * corresponds to the P(theta) gate until the next parametrized gate is created.
 */
uint32_t GATEID_Phase(fl_t theta);
/**
 * Generic single-qubit rotation gate with 3 Euler angles.
 * NOTE: This GATEID is re-used for parametrized gates. The returned ID only
 * corresponds to the U(theta,phi,lambda) gate until the next parametrized gate 
 * is created.
 */
uint32_t GATEID_U(fl_t theta, fl_t phi, fl_t lambda);



#endif
