#ifndef SYLVAN_QMDD_GATES_H
#define SYLVAN_QMDD_GATES_H

#include <stdint.h>
#include <edge_weight_storage/flt.h>

typedef uint64_t AMP; // replace with AADD_WGT?

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
    n_predef_gates
} gate_id_t;

static const uint64_t num_static_gates  = n_predef_gates+256+256; // predef gates + phase gates
static const uint64_t num_dynamic_gates = 1000; // IDs in this rage are re-used

// 2x2 gates, k := GATEID_U 
// gates[k][0] = u00 (top left)
// gates[k][1] = u01 (top right)
// gates[k][2] = u10 (bottom left)
// gates[k][3] = u11 (bottom right)
extern uint64_t gates[n_predef_gates+256+256+1000][4]; // max 2^24 gates atm
extern gate_id_t inv_gate_ids[n_predef_gates+256+256+1000]; // inv_gate_ids[i] = ID of inv(gate_i)
// TODO: add inverese gates for Rx, Ry, Rx

void qmdd_gates_init();
// The next 255 gates are reserved for parameterized phase gates.
// The reason why these are initialized beforhand instead of on-demand is that 
// we would like a (for example) pi/16 gate to always have the same unique ID 
// throughout the entire run of the circuit.
void qmdd_phase_gates_init(int n);
static inline uint32_t GATEID_Rk(int k) { return k + n_predef_gates; };
// Another 255 parameterized phase gates, but this time with negative angles.
static inline uint32_t GATEID_Rk_dag(int k){ return k + (n_predef_gates+256); };

// Reserve 'num_dynamic_gates' gate IDs (comibned) for of the following: 
// Rx, Ry, Rz. These gate IDs are re-used when a user requires more than 
// 'num_dynamic_gates' custom gates.
/**
 * Rotation around x-axis with angle 2pi*a.
 * NOTE: These gate IDs are dynamic. The returned ID (uint32_t) is only 
 * guaranteed to correspond to the Rx(a) rotation until the next generation of 
 * a custom gate id.
 */
uint32_t GATEID_Rx(fl_t a);
/**
 * Rotation around y-axis with angle 2pi*a.
 * NOTE: These gate IDs are dynamic. The returned ID (uint32_t) is only 
 * guaranteed to correspond to the Ry(a) rotation until the next generation of 
 * a custom gate id.
 */
uint32_t GATEID_Ry(fl_t a);
/**
 * Rotation around z-axis with angle 2pi*a.
 * NOTE: These gate IDs are dynamic. The returned ID (uint32_t) is only 
 * guaranteed to correspond to the Rz(a) rotation until the next generation of 
 * a custom gate id.
 */
uint32_t GATEID_Rz(fl_t a);


#endif
