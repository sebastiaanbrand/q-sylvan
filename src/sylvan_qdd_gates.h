#ifndef SYLVAN_QDD_GATES_H
#define SYLVAN_QDD_GATES_H

#include <stdint.h>
#include "amp_storage/flt.h"

static const uint64_t num_static_gates  = 525;
static const uint64_t num_dynamic_gates = 1000; // IDs in this rage are re-used 

// 2x2 gates, k := GATEID_U 
// gates[k][0] = u00 (top left)
// gates[k][1] = u01 (top right)
// gates[k][2] = u10 (bottom left)
// gates[k][3] = u11 (bottom right)
uint64_t gates[525+1000][4]; // max 2^24 gates atm

// GATE_ID's (gates are initialized in qdd_gates_init)
// currently 24 bits available for this number (see GATE_OPID)
extern const uint32_t GATEID_I;
extern const uint32_t GATEID_X;
extern const uint32_t GATEID_Y;
extern const uint32_t GATEID_Z;
extern const uint32_t GATEID_H;
extern const uint32_t GATEID_S;
extern const uint32_t GATEID_Sdag;
extern const uint32_t GATEID_T;
extern const uint32_t GATEID_Tdag;
extern const uint32_t GATEID_sqrtX;
extern const uint32_t GATEID_sqrtXdag;
extern const uint32_t GATEID_sqrtY;
extern const uint32_t GATEID_sqrtYdag;

void qdd_gates_init();
// The next 255 gates are reserved for parameterized phase gates, which are
// initialized below. (these are GATEIDs 10 265)
// The reason why these are initialized beforhand instead of on-demand is that 
// we need a (for example) pi/16 gate always to have the same unique ID 
// throughout the entire run of the circuit.
// TODO: maybe find a better way?
void init_phase_gates(int n);
static inline uint32_t GATEID_Rk(int k) { return k + 13; };
// Another 255 parameterized phase gates, but this time with negative angles.
// (GATEIDs 266 through 521)
static inline uint32_t GATEID_Rk_dag(int k){ return k + 269; };

// Reserve 1000 gate IDs (comibned) for of the following: Rx, Ry, Rz. These 
// gate IDs are re-used when a user requires more than 1000 custom gates.
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
