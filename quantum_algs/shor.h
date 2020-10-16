#include <stdbool.h>

#include "sylvan.h"

/**
 * The flollowing functions are a breakdown of the components needed for Shor
 * as in Beauregard, "Circuit for Shor's algorithm using 2n+ 3 qubits." (2002).
 */

/**
 * Beauregard (2002) Fig. 3. 
 * Addition in Fourier space. Important here to note is the endianess:
 * - Input/output statevector |x> = |q0, q1, q2>, then q0 contains the MSB.
 * - Classical bit-array a has a[0] as LSB. This is done to easier allow for 
 *   leading zeros (in a[] now trailing zeros) when enoding the number in more 
 *   bits than it needs.
 * - Carries happen from q(k) -> q(k-1), i.e. towards the MSB, so if we write
 *   the state as |q0, q1, q2> carries go to the left (as normal).
 * 
 * @param qdd A QDD encoding a state |phi(x)> = QFT|x> with |x> a z-basis state.
 * @param a A big-endian (MSB first) encoding of some integer.
 * 
 * @return A QDD encoding |phi(x + a)>, with (x+a)
 */
QDD qdd_phi_add(QDD qdd, BDDVAR first, BDDVAR last, BDDVAR c1, BDDVAR c2, bool* a); // Fig. 3
QDD qdd_phi_add_inv(QDD qdd, BDDVAR first, BDDVAR last,  BDDVAR c1, BDDVAR c2, bool* a);

/* Beauregard (2002) Fig. 5 */
QDD qdd_phi_add_mod(QDD qdd, BDDVAR c1, BDDVAR c2, uint64_t a, uint64_t N);
QDD qdd_phi_add_mod_inv(QDD qdd, BDDVAR c1, BDDVAR c2, uint64_t a, uint64_t N);

/* Beauregard (2002) Fig. 6 */
QDD qdd_cmult(QDD qdd, uint64_t a, uint64_t N);
QDD qdd_cmult_inv(QDD qdd, uint64_t a, uint64_t N);

/* Beauregard (2002) Fig. 7 */
QDD qdd_shor_ua(QDD qdd, uint64_t a, uint64_t N);

/* Beauregard (2002) Fig. 8 */
uint64_t shor_period_finding(uint64_t a, uint64_t N);


uint64_t shor_generate_a(uint64_t N);
void shor_set_globals(uint64_t a, uint64_t N);

/**
 * N is the number to factor, and 'a' is the value to use in a^x mod N. 
 * If 'a' is set to 0 a random 'a' is chosen.
 */
uint64_t run_shor(uint64_t N, uint64_t a, bool verbose);

void qdd_shor_set_testing_mode(bool on);
