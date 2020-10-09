#include "shor.h"
#include "sylvan_qdd_complex.h"

static bool testing_mode = 0; // turns on/off (expensive) sanity checks

void
qdd_shor_set_testing_mode(bool on)
{
    testing_mode = on;
}

static uint64_t 
inverse_mod(uint64_t a, uint64_t N) {
    int t = 0;
    int newt = 1;
    int r = N;
    int newr = a;
    int h;
    while(newr != 0) {
        int quotient = r / newr;
        h = t;
        t = newt;
        newt = h - quotient * newt;
        h = r;
        r = newr;
        newr = h - quotient * newr;
    }
    if(r > 1)
        printf("ERROR: a is not invertible\n");
    if(t < 0)
        t = t + N;
    return t;
}

// TODO: controls first in function definition
QDD 
qdd_phi_add(QDD qdd, BDDVAR first, BDDVAR last, BDDVAR c1, BDDVAR c2, bool* a) 
{
    LACE_ME;

    QDD res = qdd;

    int k;
    int num_qubits = (last - first) + 1;
    BDDVAR qubit;
    for (int i = 0; i < num_qubits; i++) {
        qubit = first + i;
        for (int j = 0; j < (num_qubits-i); j++) {
            if (a[j] == 1) {
                k = num_qubits-j-i; // Rk = 2pi/2^k rotation
                res = qdd_cgate2(res, GATEID_Rk(k), c1, c2, qubit);
            }
        }
    }
    return res;
}

QDD 
qdd_phi_add_inv(QDD qdd, BDDVAR first, BDDVAR last, BDDVAR c1, BDDVAR c2, bool* a) 
{
    // These are all phase gates, so they'll all commute, so this is the exact
    // same function als qdd_phi_add() but with inversed angles.
    LACE_ME;

    QDD res = qdd;

    int k;
    int num_qubits = (last - first) + 1;
    BDDVAR qubit;
    for (int i = 0; i < num_qubits; i++) {
        qubit = first + i;
        for (int j = 0; j < (num_qubits-i); j++) {
            if (a[j] == 1) {
                k = num_qubits-j-i; // Rk_dag = -2pi/2^k rotation
                res = qdd_cgate2(res, GATEID_Rk_dag(k), c1, c2, qubit);
            }
        }
    }
    return res;
}

QDD
qdd_phi_add_mod(QDD qdd, BDDVAR c1, BDDVAR c2, uint64_t a, uint64_t N)
{
    LACE_ME;
    // clear cache (this function is called with different a, and cached results
    // are not parameterized on a)
    sylvan_clear_cache();
    shor_set_globals(a, N); // set bitvalues of a/N (N says the same though)
    BDDVAR nc = QDD_INVALID_VAR; // no control

    // 1.  controlled(c1,c2) phi_add(a)
    qdd = qdd_phi_add(qdd, shor_wires.targ_first, shor_wires.targ_last, c1, c2, shor_bits_a);
    // 2.  phi_add_inv(N)
    qdd = qdd_phi_add_inv(qdd, shor_wires.targ_first, shor_wires.targ_last, nc, nc, shor_bits_N);
    // 3.  QFT_inv
    qdd = qdd_circuit_QFT_inv(qdd, shor_wires.targ_first, shor_wires.targ_last);
    // 4.  CNOT (control = carry wire? = first of phi ADD, target = helper)
    qdd = qdd_gate(qdd, GATEID_H, shor_wires.helper);
    qdd = qdd_cgate(qdd, GATEID_Z, shor_wires.helper, shor_wires.targ_first);
    qdd = qdd_gate(qdd, GATEID_H, shor_wires.helper);
    // 5.  QFT
    qdd = qdd_circuit_QFT(qdd, shor_wires.targ_first, shor_wires.targ_last);
    // 6.  controlled phi_add(N) (control = helper)
    qdd = qdd_phi_add(qdd, shor_wires.targ_first, shor_wires.targ_last, shor_wires.helper, nc, shor_bits_N);
    // 7. controlled(c1, c2) phi_add_inv(a)
    qdd = qdd_phi_add_inv(qdd, shor_wires.targ_first, shor_wires.targ_last, c1, c2, shor_bits_a);
    // 8.  QFT_inv
    qdd = qdd_circuit_QFT_inv(qdd, shor_wires.targ_first, shor_wires.targ_last);
    // 9.  X on same wire as control of CNOT in 4/10
    qdd = qdd_gate(qdd, GATEID_X, shor_wires.targ_first);
    // 10. CNOT
    qdd = qdd_gate(qdd, GATEID_H, shor_wires.helper);
    qdd = qdd_cgate(qdd, GATEID_Z, shor_wires.helper, shor_wires.targ_first);
    qdd = qdd_gate(qdd, GATEID_H, shor_wires.helper);
    // 11. X on same wire as control of CNOT in 4/10
    qdd = qdd_gate(qdd, GATEID_X, shor_wires.targ_first);
    // 12. QFT
    qdd = qdd_circuit_QFT(qdd, shor_wires.targ_first, shor_wires.targ_last);
    // 13. controlled(c1,c2) phi_add(a)
    qdd = qdd_phi_add(qdd, shor_wires.targ_first, shor_wires.targ_last, c1, c2, shor_bits_a);

    return qdd;
}

QDD
qdd_phi_add_mod_inv(QDD qdd, BDDVAR c1, BDDVAR c2, uint64_t a, uint64_t N)
{
    // Inverse of function above
    LACE_ME;
    // clear cache (this function is called with different a, and cached results
    // are not parameterized on a)
    sylvan_clear_cache();
    shor_set_globals(a, N); // set bitvalues of a/N (N says the same though)
    BDDVAR nc = QDD_INVALID_VAR; // no control

    // 13. controlled(c1,c2) phi_add_inv(a)
    qdd = qdd_phi_add_inv(qdd, shor_wires.targ_first, shor_wires.targ_last, c1, c2, shor_bits_a);
    // 12. QFT^-1
    qdd = qdd_circuit_QFT_inv(qdd, shor_wires.targ_first, shor_wires.targ_last);
    // 11. X^-1 = X
    qdd = qdd_gate(qdd, GATEID_X, shor_wires.targ_first);
    // 10. CNOT^-1 = CNOT
    qdd = qdd_gate(qdd, GATEID_H, shor_wires.helper);
    qdd = qdd_cgate(qdd, GATEID_Z, shor_wires.helper, shor_wires.targ_first);
    qdd = qdd_gate(qdd, GATEID_H, shor_wires.helper);
    // 9.  X^-1 = X
    qdd = qdd_gate(qdd, GATEID_X, shor_wires.targ_first);
    // 8.  (QFT^-1)^-1 = QFT
    qdd = qdd_circuit_QFT(qdd, shor_wires.targ_first, shor_wires.targ_last);
    // 7.  controlled(c1, c2) phi_add(a)
    qdd = qdd_phi_add(qdd, shor_wires.targ_first, shor_wires.targ_last, c1, c2, shor_bits_a);
    // 6.  controlled phi_add_inv(N) (control = helper)
    qdd = qdd_phi_add_inv(qdd, shor_wires.targ_first, shor_wires.targ_last, shor_wires.helper, nc, shor_bits_N);
    // 5.  QFT^-1
    qdd = qdd_circuit_QFT_inv(qdd, shor_wires.targ_first, shor_wires.targ_last);
    // 4. CNOT^-1 = CNOT (control = carry wire? = first of phi ADD, target = helper)
    qdd = qdd_gate(qdd, GATEID_H, shor_wires.helper);
    qdd = qdd_cgate(qdd, GATEID_Z, shor_wires.helper, shor_wires.targ_first);
    qdd = qdd_gate(qdd, GATEID_H, shor_wires.helper);
    // 3. (QFT^-1)^-1 = QFT
    qdd = qdd_circuit_QFT(qdd, shor_wires.targ_first, shor_wires.targ_last);
    // 2.  phi_add(N)
    qdd = qdd_phi_add(qdd, shor_wires.targ_first, shor_wires.targ_last, nc, nc, shor_bits_N);
    // 1.  controlled(c1,c2) phi_add_inv(a)
    qdd = qdd_phi_add_inv(qdd, shor_wires.targ_first, shor_wires.targ_last, c1, c2, shor_bits_a);

    return qdd;
}


QDD
qdd_cmult(QDD qdd, uint64_t a, uint64_t N)
{
    // this implements the _controlled_ cmult operation
    // 1. QFT on bottom register
    qdd = qdd_circuit_QFT(qdd, shor_wires.targ_first, shor_wires.targ_last);

    // 2. loop over k = {0, n-1}
    uint64_t t = a;
    //BDDVAR cs[] = {shor_wires.top, QDD_INVALID_VAR, QDD_INVALID_VAR};
    for (BDDVAR c2 = shor_wires.ctrl_last; c2 >= shor_wires.ctrl_first; c2--) {
        // 2a. double controlled phi_add_mod(a* 2^k)
        qdd = qdd_phi_add_mod(qdd, shor_wires.top, c2, t, N);
        t = (2*t) % N;
    }

    // 3. QFT^-1 on bottom register
    qdd = qdd_circuit_QFT_inv(qdd, shor_wires.targ_first, shor_wires.targ_last);

    return qdd;
}

QDD
qdd_cmult_inv(QDD qdd, uint64_t a, uint64_t N)
{
    // not quite inverse of above
    // 1. QFT on bottom register
    qdd = qdd_circuit_QFT(qdd, shor_wires.targ_first, shor_wires.targ_last);
    
    // 2. same loop over k but with phi_add_mod_inv
    uint64_t t = a;
    for (BDDVAR c2 = shor_wires.ctrl_last; c2 >= shor_wires.ctrl_first; c2--) {
        // 2a. double controlled phi_add_mod_inv(a* 2^k)
        qdd = qdd_phi_add_mod_inv(qdd, shor_wires.top, c2, t, N);
        t = (2*t) % N;
    }

    // 3. QFT^-1 on bottom register
    qdd = qdd_circuit_QFT(qdd, shor_wires.targ_first, shor_wires.targ_last);

    return qdd;
}

QDD
qdd_shor_ua(QDD qdd,  uint64_t a, uint64_t N)
{
    LACE_ME;

    // 1. controlled Cmult(a)
    qdd = qdd_cmult(qdd, a, N);

    // 2. controlled swap top/bottom registers
    BDDVAR cs[] = {shor_wires.top, QDD_INVALID_VAR, QDD_INVALID_VAR};
    for (uint32_t i = shor_wires.ctrl_first; i <= shor_wires.ctrl_last; i++) {
        qdd = qdd_ccircuit(qdd, CIRCID_swap, cs, i, shor_wires.targ_first+i);
    }

    // 3. controlled Cmult_inv(a^-1)
    uint64_t a_inv = inverse_mod(a, N);
    qdd = qdd_cmult_inv(qdd, a_inv, N);

    return qdd;
}

uint64_t
shor_period_finding(uint64_t a, uint64_t N)
{
    // Circuit (quantum period finding of f(x) = a^x mod N)
    // create QDD |0>|0..001>|0>|0..00>
    uint32_t num_qubits = 2*shor_n + 3;
    bool x[num_qubits];
    for (BDDVAR k = 0; k < num_qubits; k++) x[k] = 0;
    x[shor_wires.ctrl_last] = 1; // set the input reg. to |0...001> = |1>

    QDD qdd = qdd_create_basis_state(num_qubits, x);

    uint64_t as[2*shor_n];
    as[2*shor_n-1] = a;
    uint64_t new_a = a;
    for(int i = 2*shor_n-2; i >= 0; i--) {
        new_a = new_a * new_a;
        new_a = new_a % N;
        as[i] = new_a;
    }

    LACE_ME;

    int m_outcomes[2*shor_n];
    int m_outcome;
    double m_prob;

    for (uint32_t i = 0; i < 2*shor_n; i++) {
        if (testing_mode) assert(qdd_is_unitvector(qdd, num_qubits));

        // H on top wire
        qdd = qdd_gate(qdd, GATEID_H, shor_wires.top);

        // controlled Ua^...
        qdd = qdd_shor_ua(qdd, as[i], N);

        // phase gates based on previous measurement
        int k = 2; // First gate needs to be R^dag(2) = S^dag
        for (int j = i-1; j >= 0; j--) {
            if (m_outcomes[j] == 1)
                qdd = qdd_gate(qdd, GATEID_Rk_dag(k), shor_wires.top);
            k += 1; // R(k) is a (2*pi / 2^k) rotation
        }

        // H on top wire
        qdd = qdd_gate(qdd, GATEID_H, shor_wires.top);

        // measure q0
        //sylvan_clear_cache();
        qdd = qdd_measure_qubit(qdd, shor_wires.top, num_qubits, &m_outcome, &m_prob);
        m_outcomes[i] = m_outcome;

        // make sure q0 is in the |0> state
        if (m_outcome == 1) qdd = qdd_gate(qdd, GATEID_X, shor_wires.top);
    }

    // turn measurement outcomes into an integer
    uint64_t res = 0;
    for (uint32_t i = 0; i < 2*shor_n; i++) {
        int index = 2*shor_n-1-i;
        res = (res << 1) + m_outcomes[index];
    }
    return res;
}

void
shor_set_globals(uint64_t a, uint64_t N) 
{
    shor_n = ceil(log2(N)); // number of bits for N (not the number of qubits!)  
    uint64_t p2 = 1;
    for (uint32_t i = 0; i < shor_n; i++) { // LSB in bits[0], MSB in bits[63]
        shor_bits_a[i] = a & p2;
        shor_bits_N[i] = N & p2;
        p2 = p2 << 1;
    }
    // Set wire numbers
    shor_wires.top        = 0;
    shor_wires.ctrl_first = 1;
    shor_wires.ctrl_last  = shor_n;
    shor_wires.helper     = shor_n + 1; // easier to have this in the middle
    shor_wires.targ_first = shor_n + 2;
    shor_wires.targ_last  = 2*shor_n + 2;
}

uint64_t 
my_gcd (uint64_t a, uint64_t b) // clash with gcd in sylvan_mtbdd.c ...
{
  uint64_t c;
  while ( a != 0 ) { c = a; a = b%a;  b = c; }
  return b;
}

uint64_t modpow(uint64_t base, uint64_t exp, uint64_t modulus) {
  base %= modulus;
  uint64_t result = 1;
  while (exp > 0) {
    if (exp & 1) result = (result * base) % modulus;
    base = (base * base) % modulus;
    exp >>= 1;
  }
  return result;
}

uint64_t
shor_post_process(uint64_t N, uint64_t a, uint64_t b, uint64_t denom, bool verbose)
{
    // For b the following is true:
    // b/denom = x/r, where denom = 2^num bits, and r = the period we want.
    // This function tries to find that r.
    // Implementation from [zulehner2018advanced]
    if (b == 0) {
        if (verbose)
            printf("Factorization failed (measured 0)\n");
        return 0;
    }

    int cf_max_size = 100;
    int cf_entries = 0;
    uint64_t cf[cf_max_size];
    uint64_t old_b = b;
	uint64_t old_denom = denom;
	while(b != 0) {
        if (cf_entries >= cf_max_size) {
            printf("please hardcode cf_max_size to something bigger\n"); // I'm sorry
            exit(1);
        }
        cf[cf_entries] = (denom/b);
        cf_entries++;
		uint64_t tmp = denom % b;
		denom = b;
		b = tmp;
	}

    if (verbose) {
        printf("Continued fraction expansion of %ld/%ld = ", b, denom);
        for(int i = 0; i < cf_entries; i++) printf("%ld ", cf[i]);
        printf("\n");
    }

	for(int i=0; i < cf_entries; i++) {
		//determine candidate
		uint64_t denominator = cf[i];
		uint64_t numerator = 1;

		for(int j=i-1; j >= 0; j--) {
			uint64_t tmp = numerator + cf[j]*denominator;
			numerator = denominator;
			denominator = tmp;
		}
        if (verbose)
            printf(" Candidate %ld/%ld: ", numerator, denominator);
		if(denominator > N) {
            if (verbose) {
                printf(" denominator too large (greater than %ld)\n", N);
                printf("Factorization failed\n");
            }
            return 0;
		} else {
			double delta = (double)old_b / (double)old_denom - (double)numerator / (double) denominator;
			if(fabs(delta) < 1.0/(2.0*old_denom)) {
				if(modpow(a, denominator, N) == 1) {
                    if (verbose)
                        printf("found period = %ld\n", denominator);
					if(denominator & 1) {
                        if (verbose)
                            printf("Factorization failed (period is odd)\n");
                        return 0;
					} else {
						uint64_t f1, f2;
						f1 = modpow(a, denominator>>1, N);
						f2 = (f1+1)%N;
						f1 = (f1 == 0) ? N-1 : f1-1;
						f1 = my_gcd(f1, N);
						f2 = my_gcd(f2, N);
                        if (f1 == 1 || f1 == N) {
                            if (verbose)
                                printf("Factorization found trivial factor\n");
                            return 0;
                        }
                        if (verbose) {
                            printf("Factorization succeeded! Non-trivial factors are:\n");
                            printf(" -- gcd(%ld^(%ld/2)-1,%ld)=%ld\n", N, denominator, N, f1);
                            printf(" -- gcd(%ld^(%ld/2)+1,%ld)=%ld\n", N, denominator, N, f2);
                        }
                        return f1;
					}
					break;
				} else {
                    if (verbose)
                        printf("failed\n");
				}
			} else {
                if (verbose)
                    printf("delta is too big (%lf)\n", delta);
			}
		}
	}
    return 0;
}

uint64_t
shor_generate_a(uint64_t N)
{
    uint64_t a;
    do {
        a = rand() % N;
    } while (my_gcd(a, N) != 1 || a == 1);
    return a;
}

uint64_t
run_shor(uint64_t N, uint64_t a, bool verbose)
{
    // The classical part
    if (a == 0) a = shor_generate_a(N);

    shor_set_globals(a, N);
    
    if (verbose) {
        printf("input N        = %ld [", N);
        for (uint32_t i=0; i<shor_n; i++) printf("%d", shor_bits_N[i]);
        printf("]\n");
        printf("n (bits for N) = %d\n",  shor_n);
        printf("random a       = %ld [", a);
        for (uint32_t i=0; i<shor_n; i++) printf("%d", shor_bits_a[i]);
        printf("]\n\n");

        printf("wires:\n");
        printf("top:        %d\n", shor_wires.top);
        printf("ctrl_first: %d\n", shor_wires.ctrl_first);
        printf("ctrl_last:  %d\n", shor_wires.ctrl_last);
        printf("helper:     %d\n", shor_wires.helper);
        printf("targ_first: %d\n", shor_wires.targ_first);
        printf("targ_last:  %d\n\n", shor_wires.targ_last);
    }

    uint64_t b = shor_period_finding(a, N);
    uint64_t denom = 1 << (2*shor_n);

    return shor_post_process(N, a, b, denom, verbose);
}