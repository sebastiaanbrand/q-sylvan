
#define USING_MPREAL 1

#include "sylvan_qdd_complex.h"
#include "sylvan_qdd_complex_mpreal.h"



#if USING_MPREAL
    #define amp_abs(a)   mpreal_amp_abs(a)
    #define amp_neg(a)   mpreal_amp_neg(a)
    #define amp_add(a,b) mpreal_amp_add(a,b)
    #define amp_sub(a,b) mpreal_amp_sub(a,b)
    #define amp_mul(a,b) mpreal_amp_mul(a,b)
    #define amp_div(a,b) mpreal_amp_div(a,b)
    #define amp_to_prob(a) mpreal_amp_to_prob(a)
    #define prob_to_amp(a) mpreal_prob_to_amp(a)
    #define amp_exact_equal(a,b) mpreal_amp_exact_equal(a,b)
    #define amp_approx_equal(a,b) mpreal_amp_approx_equal(a,b)
    #define amp_epsilon_close(a,b,eps) mpreal_amp_epsilon_close(a,b,eps)
    #define amp_normalize_low(a,b) mpreal_amp_normalize_low(a,b)
    #define amp_normalize_largest(a,b) mpreal_amp_normalize_largest(a,b)
    #define amp_default_normalize(a,b) mpreal_amp_normalize_low(a,b) // default normalize is low
    #define init_amplitude_table(size,tol,backend) init_mpreal_amplitude_table(size,tol)
    #define free_amplitude_table() free_mpreal_amplitude_table()
#else
    
    #define amp_abs(a)   amp_abs(a)
    #define amp_neg(a)   amp_neg(a)
    #define amp_add(a,b) amp_add(a,b)
    #define amp_sub(a,b) amp_sub(a,b)
    #define amp_mul(a,b) amp_mul(a,b)
    #define amp_div(a,b) amp_div(a,b)
    #define amp_to_prob(a) amp_to_prob(a)
    #define prob_to_amp(a) prob_to_amp(a)
    #define amp_exact_equal(a,b) amp_exact_equal(a,b)
    #define amp_approx_equal(a,b) amp_approx_equal(a,b)
    #define amp_epsilon_close(a,b,eps) amp_epsilon_close(a,b,eps)
    #define amp_normalize_low(a,b) amp_normalize_low(a,b)
    #define amp_normalize_largest(a,b) amp_normalize_largest(a,b)
    #define amp_default_normalize(a,b) amp_normalize_largest(a,b) // default normalize is largest
    #define init_amplitude_table(size,tol,backend) init_amplitude_table(size,tol,backend)
    #define free_amplitude_table() free_amplitude_table()
#endif
