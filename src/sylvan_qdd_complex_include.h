
#define USING_MPREAL 0

#include "sylvan_qdd_complex.h"
#include "sylvan_qdd_complex_mpreal.h"



#if USING_MPREAL
    #define amp_abs(a)   mpreal_amp_abs(a)
    #define amp_neg(a)   mpreal_amp_neg(a)
    #define amp_add(a,b) mpreal_amp_add(a,b)
    #define amp_sub(a,b) mpreal_amp_sub(a,b)
    #define amp_mul(a,b) mpreal_amp_mul(a,b)
    #define amp_div(a,b) mpreal_amp_div(a,b)
    #define amp_normalize_low(a,b) mpreal_amp_normalize_low(a,b)
    #define amp_normalize_largest(a,b) mpreal_amp_normalize_largest(a,b)
#else
    #define amp_abs(a)   amp_abs(a)
    #define amp_neg(a)   amp_neg(a)
    #define amp_add(a,b) amp_add(a,b)
    #define amp_sub(a,b) amp_sub(a,b)
    #define amp_mul(a,b) amp_mul(a,b)
    #define amp_div(a,b) amp_div(a,b)
    #define amp_normalize_low(a,b) amp_normalize_low(a,b)
    #define amp_normalize_largest(a,b) amp_normalize_largest(a,b)
#endif
