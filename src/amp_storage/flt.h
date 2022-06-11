#ifndef FLT_H
#define FLT_H

/**
 * Put all the stuff related to the specific use of floats here, to make the 
 * rest of the implementation entirely (/mostly) independent of this.
 */ 

#include <math.h>
#include <quadmath.h>

// What size float to use (double or __float128 (quad) )
#define flt_quad 0

#if flt_quad
    typedef __float128 fl_t;
#else
    typedef double fl_t;
#endif

typedef struct complex_s {
  fl_t r,i;
} complex_t;


#if flt_quad
    #define flt_abs(a) fabsq(a)
    #define flt_round(a) lroundq(a)
    #define flt_cos(a) cosq(a)
    #define flt_acos(a) acosq(a)
    #define flt_sin(a) sinq(a)
    #define flt_asin(a) asinq(a)
    #define tan(a) tanq(a)
    #define atan(a) atanq(a)
    #define atan2(a,b) atan2q(a,b)
    #define flt_sqrt(a) sqrtq(a)
#else
    #define flt_abs(a) fabs(a)
    #define flt_round(a) round(a)
    #define flt_cos(a) cos(a)
    #define flt_acos(a) acos(a)
    #define flt_sin(a) sin(a)
    #define flt_asin(a) asin(a)
    #define flt_tan(a) tan(a)
    #define flt_atan(a) atan(a)
    #define flt_atan2(a,b) atan2(a,b)
    #define flt_sqrt(a) sqrt(a)
#endif


#endif // FLT