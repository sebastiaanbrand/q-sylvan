#ifndef FLT_H
#define FLT_H

/**
 * Put all the stuff related to the specific use of floats here, to make the 
 * rest of the implementation entirely (/mostly) independent of this.
 */ 

#include <math.h>

// What size float to use (double or __float128 (quad) )
#define flt_quad 0 // (setting to 1 might not work on MACOS)

#if flt_quad
    #include <quadmath.h>
    typedef __float128 fl_t;
#else
    typedef double fl_t;
#endif


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


// TODO: define complex_t (and some basic operations) in separate file
/*****************************< complex_t >************************************/
typedef struct complex_s {
  fl_t r,i;
} complex_t;

inline complex_t cmake(fl_t r, fl_t i)
{
    complex_t res;
    res.r = r;
    res.i = i;
    return res;
}
inline complex_t cmake_angle(fl_t theta, fl_t mag)
{
    complex_t c;
    c.r = flt_cos(theta) * mag;
    c.i = flt_sin(theta) * mag;
    return c;
}
inline complex_t czero() { return cmake(0.0, 0.0); }
inline complex_t cone() { return cmake(1.0, 0.0); }
inline complex_t cmone() { return cmake(-1.0, 0.0); }

inline complex_t cadd(complex_t a, complex_t b)
{
    complex_t res;
    res.r = a.r + b.r;
    res.i = a.i + b.i;
    return res;
}
inline complex_t csub(complex_t a, complex_t b)
{
    complex_t res;
    res.r = a.r - b.r;
    res.i = a.i - b.i;
    return res;
}
inline complex_t cmul(complex_t a, complex_t b)
{
    complex_t res;
    res.r = a.r * b.r - a.i * b.i;
    res.i = a.r * b.i + a.i * b.r;
    return res;
}
inline complex_t cdiv(complex_t a, complex_t b)
{
    complex_t res;
    fl_t denom;
    if (b.i == 0.0) {
        res.r = a.r / b.r;
        res.i = a.i / b.r;
    } else {
        denom = b.r * b.r + b.i * b.i;
        res.r = (a.r * b.r + a.i * b.i) / denom;
        res.i = (a.i * b.r - a.r * b.i) / denom;
    }
    return res;
}

/****************************< /complex_t >************************************/


#endif // FLT