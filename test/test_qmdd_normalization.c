#include <stdio.h>
#include <time.h>

#include "qsylvan.h"
#include <sylvan_edge_weights_complex.h>
#include "test_assert.h"

bool VERBOSE = true;

int test_wgt_norm_L2()
{
    
    AMP a, b, c, _a, _b, a_recomputed;
    complex_t wa, wb, wc, wr, wref, wref2;

    // a = 1/sqrt(2), b = 1/sqrt(2) (squares already sum to 1)
    a = complex_lookup(1.0/flt_sqrt(2.0), 0);
    b = complex_lookup(1.0/flt_sqrt(2.0), 0);
    _a = a;
    _b = b;
    c = wgt_norm_L2(&a, &b);
    a_recomputed = wgt_get_low_L2normed(b);

    weight_value(a, &wa);
    weight_value(b, &wb);
    weight_value(c, &wc);
    weight_value(a_recomputed, &wr);
    
    wref = cmake(1.0/flt_sqrt(2.0), 0); test_assert(weight_approx_eq(&wa, &wref));
    wref = cmake(1.0/flt_sqrt(2.0), 0); test_assert(weight_approx_eq(&wb, &wref));
    wref = cone();                      test_assert(weight_approx_eq(&wc, &wref));
    test_assert(a == _a);
    test_assert(b == _b);
    test_assert(c == EVBDD_ONE);
    test_assert(weight_approx_eq(&wr, &wa));
    test_assert(a_recomputed == a);

    // a = 1/sqrt(2), b = -1/sqrt(2) (squares already sum to 1)
    a = complex_lookup(1.0/flt_sqrt(2.0), 0);
    b = complex_lookup(-1.0/flt_sqrt(2.0), 0);
    _a = a;
    _b = b;
    c = wgt_norm_L2(&a, &b);
    a_recomputed = wgt_get_low_L2normed(b);

    weight_value(a, &wa);
    weight_value(b, &wb);
    weight_value(c, &wc);
    weight_value(a_recomputed, &wr);
    wref = cmake(1.0/flt_sqrt(2.0), 0); test_assert(weight_approx_eq(&wa, &wref));
    wref = cmake(-1.0/flt_sqrt(2.0),0); test_assert(weight_approx_eq(&wb, &wref));
    wref = cone();                      test_assert(weight_approx_eq(&wc, &wref));
    test_assert(a == _a);
    test_assert(b == _b);
    test_assert(c == EVBDD_ONE);
    test_assert(weight_approx_eq(&wr, &wa));
    test_assert(a_recomputed == a);

    // a = -1/sqrt(2), b = 1/sqrt(2) (squares already sum to 1, but a not in R)
    a = complex_lookup(-1.0/flt_sqrt(2.0), 0);
    b = complex_lookup(1.0/flt_sqrt(2.0), 0);
    _a = a;
    _b = b;
    c = wgt_norm_L2(&a, &b);
    a_recomputed = wgt_get_low_L2normed(b);

    weight_value(a, &wa);
    weight_value(b, &wb);
    weight_value(c, &wc);
    weight_value(a_recomputed, &wr);
    wref = cmake(1.0/flt_sqrt(2.0), 0); test_assert(weight_approx_eq(&wa, &wref));
    wref = cmake(-1.0/flt_sqrt(2.0),0); test_assert(weight_approx_eq(&wb, &wref));
    wref = cmone();                     test_assert(weight_approx_eq(&wc, &wref));
    test_assert(a == _b);
    test_assert(b == _a);
    test_assert(c == EVBDD_MIN_ONE);
    test_assert(weight_approx_eq(&wr, &wa));
    test_assert(a_recomputed == a);

    // a = 5, b = 3 (squares of a and b don't sum to 1)
    a = complex_lookup(5, 0);
    b = complex_lookup(3, 0);
    c = wgt_norm_L2(&a, &b);
    a_recomputed = wgt_get_low_L2normed(b);

    weight_value(a, &wa);
    weight_value(b, &wb);
    weight_value(c, &wc);
    weight_value(a_recomputed, &wr);

    wref  = cadd(cmul(wa, wa), cmul(wb, wb));
    wref2 = cone();
    test_assert(weight_approx_eq(&wref, &wref2));

    wref  = cmul(wa, wc);
    wref2 = cmake(5, 0);
    test_assert(weight_approx_eq(&wref, &wref2));

    wref  = cmul(wb, wc);
    wref2 = cmake(3, 0);
    test_assert(weight_approx_eq(&wref, &wref2));
    test_assert(weight_approx_eq(&wr, &wa));
    test_assert(a_recomputed == a);

    // a = 5, b = 0
    a = complex_lookup(5, 0);
    b = complex_lookup(0, 0);
    c = wgt_norm_L2(&a, &b);
    a_recomputed = wgt_get_low_L2normed(b);

    weight_value(a, &wa);
    weight_value(b, &wb);
    weight_value(c, &wc);
    weight_value(a_recomputed, &wr);
    wref = cone();          test_assert(weight_approx_eq(&wa, &wref));
    wref = czero();         test_assert(weight_approx_eq(&wb, &wref));
    wref = cmake(5, 0);     test_assert(weight_approx_eq(&wc, &wref));
    test_assert(a == EVBDD_ONE);
    test_assert(b == EVBDD_ZERO);
    test_assert(weight_approx_eq(&wr, &wa));
    test_assert(a_recomputed == a);

    // a = 0, b = 0.3
    a = complex_lookup(0, 0);
    b = complex_lookup(0.3, 0);
    c = wgt_norm_L2(&a, &b);
    a_recomputed = wgt_get_low_L2normed(b);

    weight_value(a, &wa);
    weight_value(b, &wb);
    weight_value(c, &wc);
    weight_value(a_recomputed, &wr);
    wref = czero();         test_assert(weight_approx_eq(&wa, &wref));
    wref = cone();          test_assert(weight_approx_eq(&wb, &wref));
    wref = cmake(0.3, 0);   test_assert(weight_approx_eq(&wc, &wref));
    test_assert(a == EVBDD_ZERO);
    test_assert(b == EVBDD_ONE);
    test_assert(weight_approx_eq(&wr, &wa));
    test_assert(a_recomputed == a);

    // a = 2, b = 1
    a = complex_lookup(2, 0);
    b = EVBDD_ONE;
    c = wgt_norm_L2(&a, &b);
    a_recomputed = wgt_get_low_L2normed(b);

    weight_value(a, &wa);
    weight_value(b, &wb);
    weight_value(c, &wc);
    weight_value(a_recomputed, &wr);
    wref = cmake(2.0/flt_sqrt(5.0),0);  test_assert(weight_approx_eq(&wa, &wref));
    wref = cmake(1.0/flt_sqrt(5.0),0);  test_assert(weight_approx_eq(&wb, &wref));
    wref = cmake(flt_sqrt(5.0),0);      test_assert(weight_approx_eq(&wc, &wref));
    test_assert(weight_approx_eq(&wr, &wa));
    test_assert(a_recomputed == a);

    // a = 1, b = i
    a = EVBDD_ONE;
    b = complex_lookup(0, 1);
    c = wgt_norm_L2(&a, &b);
    a_recomputed = wgt_get_low_L2normed(b);

    weight_value(a, &wa);
    weight_value(b, &wb);
    weight_value(c, &wc);
    weight_value(a_recomputed, &wr);
    wref = cmake(1.0/flt_sqrt(2.0),0);  test_assert(weight_approx_eq(&wa, &wref));
    wref = cmake(0, 1.0/flt_sqrt(2.0)); test_assert(weight_approx_eq(&wb, &wref));
    wref = cmake(flt_sqrt(2.0), 0);     test_assert(weight_approx_eq(&wc, &wref));
    test_assert(weight_approx_eq(&wr, &wa));
    test_assert(a_recomputed == a);

    // a = (1 + i)/2, b = (1 - i)/2 (magnitude 1/sqrt(2) phase difference of -i)
    a = complex_lookup(0.5, 0.5);
    b = complex_lookup(0.5, -0.5);
    c = wgt_norm_L2(&a, &b);
    a_recomputed = wgt_get_low_L2normed(b);

    weight_value(a, &wa);
    weight_value(b, &wb);
    weight_value(c, &wc);
    weight_value(a_recomputed, &wr);
    wref = cmake(1.0/flt_sqrt(2), 0);                       test_assert(weight_approx_eq(&wa, &wref));
    wref = cmake(0, -1.0/flt_sqrt(2));                      test_assert(weight_approx_eq(&wb, &wref));
    wref = cmake(1.0/flt_sqrt(2.0), 1.0/flt_sqrt(2.0));     test_assert(weight_approx_eq(&wc, &wref));
    test_assert(weight_approx_eq(&wr, &wa));
    test_assert(a_recomputed == a);

    if (VERBOSE) printf("L2 normalization:      ok\n");
    return 0;
}

int runtests() 
{
    sylvan_init_edge_weights(1LL<<11, 1LL<<11, -1, WGT_COMPLEX_128, COMP_HASHMAP);

    // test amp normalization
    if (test_wgt_norm_L2()) return 1;

    sylvan_edge_weights_free();

    return 0;
}

int main()
{
    return runtests();
}
