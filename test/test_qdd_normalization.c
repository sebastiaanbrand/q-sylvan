#include <stdio.h>
#include <time.h>

#include "qsylvan.h"
#include "test_assert.h"

bool VERBOSE = true;

int test_amp_normalize_sum()
{
    
    AMP a, b, c, _a, _b, a_recomputed;
    complex_t ca, cb, cc;

    // a = 1/sqrt(2), b = 1/sqrt(2) (squares already sum to 1)
    a = comp_lookup(comp_make(1.0/flt_sqrt(2.0), 0));
    b = comp_lookup(comp_make(1.0/flt_sqrt(2.0), 0));
    _a = a;
    _b = b;
    c = amp_normalize_sum(&a, &b);
    a_recomputed = amp_get_low_sum_normalized(b);

    test_assert(comp_approx_equal(comp_value(a), comp_make(1.0/flt_sqrt(2.0), 0)));
    test_assert(comp_approx_equal(comp_value(b), comp_make(1.0/flt_sqrt(2.0), 0)));
    test_assert(comp_approx_equal(comp_value(c), comp_one()));
    test_assert(a == _a);
    test_assert(b == _b);
    test_assert(c == C_ONE);
    test_assert(comp_approx_equal(comp_value(a_recomputed), comp_value(a)));
    test_assert(a_recomputed == a);

    // a = 1/sqrt(2), b = -1/sqrt(2) (squares already sum to 1)
    a = comp_lookup(comp_make(1.0/flt_sqrt(2.0), 0));
    b = comp_lookup(comp_make(-1.0/flt_sqrt(2.0), 0));
    _a = a;
    _b = b;
    c = amp_normalize_sum(&a, &b);
    a_recomputed = amp_get_low_sum_normalized(b);

    test_assert(comp_approx_equal(comp_value(a), comp_make(1.0/flt_sqrt(2.0), 0)));
    test_assert(comp_approx_equal(comp_value(b), comp_make(-1.0/flt_sqrt(2.0), 0)));
    test_assert(comp_approx_equal(comp_value(c), comp_one()));
    test_assert(a == _a);
    test_assert(b == _b);
    test_assert(c == C_ONE);
    test_assert(comp_approx_equal(comp_value(a_recomputed), comp_value(a)));
    test_assert(a_recomputed == a);

    // a = -1/sqrt(2), b = 1/sqrt(2) (squares already sum to 1, but a not in R)
    a = comp_lookup(comp_make(-1.0/flt_sqrt(2.0), 0));
    b = comp_lookup(comp_make(1.0/flt_sqrt(2.0), 0));
    _a = a;
    _b = b;
    c = amp_normalize_sum(&a, &b);
    a_recomputed = amp_get_low_sum_normalized(b);

    test_assert(comp_approx_equal(comp_value(a), comp_make(1.0/flt_sqrt(2.0), 0)));
    test_assert(comp_approx_equal(comp_value(b), comp_make(-1.0/flt_sqrt(2.0), 0)));
    test_assert(comp_approx_equal(comp_value(c), comp_minus_one()));
    test_assert(a == _b);
    test_assert(b == _a);
    test_assert(c == C_MIN_ONE);
    test_assert(comp_approx_equal(comp_value(a_recomputed), comp_value(a)));
    test_assert(a_recomputed == a);

    // a = 5, b = 3 (squares of a and b don't sum to 1)
    a = comp_lookup(comp_make(5, 0));
    b = comp_lookup(comp_make(3, 0));
    c = amp_normalize_sum(&a, &b);
    a_recomputed = amp_get_low_sum_normalized(b);
    ca = comp_value(a);
    cb = comp_value(b);
    cc = comp_value(c);

    test_assert(comp_approx_equal(comp_add(comp_sqr(ca), comp_sqr(cb)), comp_one()));
    test_assert(comp_approx_equal(comp_mul(ca, cc), comp_make(5, 0)));
    test_assert(comp_approx_equal(comp_mul(cb, cc), comp_make(3, 0)));
    test_assert(comp_approx_equal(comp_value(a_recomputed), comp_value(a)));
    test_assert(a_recomputed == a);

    // a = 5, b = 0
    a = comp_lookup(comp_make(5, 0));
    b = comp_lookup(comp_make(0, 0));
    c = amp_normalize_sum(&a, &b);
    a_recomputed = amp_get_low_sum_normalized(b);

    test_assert(comp_approx_equal(comp_value(a), comp_one()));
    test_assert(comp_approx_equal(comp_value(b), comp_zero()));
    test_assert(comp_approx_equal(comp_value(c), comp_make(5, 0)));
    test_assert(a == C_ONE);
    test_assert(b == C_ZERO);
    test_assert(comp_approx_equal(comp_value(a_recomputed), comp_value(a)));
    test_assert(a_recomputed == a);

    // a = 0, b = 0.3
    a = comp_lookup(comp_make(0, 0));
    b = comp_lookup(comp_make(0.3, 0));
    c = amp_normalize_sum(&a, &b);
    a_recomputed = amp_get_low_sum_normalized(b);

    test_assert(comp_approx_equal(comp_value(a), comp_zero()));
    test_assert(comp_approx_equal(comp_value(b), comp_one()));
    test_assert(comp_approx_equal(comp_value(c), comp_make(0.3, 0)));
    test_assert(a == C_ZERO);
    test_assert(b == C_ONE);
    test_assert(comp_approx_equal(comp_value(a_recomputed), comp_value(a)));
    test_assert(a_recomputed == a);

    // a = 2, b = 1
    a = comp_lookup(comp_make(2, 0));
    b = C_ONE;
    c = amp_normalize_sum(&a, &b);
    a_recomputed = amp_get_low_sum_normalized(b);

    test_assert(comp_approx_equal(comp_value(a), comp_make(2.0/flt_sqrt(5.0),0)));
    test_assert(comp_approx_equal(comp_value(b), comp_make(1.0/flt_sqrt(5.0),0)));
    test_assert(comp_approx_equal(comp_value(c), comp_make(flt_sqrt(5.0),0)));
    test_assert(comp_approx_equal(comp_value(a_recomputed), comp_value(a)));
    test_assert(a_recomputed == a);

    // a = 1, b = i
    a = C_ONE;
    b = comp_lookup(comp_make(0, 1));
    c = amp_normalize_sum(&a, &b);
    a_recomputed = amp_get_low_sum_normalized(b);

    test_assert(comp_approx_equal(comp_value(a), comp_make(1.0/flt_sqrt(2.0), 0)));
    test_assert(comp_approx_equal(comp_value(b), comp_make(0, 1.0/flt_sqrt(2.0))));
    test_assert(comp_approx_equal(comp_value(c), comp_make(flt_sqrt(2.0), 0)));
    test_assert(comp_approx_equal(comp_value(a_recomputed), comp_value(a)));
    test_assert(a_recomputed == a);

    // a = (1 + i)/2, b = (1 - i)/2 (magnitude 1/sqrt(2) phase difference of -i)
    a = comp_lookup(comp_make(0.5, 0.5));
    b = comp_lookup(comp_make(0.5, -0.5));
    c = amp_normalize_sum(&a, &b);
    a_recomputed = amp_get_low_sum_normalized(b);

    test_assert(comp_approx_equal(comp_value(a), comp_make(1.0/flt_sqrt(2), 0)));
    test_assert(comp_approx_equal(comp_value(b), comp_make(0, -1.0/flt_sqrt(2))));
    test_assert(comp_approx_equal(comp_value(c), comp_make(1.0/flt_sqrt(2.0), 1.0/flt_sqrt(2.0))));
    test_assert(comp_approx_equal(comp_value(a_recomputed), comp_value(a)));
    test_assert(a_recomputed == a);

    if (VERBOSE) printf("sum normalization:     ok\n");
    return 0;
}

int runtests() 
{
    init_amplitude_table(1LL<<11, -1, COMP_HASHMAP);

    // test amp normalization
    if (test_amp_normalize_sum()) return 1;

    free_amplitude_table();

    return 0;
}

int main()
{
    return runtests();
}
