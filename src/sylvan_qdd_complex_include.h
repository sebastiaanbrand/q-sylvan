
#define USING_MPREAL 0

#include "sylvan_qdd_complex.h"
#include "sylvan_qdd_complex_mpreal.h"



#if USING_MPREAL
    #define amp_abs(a)   mpreal_amp_abs(a)
    #define amp_neg(a)   mpreal_amp_neg(a)
    #define amp_add(a,b) mpreal_amp_add(a,b)
    #define amp_sub(a,b) mpreal_amp_sub(a,b)
    #define amp_mul(a,b) mpreal_amp_mul(a,b)
    #define amp_mul_down(a,b) mpreal_amp_mul(a,b)
    #define amp_div(a,b) mpreal_amp_div(a,b)
    #define amp_to_prob(a) mpreal_amp_to_prob(a)
    #define prob_to_amp(a) mpreal_prob_to_amp(a)
    #define amp_exact_equal(a,b) mpreal_amp_exact_equal(a,b)
    #define amp_approx_equal(a,b) mpreal_amp_approx_equal(a,b)
    #define amp_epsilon_close(a,b,eps) mpreal_amp_epsilon_close(a,b,eps)
    #define amp_normalize_low(a,b) mpreal_amp_normalize_low(a,b)
    #define amp_normalize_largest(a,b) mpreal_amp_normalize_largest(a,b)
    #define amp_normalize_low_ptr &mpreal_amp_normalize_low
    #define amp_normalize_largest_ptr &mpreal_amp_normalize_largest
    #define init_amplitude_table(size,tol,backend) init_mpreal_amplitude_table(size,tol)
    #define free_amplitude_table() free_mpreal_amplitude_table()
    #define amp_store_get_tolerance() mpreal_amp_store_get_tolerance()

    // for gc
    #define init_new_empty_table() init_new_empty_mpreal_table()
    #define delete_old_table() delete_old_mpreal_table()
    #define move_from_old_to_new(a) move_from_old_to_new_mpreal(a)
    #define get_table_entries_estimate() get_mpreal_table_num_entries() // TODO: replace w/ est
    #define count_amplitude_table_enries() get_mpreal_table_num_entries()
    #define get_table_size() get_mpreal_table_size()
#else
    
    #define amp_abs(a)   amp_abs(a)
    #define amp_neg(a)   amp_neg(a)
    #define amp_add(a,b) amp_add(a,b)
    #define amp_sub(a,b) amp_sub(a,b)
    #define amp_mul(a,b) amp_mul(a,b)
    #define amp_mul_down(a,b) amp_mul_down(a,b)
    #define amp_div(a,b) amp_div(a,b)
    #define amp_to_prob(a) amp_to_prob(a)
    #define prob_to_amp(a) prob_to_amp(a)
    #define amp_exact_equal(a,b) amp_exact_equal(a,b)
    #define amp_approx_equal(a,b) amp_approx_equal(a,b)
    #define amp_epsilon_close(a,b,eps) amp_epsilon_close(a,b,eps)
    #define amp_normalize_low(a,b) amp_normalize_low(a,b)
    #define amp_normalize_largest(a,b) amp_normalize_largest(a,b)
    #define amp_normalize_low_ptr &amp_normalize_low
    #define amp_normalize_largest_ptr &amp_normalize_largest
    #define init_amplitude_table(size,tol,backend) init_amplitude_table(size,tol,backend)
    #define free_amplitude_table() free_amplitude_table()
    #define amp_store_get_tolerance() amp_store_get_tolerance()

    // for gc
    #define init_new_empty_table() init_new_empty_table()
    #define delete_old_table() delete_old_table()
    #define move_from_old_to_new(a) move_from_old_to_new(a)
    #define get_table_entries_estimate() get_table_entries_estimate()
    #define count_amplitude_table_enries() count_amplitude_table_enries()
    #define get_table_size() get_table_size()
#endif
