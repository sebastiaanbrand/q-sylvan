/*
 * Copyright 2023 System Verification Lab, LIACS, Leiden University
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <sylvan_int.h>
#include <sylvan_mpc.h>

#include <math.h>
#include <string.h>

/**
 * helper function for hash (skip)
 */
#ifndef rotl64
static inline uint64_t
rotl64(uint64_t x, int8_t r)
{
    return ((x<<r) | (x>>(64-r)));
}
#endif

/**
 * Custom leaf operations of type dependent functions for 
 * complex numbers represented by multi precision complex types.
 * 
 *      hashing, 
 *      compare, 
 *      create, 
 *      destroy, 
 *      and print.
 * 
 */
static uint64_t
mpc_hash(const uint64_t val, const uint64_t seed)
{
    //
    // Calculate the hash based in a mpc_t complex number
    //
    // The precision of the real and imaginary parts are
    // limited, and corresponds to the long double type 
    // which has a size in bytes of 16. 
    //

    // mpc_ptr x = (mpc_ptr)(size_t)val;

    // const uint64_t prime = 1099511628211;
    // uint64_t hash = seed;
    // mp_limb_t *limbs; // Contains unsigned longs

    //
    // Calculate the hash based on the real part of the complex number val
    //
    mpfr_t real;
    mpfr_init2(real, MPC_PRECISION);
    mpc_real(real, (mpc_ptr)val, MPC_ROUNDING);

    // Convert the real part from mpc to long double type (16 x 8 bits = 128 bits) 
    long double real_limited = mpfr_get_ld(real, MPC_ROUNDING);

    // Convert the limited real part in long double (16 bytes in ARM64) to an array of bytes
    int nr_bytes_complex_parts = 16;
    unsigned char bytes[nr_bytes_complex_parts];
    assert(nr_bytes_complex_parts == sizeof(long double)); // Check if the OS supports long double, if not exit this program

    memcpy(bytes, &real_limited, nr_bytes_complex_parts);

    const uint64_t prime = 1099511628211;
    uint64_t hash = seed;
    int nr_bytes_hash = 8;
    assert(nr_bytes_hash == sizeof(uint64_t)); // Check if the OS supports long double, if not exit this program

    for(int i=0; i<nr_bytes_complex_parts; i++) {
        hash = hash ^ bytes[i];
        hash = rotl64(hash, 47);
        hash = hash * prime;
    }

    //
    // Calculate the hash further based on the imaginary part of the complex number val
    //
    mpfr_t imag;
    mpfr_init2(imag, MPC_PRECISION);
    mpc_imag(imag, (mpc_ptr)val, MPC_ROUNDING);

    // Convert the imaginary part from mpc to long double type (16 x 8 bits = 128 bits) 
    long double imag_limited = mpfr_get_ld(imag, MPC_ROUNDING);

    // Convert the limited imaginary part in long double (16 bytes in ARM64) to an array of bytes
    memcpy(bytes, &imag_limited, nr_bytes_complex_parts);

    for(int i=0; i<nr_bytes_complex_parts; i++) {
        hash = hash ^ bytes[i];
        hash = rotl64(hash, 31);
        hash = hash * prime;
    }

//
    if(false) {
        // Example conversion from long double to bytes
        long double long_double = real_limited;

        // Create an array of bytes
        unsigned char bytes[sizeof(long double)];

        // Copy the bytes of the long double to the array
        memcpy(bytes, &long_double, sizeof(long double));

        // Display the bytes in hexadecimal format
        printf("real limited to long double value: %.17Lf\n", long_double);
        printf("bytes in hexadecimal: ");
    
        for (size_t i = 0; i < sizeof(long double); ++i) {
            printf("%02x ", bytes[i]);
        }
        printf("\n");

        // Example conversion from long double to bytes
        long_double = imag_limited;

        // Copy the bytes of the long double to the array
        memcpy(bytes, &long_double, sizeof(long double));

        // Display the bytes in hexadecimal format
        printf("imag limited to long double value: %.17Lf\n", long_double);
        printf("bytes in hexadecimal: ");
    
        for (size_t i = 0; i < sizeof(long double); ++i) {
            printf("%02x ", bytes[i]);
        }
        printf("\n");
    }
//

    mpfr_clear(real);
    mpfr_clear(imag);

    return hash ^ (hash >> 32);
}

static int
mpc_equals(const uint64_t left, const uint64_t right)
{
    //
    // This function is called by the unique table when comparing a new 
    // leaf with an existing leaf.
    //

    mpc_ptr x = (mpc_ptr)(size_t)left;
    mpc_ptr y = (mpc_ptr)(size_t)right;

printf("mpc_equals analysis\n\n");

    mpc_out_str(stdout, MPC_BASE_OF_FLOAT, 3, x, MPC_ROUNDING);
    putchar('\n');

    printf("prec = %ld\n", x->re->_mpfr_prec);
    printf("sign = %d\n",  x->re->_mpfr_sign);
    printf("exp  = %ld\n", x->re->_mpfr_exp);
    printf("d    = %p\n",  x->re->_mpfr_d);

    printf("prec = %ld\n", y->re->_mpfr_prec);
    printf("sign = %d\n",  y->re->_mpfr_sign);
    printf("exp  = %ld\n", y->re->_mpfr_exp);
    printf("d    = %p\n",  y->re->_mpfr_d);

/*
        // Example conversion from long double to bytes
        long double long_double = real_limited;

        // Create an array of bytes
        unsigned char bytes[sizeof(long double)];

        // Copy the bytes of the long double to the array
        memcpy(bytes, &long_double, sizeof(long double));

        // Display the bytes in hexadecimal format
        printf("real limited to long double value: %.17Lf\n", long_double);
        printf("bytes in hexadecimal: ");
    
        for (size_t i = 0; i < sizeof(long double); ++i) {
            printf("%02x ", bytes[i]);
        }
        printf("\n");
*/

    mpc_out_str(stdout, MPC_BASE_OF_FLOAT, 3, y, MPC_ROUNDING);
    putchar('\n');

    mpfr_t real_x;
    mpfr_init2(real_x, MPC_PRECISION);
    int res = mpc_real(real_x, (mpc_ptr)x, MPC_ROUNDING);
printf("res = %d\n", res);
    long double real_x_ld = mpfr_get_ld(real_x, MPC_ROUNDING);

printf("real_x = %Lf\n", real_x_ld);

    mpfr_t imag_x;
    mpfr_init2(imag_x, MPC_PRECISION);
    res = mpc_imag(imag_x, (mpc_ptr)x, MPC_ROUNDING);
printf("res = %d\n", res);
    long double imag_x_ld = mpfr_get_ld(imag_x, MPC_ROUNDING);

printf("imag_x = %Lf\n", imag_x_ld);

    mpfr_t real_y;
    mpfr_init2(real_y, MPC_PRECISION);
    res = mpc_real(real_y, (mpc_ptr)y, MPC_ROUNDING);
printf("res = %d\n", res);
    long double real_y_ld = mpfr_get_ld(real_y, MPC_ROUNDING);

printf("real_y = %Lf\n", real_y_ld);

    mpfr_t imag_y;
    mpfr_init2(imag_y, MPC_PRECISION);
printf("--\n");
    res = mpc_imag(imag_y, (mpc_ptr)y, MPC_ROUNDING);
printf("res = %d\n", res);
    long double imag_y_ld = mpfr_get_ld(imag_y, MPC_ROUNDING);

printf("imag_y = %Lf\n", imag_y_ld);

    int result = 0;
    result = ((real_x_ld == real_y_ld) && (imag_x_ld == imag_y_ld)) ? 1 : 0;

    mpfr_clear(real_x);
    mpfr_clear(imag_x);

    mpfr_clear(real_y);
    mpfr_clear(imag_y);

    printf("sylvan_mpc.c mpc_equals() own method = %d \n", result);

    printf("sylvan_mpc.c mpc_equals() mpc_cmp() = %d \n", mpc_cmp(x, y) ? 1 : 0);

    return result;
}

static void
mpc_create(uint64_t *val)
{
    //
    // This function is called by the unique table when a leaf does not yet exist.
    // We make a copy, which will be stored in the hash table.
    //

    mpc_ptr x = (mpc_ptr)malloc(sizeof(mpc_t));
    mpc_init2(x, MPC_PRECISION);
    mpc_set(x, *(mpc_ptr*)val, MPC_ROUNDING);
    *(mpc_ptr*)val = (mpc_ptr)x;

    return;
}

static void
mpc_destroy(uint64_t val)
{
    //
    // This function is called by the unique table
    // when a leaf is removed during garbage collection. 
    //

    mpc_clear((mpc_ptr)val);
    free((void*)val);

    return;
}

/*
static char*
mpc_to_str(int comp, uint64_t val, char *buf, size_t buflen)
{
    //
    // Generate a string of real and imaginary part enclosed 
    // in a pair of parentheses of the complex number val.
    //

    mpc_ptr op = (mpc_ptr)val;

    buf = mpc_get_str(MPC_BASE_OF_FLOAT, buflen, op, (mpc_rnd_t)MPC_ROUNDING);

    (void)comp;  // TODO: make pretty!

    return NULL;
}
*/

/*
static int
mpc_write_binary(FILE* out, uint64_t val)
{
    //mpc_ptr op = (mpc_ptr)val;

    mpc_t x;
    mpc_init2(x, MPC_PRECISION);
    mpc_set(x, *(mpc_ptr*)val, MPC_ROUNDING);

    mpc_out_str(out, MPC_BASE_OF_FLOAT, MPC_MAXLENGTH_FILESTRING, x, MPC_ROUNDING);

    return 0;
}

static int
mpc_read_binary(FILE* in, uint64_t *val)
{
    mpc_t x;
    mpc_init2(x, MPC_PRECISION);

    mpc_inp_str(x, in, (long unsigned int*)MPC_MAXLENGTH_FILESTRING, MPC_BASE_OF_FLOAT, MPC_ROUNDING);

    mpc_get(x, *(mpc_ptr*)val, MPC_ROUNDING);

    return 0;
}
*/

/**
 * Initialize mpc custom leaves
 */
uint32_t
mpc_init()
{
    // Register custom leaf type callback functions
    uint32_t mpc_type = sylvan_mt_create_type();

    // Basic functions
    sylvan_mt_set_hash(mpc_type, mpc_hash);
    sylvan_mt_set_equals(mpc_type, mpc_equals);
    sylvan_mt_set_create(mpc_type, mpc_create);
    sylvan_mt_set_destroy(mpc_type, mpc_destroy);

    // Printing / storage functions
    sylvan_mt_set_to_str(mpc_type, NULL);         // mpc_to_str);         Switched out for now.
    sylvan_mt_set_write_binary(mpc_type, NULL);   // mpc_write_binary);
    sylvan_mt_set_read_binary(mpc_type, NULL);    // mpc_read_binary);

    return mpc_type;
}

/**
 * Create mpc leaf
 */
MTBDD
mtbdd_mpc(mpc_t val)
{
    uint32_t mpc_type = MPC_TYPE;
    uint32_t g_mpc_type = MPC_TYPE;

printf("sylvan_mpc.c mtbdd_mpc(val) g_mpc_type = %d\n", g_mpc_type);

    return mtbdd_makeleaf(mpc_type, (size_t)val);
}

/**
 * Assign a complex number
*/
void 
mpc_assign(mpc_ptr complexnumber, double real, double imag)
{
    mpc_init2(complexnumber, MPC_PRECISION);
    mpc_set_d_d(complexnumber, real, imag, MPC_ROUNDING);
    return;
}


/**
 * Compare mpc leafs
*/
int
mpc_compare(const uint64_t left, const uint64_t right)
{
    //
    // This function is called by the unique table when comparing a new 
    // leaf with an existing leaf.
    //

    mpc_ptr x = (mpc_ptr)(size_t)left;
    mpc_ptr y = (mpc_ptr)(size_t)right;

printf("sylvan_mpc.c mpc_compare()\n\n");

    //mpc_out_str(stdout, MPC_BASE_OF_FLOAT, 3, x, MPC_ROUNDING);
    //putchar('\n');

    //mpfr_t re = y->re;
    //mpfr_t im = y->im;

    printf("prec = %ld\n", x->re->_mpfr_prec);
    printf("sign = %d\n",  x->re->_mpfr_sign);
    printf("exp  = %ld\n", x->re->_mpfr_exp);
    printf("d    = %p\n",  x->re->_mpfr_d);

    printf("prec = %ld\n", y->re->_mpfr_prec);
    printf("sign = %d\n",  y->re->_mpfr_sign);
    printf("exp  = %ld\n", y->re->_mpfr_exp);
    printf("d    = %p\n",  y->re->_mpfr_d);

printf("sylvan_mpc.c mpc_compare() = %d \n", mpc_cmp(x,y));

    return !mpc_cmp(x, y);  // mpc_cmp == 0 if x == y
}


/**
 * Operation "plus" for two mpc MTBDDs
 * Interpret partial function as "0"
 */
TASK_IMPL_2(MTBDD, mpc_op_plus, MTBDD*, pa, MTBDD*, pb)
{
    MTBDD a = *pa, b = *pb;

printf("sylvan_mpc.c mpc_op_plus()\n");

    // Check for partial functions
    if (a == mtbdd_false) return b;
    if (b == mtbdd_false) return a;

    // If both leaves, compute plus
    if (mtbdd_isleaf(a) && mtbdd_isleaf(b)) {

printf("sylvan_mpc.c mtbdd_gettype(a) = %d\n", mtbdd_gettype(a));

        assert(mtbdd_gettype(a) == MPC_TYPE && mtbdd_gettype(b) == MPC_TYPE);

        mpc_ptr ma = (mpc_ptr)mtbdd_getvalue(a);
        mpc_ptr mb = (mpc_ptr)mtbdd_getvalue(b);

        mpc_t x;
        mpc_init2(x, MPC_PRECISION);
        mpc_add(x, ma, mb, MPC_ROUNDING);

    printf("prec = %ld\n", x->re->_mpfr_prec);
    printf("sign = %d\n",  x->re->_mpfr_sign);
    printf("exp  = %ld\n", x->re->_mpfr_exp);
    printf("d    = %p\n",  x->re->_mpfr_d);

        MTBDD result = mtbdd_mpc((mpc_ptr)x);
        mpc_clear(x);
        return result;
    }

    // Commutative, so swap a,b for better cache performance
    if (a < b) {
        *pa = b;
        *pb = a;
    }

    return mtbdd_invalid;
}

/**
 * Operation "times" for two mpc MTBDDs
 * Interpret partial function as "0"
 */
TASK_IMPL_2(MTBDD, mpc_op_times, MTBDD*, pa, MTBDD*, pb)
{
    MTBDD a = *pa, b = *pb;

    /* Check for partial functions and for Boolean (filter) */
    if (a == mtbdd_false || b == mtbdd_false) return mtbdd_false;

    // Check for partial functions
    if (a == mtbdd_true) return b;
    if (b == mtbdd_true) return a;

    // If both leaves, compute multiplication
    if (mtbdd_isleaf(a) && mtbdd_isleaf(b)) {

        assert(mtbdd_gettype(a) == MPC_TYPE && mtbdd_gettype(b) == MPC_TYPE);

        mpc_ptr ma = (mpc_ptr)mtbdd_getvalue(a);
        mpc_ptr mb = (mpc_ptr)mtbdd_getvalue(b);

        mpc_t x;
        mpc_init2(x, MPC_PRECISION);
        mpc_mul(x, ma, mb, MPC_ROUNDING);
        MTBDD result = mtbdd_mpc((mpc_ptr)x);
        mpc_clear(x);
        return result;
    }

    // Commutative, so swap a,b for better cache performance <-- Also for times?
    if (a < b) {
        *pa = b;
        *pb = a;
    }

    return mtbdd_invalid;
}

/**
 * Operation "minus" for two mpc MTBDDs
 * Interpret partial function as "0"
 *
TASK_IMPL_2(MTBDD, mpc_op_minus, MTBDD*, pa, MTBDD*, pb)
{
    MTBDD a = *pa, b = *pb;

    // Check for partial functions
    if (a == mtbdd_false) return gmp_neg(b);
    if (b == mtbdd_false) return a;

    // If both leaves, compute plus
    if (mtbdd_isleaf(a) && mtbdd_isleaf(b)) {

        assert(mtbdd_gettype(a) == mpc_type && mtbdd_gettype(b) == mpc_type);

        mpc_ptr ma = (mpc_ptr)mtbdd_getvalue(a);
        mpc_ptr mb = (mpc_ptr)mtbdd_getvalue(b);

        mpc_t x;
        mpq_init2(x, MPC_PRECISION);
        mpc_sub(x, ma, mb, MPC_ROUNDING);
        MTBDD result = mtbdd_mpc((mpc_ptr)x);
        mpq_clear(x);
        return result;
    }

    return mtbdd_invalid;
}

**
 * Operation "times" for two mpq MTBDDs.
 * One of the parameters can be a BDD, then it is interpreted as a filter.
 * For partial functions, domain is intersection
 *
TASK_IMPL_2(MTBDD, gmp_op_times, MTBDD*, pa, MTBDD*, pb)
{
    MTBDD a = *pa, b = *pb;

    // Check for partial functions and for Boolean (filter) 
    if (a == mtbdd_false || b == mtbdd_false) return mtbdd_false;

    // If one of Boolean, interpret as filter
    if (a == mtbdd_true) return b;
    if (b == mtbdd_true) return a;

    // Handle multiplication of leaves
    if (mtbdd_isleaf(a) && mtbdd_isleaf(b)) {
        assert(mtbdd_gettype(a) == gmp_type && mtbdd_gettype(b) == gmp_type);

        mpq_ptr ma = (mpq_ptr)mtbdd_getvalue(a);
        mpq_ptr mb = (mpq_ptr)mtbdd_getvalue(b);

        // compute result
        mpq_t mres;
        mpq_init(mres);
        mpq_mul(mres, ma, mb);
        MTBDD res = mtbdd_gmp(mres);
        mpq_clear(mres);
        return res;
    }

    // Commutative, so make "a" the lowest for better cache performance
    if (a < b) {
        *pa = b;
        *pb = a;
    }

    return mtbdd_invalid;
}

**
 * Operation "divide" for two mpq MTBDDs.
 * For partial functions, domain is intersection
 *
TASK_IMPL_2(MTBDD, gmp_op_divide, MTBDD*, pa, MTBDD*, pb)
{
    MTBDD a = *pa, b = *pb;

    // Check for partial functions
    if (a == mtbdd_false || b == mtbdd_false) return mtbdd_false;

    // Handle division of leaves
    if (mtbdd_isleaf(a) && mtbdd_isleaf(b)) {
        assert(mtbdd_gettype(a) == gmp_type && mtbdd_gettype(b) == gmp_type);

        mpq_ptr ma = (mpq_ptr)mtbdd_getvalue(a);
        mpq_ptr mb = (mpq_ptr)mtbdd_getvalue(b);

        // compute result
        mpq_t mres;
        mpq_init(mres);
        mpq_div(mres, ma, mb);
        MTBDD res = mtbdd_gmp(mres);
        mpq_clear(mres);
        return res;
    }

    return mtbdd_invalid;
}

**
 * Operation "min" for two mpq MTBDDs.
 *
TASK_IMPL_2(MTBDD, gmp_op_min, MTBDD*, pa, MTBDD*, pb)
{
    MTBDD a = *pa, b = *pb;

    // Handle partial functions
    if (a == mtbdd_false) return b;
    if (b == mtbdd_false) return a;

    // Handle trivial case
    if (a == b) return a;

    // Compute result for leaves
    if (mtbdd_isleaf(a) && mtbdd_isleaf(b)) {
        assert(mtbdd_gettype(a) == gmp_type && mtbdd_gettype(b) == gmp_type);

        mpq_ptr ma = (mpq_ptr)mtbdd_getvalue(a);
        mpq_ptr mb = (mpq_ptr)mtbdd_getvalue(b);
        int cmp = mpq_cmp(ma, mb);
        return cmp < 0 ? a : b;
    }

    // For cache performance
    if (a < b) {
        *pa = b;
        *pb = a;
    }

    return mtbdd_invalid;
}

**
 * Operation "max" for two mpq MTBDDs.
 *
TASK_IMPL_2(MTBDD, gmp_op_max, MTBDD*, pa, MTBDD*, pb)
{
    MTBDD a = *pa, b = *pb;

    // Handle partial functions
    if (a == mtbdd_false) return b;
    if (b == mtbdd_false) return a;

    // Handle trivial case
    if (a == b) return a;

    // Compute result for leaves
    if (mtbdd_isleaf(a) && mtbdd_isleaf(b)) {
        assert(mtbdd_gettype(a) == gmp_type && mtbdd_gettype(b) == gmp_type);

        mpq_ptr ma = (mpq_ptr)mtbdd_getvalue(a);
        mpq_ptr mb = (mpq_ptr)mtbdd_getvalue(b);
        int cmp = mpq_cmp(ma, mb);
        return cmp > 0 ? a : b;
    }

    // For cache performance
    if (a < b) {
        *pa = b;
        *pb = a;
    }

    return mtbdd_invalid;
}

**
 * Operation "neg" for one mpq MTBDD
 *
TASK_IMPL_2(MTBDD, gmp_op_neg, MTBDD, dd, size_t, p)
{
    // Handle partial functions
    if (dd == mtbdd_false) return mtbdd_false;

    // Compute result for leaf
    if (mtbdd_isleaf(dd)) {
        assert(mtbdd_gettype(dd) == gmp_type);

        mpq_ptr m = (mpq_ptr)mtbdd_getvalue(dd);

        mpq_t mres;
        mpq_init(mres);
        mpq_neg(mres, m);
        MTBDD res = mtbdd_gmp(mres);
        mpq_clear(mres);
        return res;
    }

    return mtbdd_invalid;
    (void)p;
}

**
 * Operation "abs" for one mpq MTBDD
 *
TASK_IMPL_2(MTBDD, gmp_op_abs, MTBDD, dd, size_t, p)
{
    // Handle partial functions
    if (dd == mtbdd_false) return mtbdd_false;

    // Compute result for leaf
    if (mtbdd_isleaf(dd)) {
        assert(mtbdd_gettype(dd) == gmp_type);

        mpq_ptr m = (mpq_ptr)mtbdd_getvalue(dd);

        mpq_t mres;
        mpq_init(mres);
        mpq_abs(mres, m);
        MTBDD res = mtbdd_gmp(mres);
        mpq_clear(mres);
        return res;
    }

    return mtbdd_invalid;
    (void)p;
}

**
 * The abstraction operators are called in either of two ways:
 * - with k=0, then just calculate "a op b"
 * - with k<>0, then just calculate "a := a op a", k times
 *
TASK_IMPL_3(MTBDD, gmp_abstract_op_plus, MTBDD, a, MTBDD, b, int, k)
{
    if (k==0) {
        return mtbdd_apply(a, b, TASK(gmp_op_plus));
    } else {
        MTBDD res = a;
        for (int i=0; i<k; i++) {
            mtbdd_refs_push(res);
            res = mtbdd_apply(res, res, TASK(gmp_op_plus));
            mtbdd_refs_pop(1);
        }
        return res;
    }
}

TASK_IMPL_3(MTBDD, gmp_abstract_op_times, MTBDD, a, MTBDD, b, int, k)
{
    if (k==0) {
        return mtbdd_apply(a, b, TASK(gmp_op_times));
    } else {
        MTBDD res = a;
        for (int i=0; i<k; i++) {
            mtbdd_refs_push(res);
            res = mtbdd_apply(res, res, TASK(gmp_op_times));
            mtbdd_refs_pop(1);
        }
        return res;
    }
}

TASK_IMPL_3(MTBDD, gmp_abstract_op_min, MTBDD, a, MTBDD, b, int, k)
{
    if (k == 0) {
        return mtbdd_apply(a, b, TASK(gmp_op_min));
    } else {
        // nothing to do: min(a, a) = a
        return a;
    }
}

TASK_IMPL_3(MTBDD, gmp_abstract_op_max, MTBDD, a, MTBDD, b, int, k)
{
    if (k == 0) {
        return mtbdd_apply(a, b, TASK(gmp_op_max));
    } else {
        // nothing to do: max(a, a) = a
        return a;
    }
}

**
 * Convert to Boolean MTBDD, terminals >= value (double) to True, or False otherwise.
 *
TASK_2(MTBDD, gmp_op_threshold_d, MTBDD, a, size_t, svalue)
{
    // Handle partial function
    if (a == mtbdd_false) return mtbdd_false;

    // Compute result
    if (mtbdd_isleaf(a)) {
        assert(mtbdd_gettype(a) == gmp_type);

        double value = *(double*)&svalue;
        mpq_ptr ma = (mpq_ptr)mtbdd_getvalue(a);
        return mpq_get_d(ma) >= value ? mtbdd_true : mtbdd_false;
    }

    return mtbdd_invalid;
}

**
 * Convert to Boolean MTBDD, terminals > value (double) to True, or False otherwise.
 *
TASK_2(MTBDD, gmp_op_strict_threshold_d, MTBDD, a, size_t, svalue)
{
    // Handle partial function
    if (a == mtbdd_false) return mtbdd_false;

    // Compute result
    if (mtbdd_isleaf(a)) {
        assert(mtbdd_gettype(a) == gmp_type);

        double value = *(double*)&svalue;
        mpq_ptr ma = (mpq_ptr)mtbdd_getvalue(a);
        return mpq_get_d(ma) > value ? mtbdd_true : mtbdd_false;
    }

    return mtbdd_invalid;
}

TASK_IMPL_2(MTBDD, gmp_threshold_d, MTBDD, dd, double, d)
{
    return mtbdd_uapply(dd, TASK(gmp_op_threshold_d), *(size_t*)&d);
}

TASK_IMPL_2(MTBDD, gmp_strict_threshold_d, MTBDD, dd, double, d)
{
    return mtbdd_uapply(dd, TASK(gmp_op_strict_threshold_d), *(size_t*)&d);
}

**
 * Operation "threshold" for mpq MTBDDs.
 * The second parameter must be a mpq leaf.
 *
TASK_IMPL_2(MTBDD, gmp_op_threshold, MTBDD*, pa, MTBDD*, pb)
{
    MTBDD a = *pa, b = *pb;

    // Check for partial functions
    if (a == mtbdd_false) return mtbdd_false;

    // Handle comparison of leaves
    if (mtbdd_isleaf(a)) {
        assert(mtbdd_gettype(a) == gmp_type);

        mpq_ptr ma = (mpq_ptr)mtbdd_getvalue(a);
        mpq_ptr mb = (mpq_ptr)mtbdd_getvalue(b);
        int cmp = mpq_cmp(ma, mb);
        return cmp >= 0 ? mtbdd_true : mtbdd_false;
    }

    return mtbdd_invalid;
}

**
 * Operation "strict threshold" for mpq MTBDDs.
 * The second parameter must be a mpq leaf.
 *
TASK_IMPL_2(MTBDD, gmp_op_strict_threshold, MTBDD*, pa, MTBDD*, pb)
{
    MTBDD a = *pa, b = *pb;

    // Check for partial functions
    if (a == mtbdd_false) return mtbdd_false;

    // Handle comparison of leaves
    if (mtbdd_isleaf(a)) {
        assert(mtbdd_gettype(a) == gmp_type);

        mpq_ptr ma = (mpq_ptr)mtbdd_getvalue(a);
        mpq_ptr mb = (mpq_ptr)mtbdd_getvalue(b);
        int cmp = mpq_cmp(ma, mb);
        return cmp > 0 ? mtbdd_true : mtbdd_false;
    }

    return mtbdd_invalid;
}

// TODO in another file

**
 * Operation "times" for two mpc MTBDDs (leafs are complex numbers).
 * One of the parameters can be a BDD, then it is interpreted as a filter.
 * For partial functions, domain is intersection
 *
TASK_IMPL_2(MTBDD, gmp_op_times, MTBDD*, pa, MTBDD*, pb)
{
    MTBDD a = *pa, b = *pb;

    // Check for partial functions and for Boolean (filter)
    if (a == mtbdd_false || b == mtbdd_false) return mtbdd_false;

    // If one of Boolean, interpret as filter
    if (a == mtbdd_true) return b;
    if (b == mtbdd_true) return a;

    // Handle multiplication of leaves
    if (mtbdd_isleaf(a) && mtbdd_isleaf(b)) {
        assert(mtbdd_gettype(a) == gmp_type && mtbdd_gettype(b) == gmp_type);

        mpc_ptr ma = (mpc_ptr)mtbdd_getvalue(a);
        mpc_ptr mb = (mpc_ptr)mtbdd_getvalue(b);

        // compute result
        mpf_t mres[6];
        mpf_init(mres);
        mpf_mul(mres[0], ma->real, mb->real); // TODO naming vars
        mpf_mul(mres[1], ma->real, mb->imag);
        mpf_mul(mres[2], ma->imag, mb->real);
        mpf_mul(mres[3], ma->imag, mb->imag);

        mpf_min(mres[4], mres[0], mres[3]); // (a + ib)(c + id) = ac-bd + i(bc+ad)
        mpf_plus(mres[5], mres[1], mres[2]);

        MTBDD res = mtbdd_gmp(mres); // TODO how to allocate memory?
        for(int i=0; i<6; i++)
            mpq_clear(mres[i]);
        return res; // node in decision diagram
    }

    // Commutative, so make "a" the lowest for better cache performance
    if (a < b) {
        *pa = b;
        *pb = a;
    }

    return mtbdd_invalid;
}

**
 * Multiply <a> and <b>, and abstract variables <vars> using summation.
 * This is similar to the "and_exists" operation in BDDs.
 *
TASK_IMPL_3(MTBDD, gmp_and_abstract_plus, MTBDD, a, MTBDD, b, MTBDD, v)
{
    // Check terminal cases

    // If v == true, then <vars> is an empty set
    if (v == mtbdd_true) return mtbdd_apply(a, b, TASK(gmp_op_times));

    // Try the times operator on a and b
    MTBDD result = CALL(gmp_op_times, &a, &b);
    if (result != mtbdd_invalid) {
        
        // Times operator successful, store reference (for garbage collection)
        mtbdd_refs_push(result);

        // ... and perform abstraction
        result = mtbdd_abstract(result, v, TASK(gmp_abstract_op_plus));
        mtbdd_refs_pop(1);

        // Note that the operation cache is used in mtbdd_abstract
        return result;
    }

    // Maybe perform garbage collection
    sylvan_gc_test();

    // Count operation
    sylvan_stats_count(MTBDD_AND_ABSTRACT_PLUS);

    // Check cache. Note that we do this now, since the times operator might swap a and b (commutative)
    if (cache_get3(CACHE_MTBDD_AND_ABSTRACT_PLUS, a, b, v, &result)) {
        sylvan_stats_count(MTBDD_AND_ABSTRACT_PLUS_CACHED);
        return result;
    }

    // Now, v is not a constant, and either a or b is not a constant

    // Get top variable
    int la = mtbdd_isleaf(a);
    int lb = mtbdd_isleaf(b);
    mtbddnode_t na = la ? 0 : MTBDD_GETNODE(a);
    mtbddnode_t nb = lb ? 0 : MTBDD_GETNODE(b);
    uint32_t va = la ? 0xffffffff : mtbddnode_getvariable(na);
    uint32_t vb = lb ? 0xffffffff : mtbddnode_getvariable(nb);
    uint32_t var = va < vb ? va : vb;

    mtbddnode_t nv = MTBDD_GETNODE(v);
    uint32_t vv = mtbddnode_getvariable(nv);

    if (vv < var) {
        // Recursive, then abstract result
        result = CALL(gmp_and_abstract_plus, a, b, node_gethigh(v, nv));
        mtbdd_refs_push(result);
        result = mtbdd_apply(result, result, TASK(gmp_op_plus));
        mtbdd_refs_pop(1);
    } else {
        // Get cofactors
        MTBDD alow, ahigh, blow, bhigh;
        alow  = (!la && va == var) ? node_getlow(a, na)  : a;
        ahigh = (!la && va == var) ? node_gethigh(a, na) : a;
        blow  = (!lb && vb == var) ? node_getlow(b, nb)  : b;
        bhigh = (!lb && vb == var) ? node_gethigh(b, nb) : b;

        if (vv == var) {
            // Recursive, then abstract result
            mtbdd_refs_spawn(SPAWN(gmp_and_abstract_plus, ahigh, bhigh, node_gethigh(v, nv)));
            MTBDD low = mtbdd_refs_push(CALL(gmp_and_abstract_plus, alow, blow, node_gethigh(v, nv)));
            MTBDD high = mtbdd_refs_push(mtbdd_refs_sync(SYNC(gmp_and_abstract_plus)));
            result = CALL(mtbdd_apply, low, high, TASK(gmp_op_plus));
            mtbdd_refs_pop(2);
        } else { // vv > v
            // Recursive, then create node
            mtbdd_refs_spawn(SPAWN(gmp_and_abstract_plus, ahigh, bhigh, v));
            MTBDD low = mtbdd_refs_push(CALL(gmp_and_abstract_plus, alow, blow, v));
            MTBDD high = mtbdd_refs_sync(SYNC(gmp_and_abstract_plus));
            mtbdd_refs_pop(1);
            result = mtbdd_makenode(var, low, high);
        }
    }

    // Store in cache
    if (cache_put3(CACHE_MTBDD_AND_ABSTRACT_PLUS, a, b, v, result)) {
        sylvan_stats_count(MTBDD_AND_ABSTRACT_PLUS_CACHEDPUT);
    }

    return result;
}

**
 * Multiply <a> and <b>, and abstract variables <vars> by taking the maximum.
 *
TASK_IMPL_3(MTBDD, gmp_and_abstract_max, MTBDD, a, MTBDD, b, MTBDD, v)
{
    // Check terminal cases 

    // If v == true, then <vars> is an empty set
    if (v == mtbdd_true) return mtbdd_apply(a, b, TASK(gmp_op_times));

    // Try the times operator on a and b
    MTBDD result = CALL(gmp_op_times, &a, &b);
    if (result != mtbdd_invalid) {

        // Times operator successful, store reference (for garbage collection)
        mtbdd_refs_push(result);
        
        // ... and perform abstraction
        result = mtbdd_abstract(result, v, TASK(gmp_abstract_op_max));
        mtbdd_refs_pop(1);

        // Note that the operation cache is used in mtbdd_abstract
        return result;
    }

    // Now, v is not a constant, and either a or b is not a constant

    // Get top variable
    int la = mtbdd_isleaf(a);
    int lb = mtbdd_isleaf(b);
    mtbddnode_t na = la ? 0 : MTBDD_GETNODE(a);
    mtbddnode_t nb = lb ? 0 : MTBDD_GETNODE(b);
    uint32_t va = la ? 0xffffffff : mtbddnode_getvariable(na);
    uint32_t vb = lb ? 0xffffffff : mtbddnode_getvariable(nb);
    uint32_t var = va < vb ? va : vb;

    mtbddnode_t nv = MTBDD_GETNODE(v);
    uint32_t vv = mtbddnode_getvariable(nv);

    while (vv < var) {
        // we can skip variables, because max(r,r) = r
        v = node_high(v, nv);
        if (v == mtbdd_true) return mtbdd_apply(a, b, TASK(gmp_op_times));
        nv = MTBDD_GETNODE(v);
        vv = mtbddnode_getvariable(nv);
    }

    // Maybe perform garbage collection
    sylvan_gc_test();

    // Count operation
    sylvan_stats_count(MTBDD_AND_ABSTRACT_MAX);

    // Check cache. Note that we do this now, since the times operator might swap a and b (commutative)
    if (cache_get3(CACHE_MTBDD_AND_ABSTRACT_MAX, a, b, v, &result)) {
        sylvan_stats_count(MTBDD_AND_ABSTRACT_MAX_CACHED);
        return result;
    }

    // Get cofactors
    MTBDD alow, ahigh, blow, bhigh;
    alow  = (!la && va == var) ? node_getlow(a, na)  : a;
    ahigh = (!la && va == var) ? node_gethigh(a, na) : a;
    blow  = (!lb && vb == var) ? node_getlow(b, nb)  : b;
    bhigh = (!lb && vb == var) ? node_gethigh(b, nb) : b;

    if (vv == var) {
        // Recursive, then abstract result
        mtbdd_refs_spawn(SPAWN(gmp_and_abstract_max, ahigh, bhigh, node_gethigh(v, nv)));
        MTBDD low = mtbdd_refs_push(CALL(gmp_and_abstract_max, alow, blow, node_gethigh(v, nv)));
        MTBDD high = mtbdd_refs_push(mtbdd_refs_sync(SYNC(gmp_and_abstract_max)));
        result = CALL(mtbdd_apply, low, high, TASK(gmp_op_max));
        mtbdd_refs_pop(2);

    } else { // vv > v
        // Recursive, then create node
        mtbdd_refs_spawn(SPAWN(gmp_and_abstract_max, ahigh, bhigh, v));
        MTBDD low = mtbdd_refs_push(CALL(gmp_and_abstract_max, alow, blow, v));
        MTBDD high = mtbdd_refs_sync(SYNC(gmp_and_abstract_max));
        mtbdd_refs_pop(1);
        result = mtbdd_makenode(var, low, high);
    }

    // Store in cache
    if (cache_put3(CACHE_MTBDD_AND_ABSTRACT_MAX, a, b, v, result)) {
        sylvan_stats_count(MTBDD_AND_ABSTRACT_MAX_CACHEDPUT);
    }

    return result;
}

*/

