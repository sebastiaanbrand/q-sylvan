/***************************************************************

Complex number defnitions and routines for QMDD using
doubles for the real and imaginary part of a complex number.

January 28, 2008
Michael Miller
University of Victoria
Victoria, BC 
CANADA V8W 3P6
mmiller@cs.uvic.ca

****************************************************************/

#include <stdio.h>

#include "sylvan_qdd_complex.h"
#include "util/cmap.h"


/**********************************************
Compute trig functions for angle factor*Pi/div
Note use of cosl and sinl for long double computation
**********************************************/
#define qdd_cos(fac,div) cosl((long double)(fac)*Pi/(long double)(div))
#define qdd_sin(fac,div) sinl((long double)(fac)*Pi/(long double)(div))


static long double Pi;    // set value of global Pi


int     SIZE;
cmap_t *ctable;
cmap_t *ctable_old;



// return complex conjugate
complex_t
Conj (complex_t c)
{
    c.i = -c.i;
    return c;
}

// make a complex value
complex_t
Cmake (double r, double i)
{
    complex_t c;
    c.r = r;
    c.i = i;
    return (c);
}

complex_t
CmakeAngle (double theta)
{
    complex_t c;
    c.r = cos(theta);
    c.i = sin(theta);
    return c;
}

complex_t
CmakeOne (void)
{
    return Cmake (1.0, 0.0);
}

complex_t
CmakeZero (void)
{
    return Cmake (0.0, 0.0);
}

complex_t
CmakeMOne (void)
{
    return Cmake (-1.0, 0.0);
}

void
Cprint(complex_t c)
{
    if(c.r >= 0)
        printf(" ");
    printf("%.60f", c.r);
    if (c.i > 0)
        printf("+%fi", c.i);
    if (c.i < 0)
        printf("%fi", c.i);
}

void
Cprint_bitvalues(AMP i)
{
    print_bitvalues(ctable, i);
}

// returns the complex number equal to (a+b*sqrt(2))/c
// required to be compatible with quadratic irrational-based
// complex number package
double
Qmake (int a, int b, int c)
{
    return (((float) a + ((float) b) * sqrt (2.0)) / (float) (c));
}


// lookup a complex value in the complex value table
// if not found add it
// this routine uses linear searching
AMP
Clookup (complex_t c)
{
    uint64_t res;
    cmap_find_or_put(ctable, &c, &res);
    return (AMP) res;
}

complex_t
Cvalue (AMP i)
{
    complex_t * res;
    res = cmap_get(ctable, i);
    return *res;
}

complex_t
Cvalue_old (AMP i)
{
    complex_t * res;
    res = cmap_get(ctable_old, i);
    return *res;
}



bool comp_exact_equal(complex_t a, complex_t b)
{
    return (a.r == b.r && a.i == b.i);
}

bool comp_approx_equal(complex_t a, complex_t b)
{
    return comp_epsilon_close(a, b, TOLERANCE);
}

bool comp_epsilon_close(complex_t a, complex_t b, double epsilon)
{
    return ( (fabs(a.r - b.r) < epsilon) && (fabs(a.i - b.i) < epsilon) );
}



/* Arithmetic operations on AMPs */

AMP
amp_abs(AMP a)
{
    // special cases
    if (a == C_ZERO || a == C_ONE) return a;

    complex_t ca, cr;
    AMP res;

    ca = Cvalue(a);
    cr = comp_abs(ca);

    res = Clookup(cr);
    return res;
}

AMP
amp_neg(AMP a)
{
    // special cases
    if (a == C_ZERO) return C_ZERO;

    complex_t ca, cr;
    AMP res;

    ca = Cvalue(a);
    cr = comp_neg(ca);

    res = Clookup(cr);
    return res;
}

AMP
amp_add(AMP a, AMP b)
{
    // special cases
    if (a == C_ZERO) return b;
    if (b == C_ZERO) return a;

    complex_t ca, cb, cr;
    AMP res;

    ca = Cvalue(a);
    cb = Cvalue(b);
    cr = comp_add(ca, cb);

    res = Clookup(cr);
    return res;
}

AMP
amp_sub(AMP a, AMP b)
{
    // special cases
    if (b == C_ZERO) return a;
    if (a == C_ZERO) return amp_neg(b);

    complex_t ca, cb, cr;
    AMP res;

    ca = Cvalue(a);
    cb = Cvalue(b);
    cr = comp_sub(ca, cb);

    res = Clookup(cr);
    return res;
}

AMP
amp_mul(AMP a, AMP b)
{
    // special cases
    if (a == C_ONE) return b;
    if (b == C_ONE) return a;
    if (a == C_ZERO || b == C_ZERO) return C_ZERO;

    complex_t ca, cb, cr;
    AMP res;

    ca = Cvalue(a);
    cb = Cvalue(b);
    cr = comp_mul(ca, cb);

    res = Clookup(cr);
    return res;
}

AMP
amp_div(AMP a, AMP b)
{
    // special cases
    if (a == b)      return C_ONE;
    if (a == C_ZERO) return C_ZERO;
    if (b == C_ONE)  return a;

    complex_t ca, cb, cr;
    AMP res;

    ca = Cvalue(a);
    cb = Cvalue(b);
    cr = comp_div(ca, cb);

    res = Clookup(cr);
    return res;
}



/* Arithmetic operations on complex structs */

complex_t
comp_abs(complex_t a)
{
    complex_t res;
    res.r = sqrt( (a.r*a.r) + (a.i*a.i) );
    res.i = 0.0;
    return res;
}

complex_t
comp_neg(complex_t a)
{
    complex_t res;
    res.r = -a.r;
    res.i = -a.i;
    return res;
}

complex_t
comp_add(complex_t a, complex_t b)
{
    complex_t res;
    res.r = a.r + b.r;
    res.i = a.i + b.i;
    return res;
}

complex_t
comp_sub(complex_t a, complex_t b)
{
    complex_t res;
    res.r = a.r - b.r;
    res.i = a.i - b.i;
    return res;
}

complex_t
comp_mul(complex_t a, complex_t b)
{
    complex_t res;
    res.r = a.r * b.r - a.i * b.i;
    res.i = a.r * b.i + a.i * b.r;
    return res;
}

complex_t
comp_div(complex_t a, complex_t b)
{
    complex_t res;
    double denom;
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


void
init_amplitude_table(size_t size)
{
    SIZE   = size;
    ctable = cmap_create(size);
    
    // TODO: treat 0 and 1 seperately and don't put them in table.
    C_ONE  = Clookup(CmakeOne());
    C_ZERO = Clookup(CmakeZero());
    
    // initialize 2x2 gates (complex values from gates currently stored in 
    // same table as complex amplitude values)
    uint32_t k;

    k = GATEID_I;
    gates[k][0] = C_ONE;  gates[k][1] = C_ZERO;
    gates[k][2] = C_ZERO; gates[k][3] = C_ONE;

    
    k = GATEID_X;
    gates[k][0] = C_ZERO; gates[k][1] = C_ONE;
    gates[k][2] = C_ONE;  gates[k][3] = C_ZERO;

    
    k = GATEID_Y;
    gates[k][0] = C_ZERO; gates[k][1] = Clookup(Cmake(0.0, -1.0));
    gates[k][2] = Clookup(Cmake(0.0, 1.0));  gates[k][3] = C_ZERO;

    
    k = GATEID_Z;
    gates[k][0] = C_ONE;  gates[k][1] = C_ZERO;
    gates[k][2] = C_ZERO; gates[k][3] = Clookup(Cmake(-1.0, 0.0));

    k = GATEID_H;
    gates[k][0] = gates[k][1] = gates[k][2] = Clookup(Cmake(Qmake(0,1,2),0));
    gates[k][3] = Clookup(Cmake(Qmake(0,-1,2),0));

    k = GATEID_S;
    gates[k][0] = C_ONE;  gates[k][1] = C_ZERO;
    gates[k][2] = C_ZERO; gates[k][3] = Clookup(Cmake(0.0, 1.0));

    k = GATEID_T;
    gates[k][0] = C_ONE;  gates[k][1] = C_ZERO;
    gates[k][2] = C_ZERO; gates[k][3] = Clookup(Cmake(1.0/sqrt(2.0), 1.0/sqrt(2.0)));

    k = GATEID_Tdag;
    gates[k][0] = C_ONE;  gates[k][1] = C_ZERO;
    gates[k][2] = C_ZERO; gates[k][3] = Clookup(Cmake(1.0/sqrt(2.0), -1.0/sqrt(2.0)));

    k = GATEID_sqrtX;
    gates[k][0] = Clookup(Cmake(0.5, 0.5)); gates[k][1] = Clookup(Cmake(0.5,-0.5));
    gates[k][2] = Clookup(Cmake(0.5,-0.5)); gates[k][3] = Clookup(Cmake(0.5, 0.5));

    k = GATEID_sqrtY;
    gates[k][0] = Clookup(Cmake(0.5, 0.5)); gates[k][1] = Clookup(Cmake(-0.5,-0.5));
    gates[k][2] = Clookup(Cmake(0.5, 0.5)); gates[k][3] = Clookup(Cmake(0.5, 0.5));


    Pi = 2.0 * acos(0.0);

    init_phase_gates(255);
}

void
init_phase_gates(int n)
{
    // amplitude table needs to be initialized already
    
    // add gate R_k to gates table
    // (note that R_0 = I, R_1 = Z, R_2 = S, R_4 = T)
    uint32_t gate_id;
    double angle;
    complex_t cartesian;
    for (int k=0; k<=n; k++) {
        // forward rotation
        angle = 2*Pi / pow(2.0, (double) k);
        cartesian = CmakeAngle(angle);
        gate_id = GATEID_Rk(k);
        gates[gate_id][0] = C_ONE;  gates[gate_id][1] = C_ZERO;
        gates[gate_id][2] = C_ZERO; gates[gate_id][3] = Clookup(cartesian);

        // backward rotation
        angle = -2*Pi / pow(2.0, (double) k);
        cartesian = CmakeAngle(angle);
        gate_id = GATEID_Rk_dag(k);
        gates[gate_id][0] = C_ONE;  gates[gate_id][1] = C_ZERO;
        gates[gate_id][2] = C_ZERO; gates[gate_id][3] = Clookup(cartesian);
    }
}

uint64_t
count_amplitude_table_enries()
{
    return cmap_count_entries(ctable);
}

void
free_amplitude_table()
{
    cmap_free(ctable);
}

void
init_new_empty_table()
{
    // point old to current (full) ctable
    ctable_old = ctable;

    // re-init new (empty) ctable
    init_amplitude_table(SIZE);
}

void
delete_old_table()
{
    // delete  old (full) table
    cmap_free(ctable_old);
}

AMP
move_from_old_to_new(AMP a)
{
    complex_t c = Cvalue_old(a);
    return Clookup(c);
}
