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

#include "sylvan_qdd_int.h"
#include "util/cmap.h"


#include <stdio.h>

/**********************************************
Compute trig functions for angle factor*Pi/div
Note use of cosl and sinl for long double computation
**********************************************/
#define qdd_cos(fac,div) cosl((long double)(fac)*Pi/(long double)(div))
#define qdd_sin(fac,div) sinl((long double)(fac)*Pi/(long double)(div))


static long double Pi;    // set value of global Pi


int      SIZE;
cmap_t * ctable;
cmap_t * ctable_old;



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
Cprint_bitvalues(cint i)
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
cint
Clookup (complex_t c)
{
    uint64_t res;
    cmap_find_or_put(ctable, &c, &res);
    return (cint) res;
}

complex_t
Cvalue (cint i)
{
    complex_t * res;
    res = cmap_get(ctable, i);
    return *res;
}

complex_t
Cvalue_old (cint i)
{
    complex_t * res;
    res = cmap_get(ctable_old, i);
    return *res;
}

// computes angle for polar coordinate representation of Cvalue(a)
/*
long double
angle (cint a)
{
    complex_t ca;
    ca = Cvalue (a);
    if (ca.i >= 0 - Ctol)
        return (acos (ca.r / ca.m));
    else return (2 * Pi - acos (ca.r / ca.m));
} */

/*
cint
Cgt (cint a, cint b)
// returns 1 if |a|>|b|
// returns 0 if |b|>|a|
// returns angle(a)<angle(b)
// where angle is the angle in polar coordinate representation
{
    complex_t ca, cb;
    if (a == b) return 0;
    ca = Cvalue (a);
    cb = Cvalue (b);

    /// BETA: 121017
    // returns the smaller nonzero value
    if (a == 0) return (1);
    if (b == 0) return (0);

    //int c=a;
    // a=b;
    // b=c;

    // BETA END 
    if (ca.m > (cb.m + Ctol)) return (1);
    if (cb.m > (ca.m + Ctol)) return (0);

    //CHANGED by pN 120831
    return ((ca.a + Ctol) < cb.a);
    return (0);
} */

/*
cint
Cgt_new (cint a, cint b)
{
    complex_t ca, cb;
    if (a == b) return (0);
    ca = Cvalue (a);
    cb = Cvalue (b);
    if ((ca.a + Ctol) < cb.a) return (1);
    return (ca.m > (cb.m + Ctol));
} */

/*
cint
Clt (cint a, cint b)
// analogous to Cgt
{
    complex_t ca, cb;
    if (a == b) return (0);
    ca = Cvalue (a);
    cb = Cvalue (b);
    if (ca.m < (cb.m + Ctol)) return (1);
    if (cb.m < (ca.m + Ctol)) return (0);
    return ((angle (a) + Ctol) > angle (b));
} */

bool CexactEqual(complex_t a, complex_t b)
{
    return (a.r == b.r && a.i == b.i);
}

bool CapproxEqual(complex_t a, complex_t b)
{
    return CepsilonClose(a, b, TOLERANCE);
}

bool CepsilonClose(complex_t a, complex_t b, double epsilon)
{
    return ( (fabs(a.r - b.r) < epsilon) && (fabs(a.i - b.i) < epsilon) );
}


// basic operations on complex values
// meanings are self-evident from the names
// NOTE arguments are the indices to the values 
// in the complex value table not the values themselves

cint
Cnegative (cint a)
{
    complex_t c = Cvalue (a);
    c.r = -c.r;
    c.i = -c.i;
    return Clookup(c);
}

cint
Cadd (cint ai, cint bi)
{
    complex_t a, b, r;
    int t;

    if (ai == 0) return (bi); // identity cases
    if (bi == 0) return (ai);

    // TODO: dynamic programming
//    if (0 <= (t = cta[ai][bi])) return t; // look in computation table

    a = Cvalue (ai); // if new compute result
    b = Cvalue (bi);
    r.r = a.r + b.r;
    r.i = a.i + b.i;

    t = //cta[ai][bi] = cta[bi][ai] = ///TODO
                    Clookup (r); // save result
    return (t);
}

cint
Csub (cint ai, cint bi)
{
    complex_t a, b, r;
    int t;

    if (bi == 0) return (ai); // identity case

    // TODO: dynamic programming
    //if (0 <= (t = cts[ai][bi])) return t; // look in computation table

    a = Cvalue (ai);  // if new compute result
    b = Cvalue (bi);
    r.r = a.r - b.r;
    r.i = a.i - b.i;

    t //= cts[ai][bi] // TODO
       = Clookup (r); // save result
    return t;
}

cint
Cmul (cint ai, cint bi)
{
    complex_t a, b, r;
    int t;

    // What's best imo is to have these checks as optimization, but have the 
    // code so that even without treating 0 and 1 as special cases it still
    // runs correctly.
    if (ai == C_ONE) return bi; // identity cases
    if (bi == C_ONE) return ai;
    if (ai == C_ZERO || bi == C_ZERO) return C_ZERO;

    // TODO: dynamic programming
    //if (0 <= (t = ctm[ai][bi])) return (t); // look in computation table

    a = Cvalue (ai); // if new compute result
    b = Cvalue (bi);
    r.r = a.r * b.r - a.i * b.i;
    r.i = a.r * b.i + a.i * b.r;

    t //= ctm[ai][bi] = ctm[bi][ai] // TODO
      = Clookup (r); // save result
    return t;
}

cint
CintMul (cint a, cint bi)
{
    complex_t r = Cvalue (bi);
    r.r *= a;
    r.i *= a;
    return Clookup(r);
}

cint
Cdiv (cint ai, cint bi)
{
    complex_t a, b, r;
    int t;
    double d;

    if (ai == bi)     return C_ONE;
    if (ai == C_ZERO) return C_ZERO;
    if (bi == C_ONE)  return ai;

    // TODO: dynamic programming
    //if (0 <= (t = ctd[ai][bi])) return (t); // check computation table

    a = Cvalue (ai); // if new compute result
    b = Cvalue (bi);
    if (b.i == 0.0) {
        r.r = a.r / b.r;
        r.i = a.i / b.r;
    } else {
        d = b.r * b.r + b.i * b.i;
        r.r = (a.r * b.r + a.i * b.i) / d;
        r.i = (a.i * b.r - a.r * b.i) / d;
    }
    t = ///ctd[ai][bi] = // TODO
                    Clookup (r); // save result
    return t;
}

/// by PN: returns the absolut value of a complex number
/*
cint
CAbs (cint a)
{
    int b;
    complex_t r, s;

    if (a < 2) return a; // trivial cases 0/1

    s = Cvalue (a);
    //printf("CAbs: "); Cprint(s); printf(" is ");
    r.r = s.m;
    r.i = 0;
    b = Clookup (r);
    //Cprint(r);   printf("\n");
    return b;
} */

///by PN: returns whether a complex number has norm 1
/*
cint
CUnit (cint a)
{
    /// BETA 121017

    if (a < 2) return a;

    complex_t ca = Cvalue (a);

    return !(ca.m < 1 - Ctol);
} */


// TODO: put in header
//typedef complex_t qdd_matrix[2][2];

//qdd_matrix Nm, Vm, VPm, Sm, Rm, Hm, Zm, ZEROm, Qm;



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

cint
move_from_old_to_new(cint a)
{
    complex_t c = Cvalue_old(a);
    return Clookup(c);
}
