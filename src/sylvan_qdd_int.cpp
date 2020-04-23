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

#include <stdio.h>
#include <map>

/**********************************************
Compute trig functions for angle factor*Pi/div
Note use of cosl and sinl for long double computation
**********************************************/
#define qdd_cos(fac,div) cosl((long double)(fac)*Pi/(long double)(div))
#define qdd_sin(fac,div) sinl((long double)(fac)*Pi/(long double)(div))

/***********************************************
Tolerance for testing equality of complex values
***********************************************/
static long double Ctol = 1.0e-10;
static long double Pi;    // set value of global Pi


#define Dzero 0.0
#define Ceq(x,y) ((fabs((x.r)-(y.r))<Ctol)&&(fabs((x.i)-(y.i))<Ctol))




struct complex_cmp
{
    bool operator() ( complex_t x, complex_t y ) const {
    	return Ccomp(x,y);
    }
};


std::map<cint, complex_t> Ctable;
std::map<complex_t, cint, complex_cmp> Ctable2;

bool
Ccomp(complex_t x, complex_t y)
{ // TODO: do we need to use x.m and x.a instead?
    if(fabs((x.r)-(y.r)) < Ctol){ // real parts are "equal"
        if(fabs(x.i)-(y.i) < Ctol) { // im parts are "equal"
            return 0;
        }
        else{
            return x.i > y.i;
        }
    }
    else{
        return x.r > y.r;
    }
}

// return complex conjugate
complex_t
Conj (complex_t c)
{
    c.i = -c.i;
    return c;
}

// make a complex value
complex_t
Cmake (long double r, long double i)
{
    complex_t c;
    c.r = r;
    c.i = i;
    return (c);
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
    printf("%Lf", c.r);
    if (c.i > 0)
        printf("+%Lfi", c.i);
    if (c.i < 0)
        printf("%Lfi", c.i);
}

// returns the complex number equal to (a+b*sqrt(2))/c
// required to be compatible with quadratic irrational-based
// complex number package
long double
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
    cint i = 0;

    // TODO: use something better than copy-paste from QMDD code ?

    // First try to look up c
    std::map<complex_t, cint, complex_cmp>::iterator it;
    it = Ctable2.find(c);
    if(it != Ctable2.end()){
        return it->second;
    }

    // If not found, add a new entry to the table.
    if(Ctentries >= 1048576) {
        printf("ERROR, too many complex numbers\n");
        exit(-4);
    }

    Ctentries++;
    i=Ctentries-1;

    complex_t new_c;
    new_c.r = c.r;
    new_c.i = c.i;
    Ctable[i] = new_c;
    Ctable2[new_c] = i;
    

    // what about deleting values?

    return i;
    (void) c;
}

complex_t
Cvalue (cint i)
{
    //complex_t c = (complex_t) { .r = 0, .i = 0, .m = 0, .a = 0  };

    // it might be better to not call this function on CZRO and CONE in the
    // first place.
    //if(i == CZRO) return CmakeZero();
    //if(i == CONE) return CmakeOne();
    // added Ctable[CZRO] = CmakeZero() and Ctable[CONE] = CmakeOne()
    complex_t c = Ctable[i];

    return c;
    (void) i;
}

// computes angle for polar coordinate representation of Cvalue(a)
long double
angle (cint a)
{
    complex_t ca;
    ca = Cvalue (a);
    if (ca.i >= 0 - Ctol)
        return (acos (ca.r / ca.m));
    else return (2 * Pi - acos (ca.r / ca.m));
}

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

    /*int c=a;
     a=b;
     b=c;*/

    /// BETA END */
    if (ca.m > (cb.m + Ctol)) return (1);
    if (cb.m > (ca.m + Ctol)) return (0);

    //CHANGED by pN 120831
    return ((ca.a + Ctol) < cb.a);
    return (0);
}

cint
Cgt_new (cint a, cint b)
{
    complex_t ca, cb;
    if (a == b) return (0);
    ca = Cvalue (a);
    cb = Cvalue (b);
    if ((ca.a + Ctol) < cb.a) return (1);
    return (ca.m > (cb.m + Ctol));
}

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
    if (ai == 1) return bi; // identity cases
    if (bi == 1) return ai;
    if (ai == 0 || bi == 0) return 0;

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
    long double d;

    if (ai == bi) return (1); // equal case
    if (ai == 0) return (0); // identity cases
    if (bi == 1) return (ai);

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
}

///by PN: returns whether a complex number has norm 1
cint
CUnit (cint a)
{
    /// BETA 121017

    if (a < 2) return a;

    complex_t ca = Cvalue (a);

    return !(ca.m < 1 - Ctol);
}


// TODO: put in header
//typedef complex_t qdd_matrix[2][2];

//qdd_matrix Nm, Vm, VPm, Sm, Rm, Hm, Zm, ZEROm, Qm;



void
qdd_complex_init()
{
    //complex_t v, vc;

    // Clear any existing content
    Ctable.clear();
    Ctable2.clear();
    Ctentries = 0;

    // Set complex(0) and complex(1) to indices 0 and 1 in Ctable
    Ctable[C_ZERO] = CmakeZero();  //Ctable[0] = 0.0 + 0.0i
    Ctable[C_ONE]  = CmakeOne();   //Ctable[1] = 1.0 + 0.0i
    Ctable2[CmakeZero()] = C_ZERO; //Ctable2[0.0 + 0.0i] = 0
    Ctable2[CmakeOne()]  = C_ONE;  //Ctable2[1.0 + 0.0i] = 1
    Ctentries = 2;

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


    Pi = 2.0 * acos(0.0);
}
  
