#include "sylvan_qdd_complex_mpreal.h"

#include <mpreal.h>
#include <map>

int
mpreal_test()
{
    std::map<int, double> *cpp_map = new std::map<int,double>;
    delete cpp_map;
    return 42;
}