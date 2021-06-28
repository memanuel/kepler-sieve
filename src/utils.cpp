/** @file utils.cpp
 *  @brief Implmentation of general purpose utilities.
 *
 *  @author Michael S. Emanuel
 *  @date 2021-06-24
 */

// *****************************************************************************
// Local dependencies
#include "utils.hpp"

// *****************************************************************************
// Put functions into the namespace ks (for Kepler Sieve)
namespace ks {

// *****************************************************************************
void print_stars()
{
    print("********************************************************************************\n");
}

// *****************************************************************************
void print_newline()
{
    print("\n");
}

// *****************************************************************************
const string test_message(bool is_ok)
{
    if (is_ok)
        return string("PASS");
    else
        return string("FAIL");
}

// *****************************************************************************
void report_test(const string test_name, bool is_ok)
{
    print("{:s}:\n", test_name);
    print("**** {:s} ****\n", test_message(is_ok));
}

// *****************************************************************************
/// Calculate Cartesian squared distance between two 3-vectors
double norm2(const double *v0, const double *v1)
{
    // Get the three distance components out
    double dx = v1[0]-v0[0];
    double dy = v1[1]-v0[1];
    double dz = v1[2]-v0[2];
    return sqr(dx) + sqr(dy) + sqr(dz);
}

// *****************************************************************************
/// Calculate Cartesian squared distance between two 3-vectors
double norm(const double *v0, const double *v1)
{
    // Delegate to norm2 and take square root
    return sqrt(norm2(v0, v1));
}

// *****************************************************************************
}; // namespace
