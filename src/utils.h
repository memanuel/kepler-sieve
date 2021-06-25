/*****************************************************************************
 * General purpose utilities.
 * 
 * Michael S. Emanuel
 * 2021-06-24
 * ****************************************************************************/

#pragma once

// *****************************************************************************
// Included libraries
#include <cmath>
#include <iostream>
#include <string>
#include <boost/format.hpp>

// *****************************************************************************
// Names used
using std::string;

// *****************************************************************************
// Put all functions into namespace ks
namespace ks {

// *****************************************************************************
//*Print a row of 80 stars to the console
void print_stars();

// *****************************************************************************
//*Print a single newline to the console
void print_newline();

// *****************************************************************************
//*Message describing test results (PASS or FAIL)
const string test_message(bool is_ok);

// *****************************************************************************
//*Report results of a test
void report_test(const string test_name, bool is_ok);

// *****************************************************************************
//*Square a number (double)
inline double sqr(double x)
{
    return x*x;
}

//*Square a number (float)
inline float sqr(float x)
{
    return x*x;
}

//*Square a number (integer)
inline int sqr(int x)
{
    return x*x;
}

// *****************************************************************************
//* Calculate Cartesian squared distance between two 3-vectors
double norm2(const double *v0, const double *v1);

// *****************************************************************************
//* Calculate Cartesian squared distance between two 3-vectors
double norm(const double *v0, const double *v1);

// *****************************************************************************
} // Namespace ks
