/** @file utils.cpp
 *  @brief General purpose utilities.
 *
 *  @author Michael S. Emanuel
 *  @date 2021-06-24
 */

// *****************************************************************************
#pragma once

// *****************************************************************************
// Library dependencies
#include <cmath>
#include <string>
    using std::string;
#include <iostream>
    using std::cout;
    using std::flush;
#include <fmt/format.h>
#include <fmt/chrono.h>
    using fmt::print;
    using fmt::format;

// *****************************************************************************
namespace ks {

// *****************************************************************************
/// Print a row of 80 stars to the console
void print_stars();

// *****************************************************************************
/// Print a single newline to the console
void print_newline();

// *****************************************************************************
/// Convert a time in seconds to a string in the form hh:mm:ss.s
string time2hms(double t);

// *****************************************************************************
/// Flush buffer so text on console appears
void flush_console();

// *****************************************************************************
/// Message describing test results (PASS or FAIL)
const string test_message(bool is_ok);

// *****************************************************************************
/// Report results of a test
void report_test(const string test_name, bool is_ok);

// *****************************************************************************
/// Test if two numbers are close together in absolute distance
bool is_close_abs(double x, double y, double tol);

// *****************************************************************************
/// Test if two numbers are close together in relative distance
bool is_close_rel(double x, double y, double rel_tol);

// *****************************************************************************
/// Square a number (double)
inline double sqr(double x)
{
    return x*x;
}

/// Square a number (float)
inline float sqr(float x)
{
    return x*x;
}

/// Square a number (integer)
inline int sqr(int x)
{
    return x*x;
}

// *****************************************************************************
/// Sign of a number (double)
inline double sign(double x)
{
    return x > 0 ? 1.0 : ( x < 0 ? -1.0 : 0.0);
}

/// Sign of a number (float)
inline float sign(float x)
{
    return x > 0 ? 1.0 : ( x < 0 ? -1.0 : 0.0);
}

// *****************************************************************************
/// Calculate Cartesian squared distance between two 3-vectors represented as arrays
double norm2(const double *v0, const double *v1);

// *****************************************************************************
/// Calculate Cartesian squared distance between two 3-vectors represented as arrays
double norm(const double *v0, const double *v1);

// *****************************************************************************
} // Namespace ks
