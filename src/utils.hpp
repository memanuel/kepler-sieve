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
#include <concepts>
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
// Create a concept for a numeric type
template<typename T>
concept numeric = std::integral<T> or std::floating_point<T>;

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
/// Square a number (any numeric type)
template<typename T>
requires numeric<T>
inline T sqr(T x)
{
    return x*x;
} 

/// Specialize template of sqr for double
template<>
inline double sqr(double x)
{
    return x*x;
} 

// *****************************************************************************
/// Cube a number (any numeric type)
template<typename T>
requires numeric<T>
inline T cube(T x)
{
    return x*x*x;
} 

/// Specialize template of cube for double
template<>
inline double cube(double x)
{
    return x*x*x;
} 

// *****************************************************************************
/// Sign of a number (any numeric type)
template<typename T>
requires numeric<T>
inline T sign(T x)
{
    return x > 0 ? 1 : ( x < 0 ? -1 : 0);
} 

/// Specialize template of sign for double
template<>
inline double sign(double x)
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
