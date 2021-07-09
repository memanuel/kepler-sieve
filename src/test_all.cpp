/** @file   test_all.cpp
 *  @brief  Run full suite of C++ tests for Kepler Sieve project.
 *
 *  @author Michael S. Emanuel
 *  @date   2021-07-01
 * 
 * Example call:
 * ./test_all.x
 */

// *****************************************************************************
// Library dependencies
#include <cstdlib>
    using std::system;
#include <fmt/format.h>
    using fmt::print;

// Local dependencies
#include "utils.hpp"
    using ks::report_test;

// *****************************************************************************
int main()
{
    // Accumulate overall test result
    int result;
    bool is_ok, is_ok_all=true;

    // Run DB test
    result = system("./db_test.x > test/db_test.txt");
    is_ok = (result==0);
    is_ok_all = is_ok_all && is_ok;
    report_test("db_test", is_ok);

    // Run sky patch test
    result = system("./sky_patch_test.x > test/sky_patch_test.txt");
    is_ok = (result==0);
    is_ok_all = is_ok_all && is_ok;
    report_test("sky_patch_test", is_ok);

    // Run detection test
    result = system("./detection_test.x > test/detection_test.txt");
    is_ok = (result==0);
    is_ok_all = is_ok_all && is_ok;
    report_test("detection_test", is_ok);

    // Run rebound test
    result = system("./rebound_test.x > test/rebound_test.txt");
    is_ok = (result==0);
    is_ok_all = is_ok_all && is_ok;
    report_test("rebound_test", is_ok);

    // // TODO
    // // Run planet element test
   
    // Run detection near element test
    result = system("./detection_near_elt_test.x > test/detection_near_elt_test.txt");
    is_ok = (result==0);
    is_ok_all = is_ok_all && is_ok;
    report_test("detection_near_elt_test", is_ok);

    // Report overall test suite results
    report_test("\nKepler Sieve Test Suite:", is_ok_all);

    // Normal program exit; return 1 to signal test failure
    return is_ok_all ? 0 : 1;
}
