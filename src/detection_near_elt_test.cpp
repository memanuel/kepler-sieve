/** @file detection_near_elt_test.cpp
 *  @brief Test harness for functions used in finding detections near candidate orbital elements.
 *
 *  @author Michael S. Emanuel
 *  @date 2021-07-06
 * 
 * Example call:
 * ./detection_near_elt_test.x
 * ****************************************************************************/

// *****************************************************************************
// Library dependencies
#include <fmt/format.h>
    using fmt::print;

// Local dependencies
#include "db_utils.hpp"
    using ks::db_conn_type;
    using ks::get_db_conn;
    using ks::sp_run;

#include "utils.hpp"
    using ks::report_test;

#include "Detection.hpp"
    using ks::Detection;
    using ks::DetectionTable;
    using ks::DetectionTime;
    using ks::DetectionTimeTable;

#include "SkyPatchNeighbor.hpp"
    using ks::SkyPatch;

// *****************************************************************************
// Declare functions defined in this module
void test_all();
int main(int argc, char* argv[]);

// *****************************************************************************
void test_all()
{
    // Establish DB connection
    db_conn_type conn = get_db_conn();

    // Inputs to build DetectionTable
    int d0 = 0;
    int d1 = 1000;
    bool progbar = true;

    // Initialize DetectionTimeTable
    DetectionTimeTable dtt = DetectionTimeTable(conn);
    print("Loaded DetectionTimeTable.\n");

    // Initialize DetectionTable
    DetectionTable dt = DetectionTable(conn, d0, d1, progbar);
    print("Loaded DetectionTable.\n");

    // Close DB connection
    conn->close();
}

// *****************************************************************************
int main(int argc, char* argv[])
{
    test_all();
}
