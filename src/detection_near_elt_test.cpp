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
#include "DetectionTime.hpp"
    using ks::DetectionTime;
    using ks::DetectionTimeTable;
#include "Detection.hpp"
    using ks::Detection;
    using ks::DetectionTable;
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
    int d1 = 1000000;

    // Initialize DetectionTimeTable
    DetectionTimeTable dtt = DetectionTimeTable();
    print("Loaded DetectionTimeTable with {:d} detection times.\n", dtt.N());

    // Initialize DetectionTable
    DetectionTable dt = DetectionTable(d0, d1);
    // dt.load(conn, progbar);
    dt.load();
    print("Loaded DetectionTable with detection_id in [{:d}, {:d}).\n", d0, d1);

    // Close DB connection
    conn->close();
}

// *****************************************************************************
int main(int argc, char* argv[])
{
    test_all();
}
