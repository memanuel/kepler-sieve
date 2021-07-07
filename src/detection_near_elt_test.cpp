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
#include "OrbitalElement.hpp"
    using ks::OrbitalElement;
    using ks::print_orbital_element;
#include "DetectionTime.hpp"
    using ks::DetectionTime;
    using ks::DetectionTimeTable;
#include "Detection.hpp"
    using ks::Detection;
    using ks::DetectionTable;
#include "CandidateElement.hpp"
    using ks::CandidateElement;
#include "db_utils.hpp"
    using ks::db_conn_type;
    using ks::get_db_conn;
    using ks::sp_run;
#include "utils.hpp"
    using ks::report_test;

// *****************************************************************************
// Declare functions defined in this module
int main(int argc, char* argv[]);
void test_all();
void test_load_detection(db_conn_type& conn);

// *****************************************************************************
void test_load_detection(DetectionTable& dt)
{
    Detection d = dt[10];
    print("\nExample Detection:\n");
    print("d.detection_id = {:d}\n", d.detection_id);
    print("d.sky_patch_id = {:d}\n", d.sky_patch_id);
    print("d.time_id = {:d}\n", d.time_id);
    print("d.ux      = {:+8.6f}\n", d.ux);
    print("d.uy      = {:+8.6f}\n", d.uy);
    print("d.uz      = {:+8.6f}\n", d.uz);
    print("d.mag     = {:8.4f}\n", d.mag);
}


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
    dt.load();
    print("Loaded DetectionTable with detection_id in [{:d}, {:d}).\n", d0, d1);
    // test_load_detection(dt);

    // Build orbital elements for asteroid 3 (Juno), which has 8 hits
    // Values copy / pasted from CALL KS.GetAsteroidElements(3, 4, 59000, 59000);
    OrbitalElement elt
    {
        .mjd    =59000.0,
        .a      =2.6682852999999986,
        .e      =0.25693643000000027,
        .inc    =0.22673642125828372,
        .Omega  =2.9644675653853,
        .omega  =-1.9536135288017569,
        .f      =46.517431678528794,
        .M      = 46.171557086433715
    };
    print("\nConstructed OrbitalElement for Juno @ {:f}\n", elt.mjd);
    print_orbital_element(elt);

    // Build CandidateElement for these elements
    CandidateElement celt(elt, dtt);
    print("\nConstructed CandidateOrbitalElement from OrbitalElement.\n");

    // Close DB connection
    conn->close();
}

// *****************************************************************************
int main(int argc, char* argv[])
{
    test_all();
}
