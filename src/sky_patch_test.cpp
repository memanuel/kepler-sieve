/*****************************************************************************
 * Test harness for SkyPatch and sky_patch.
 * 
 * Michael S. Emanuel
 * 2021-06-24
 * ****************************************************************************/

// *****************************************************************************
// Included libraries
#include <cmath>
#include <iostream>
#include <fmt/core.h>
#include <boost/format.hpp>

// Local dependencies
#include "utils.h"
#include "astro_utils.h"
#include "sky_patch.h"

// *****************************************************************************
// Names used
using std::cout;
using std::min_element;
using std::max_element;
using boost::format;
using fmt::print;
using ks::CubeFace;
using ks::SkyPatch;
using ks::sky_patch::N;
using ks::sky_patch::N_sp;
using ks::sky_patch::N_spn;
using ks::SkyPatch_from_id;
using ks::fij2spid;
using ks::spn_type;
using ks::spnd_type;
using ks::make_sky_patch_neighbor_table;
using ks::make_sky_patch_neighbor_dist_table;
using ks::sqr;
using ks::print_stars;
using ks::print_newline;
using ks::report_test;
using ks::dist2sec;

// *****************************************************************************
void test_cube_face()
{
    // Initialize a CubeFace
    CubeFace cf = CubeFace(0);
    // Report results
    print_stars();
    cout << "CubeFace:\n";
    // String attributes
    print("description: {} \n", cf.description());
    print("code : {}\n", cf.str());
    print("alpha: {}\n", cf.alpha());
    print("beta : {}\n", cf.beta());
    print("gamma: {}\n", cf.gamma());
    // Integer and float attributes
    print("id   : {}\n", cf.id);
    print("k1   : {}\n", cf.k1());
    print("k2   : {}\n", cf.k2());
    print("k3   : {}\n", cf.k3());
    print("c    :{:+3.1f}\n", cf.c());
}

// *****************************************************************************
void test_sky_patch()
{
    // Set SkyPatch integer values for test
    int8_t f_ = 0;
    int16_t i = 512;
    int16_t j = 1024+512;
    int32_t sky_patch_id = fij2spid(f_, i, j);
    // Initialize a SkyPatch
    SkyPatch sp = SkyPatch(f_, i, j);
    // Report results
    print_newline();
    print_stars();
    cout << "SkyPatch:\n";
    cout << sp.str();
    // Calculate sky patch integer attributes
    cout << "\nInteger attributes\n";
    cout << format("f:   %d\n") % static_cast<int>(sp.f.id);
    cout << format("i:   %d\n") % static_cast<int>(sp.i);
    cout << format("j:   %d\n") % static_cast<int>(sp.j);
    cout << format("id   %d\n") % sp.id();
    // Test that recovered ID matches
    bool is_ok_id = (sky_patch_id == sp.id());
    string test_name = (format("SkyPatch: recovered ID matches calculated ID (%s)") % sky_patch_id).str();
    report_test(test_name, is_ok_id);
    
    // Calculate local sky patch coordinates (u, v, w) on unit sphere
    double u = sp.u();
    double v = sp.v();
    double w = sp.w();
    // Display (u, v, w)
    cout << "\nLocal coordinates of midpoint on unit sphere\n";
    cout << format("u:   %+10.8f\n") % u;
    cout << format("v:   %+10.8f\n") % v;
    cout << format("w:   %+10.8f\n") % w;

    // Calculate global sky patch coordinates (x, y, z) on unit sphere
    double x = sp.x();
    double y = sp.y();
    double z = sp.z();
    // Display (x, y, z)
    cout << "Global coordinates of midpoint on unit sphere\n";
    cout << format("x:   %+10.8f\n") % x;
    cout << format("y:   %+10.8f\n") % y;
    cout << format("z:   %+10.8f\n") % z;

    // Check that (u, v, w) is really on the sphere
    double r_sph = sqrt(sqr(u) + sqr(v) + sqr(w));
    cout << "\nRecovered radius of point (u, v, w) on sphere:\n";
    cout << format("r    : %8.6f\n") % r_sph;
    double err = fabs(1.0 - r_sph);
    cout << format("Error: %8.6e\n") % err;
    bool is_ok_sphere = (err < 1.0E-15);
    report_test("SkyPatch: point(u, v, w) on unit sphere", is_ok_sphere);

    // Create a second SkyPatch from the integer ID and check that it's the same as the first
    SkyPatch sp2 = SkyPatch_from_id(sky_patch_id);
    cout << format("\nSkyPatch built using this ID (%d):\n") % sky_patch_id;
    cout << format("(%+8.6f, %+8.6f, %+8.6f)\n") % sp2.u() % sp2.v() % sp2.w();
    cout << format("id:   %d\n") % sp2.id();
    bool is_ok_from_id = (sp2.id() == sky_patch_id);
    report_test("SkyPatch: instance built from ID matches input sky_patch_id", is_ok_from_id);

    // Combined test results
    bool is_ok = (is_ok_id && is_ok_sphere && is_ok_from_id);
    report_test("\nSkyPatch: overall test results", is_ok);
}

// *****************************************************************************
spn_type test_sky_patch_neighbor()
{
    // Build the SkyPatchNeighbor table
    print_newline();
    print_stars();
    cout << format("Building SkyPatch neighbors for N = %d...\n") % N;
    spn_type spn = make_sky_patch_neighbor_table();
    cout << format("Completed SkyPatch neighbor table spn.\n");

    // Count the number of nonzero neighbors
    int neighbor_count = 0;
    int missing_count = 0;
    for (int i=0; i<N_spn; i++)
    {
        if (spn[i]>= 0) 
        {
            neighbor_count++;
        }
        else
        {
            missing_count++;
        }
    }
    // Report number of nonzero neighbors
    cout << format("SkyPatchNeighbor table has %d valid entries and %d holes.\n") % neighbor_count % missing_count;

    // Initialize a starting SkyPatch
    int8_t f = 0;
    int16_t i = 0;
    int16_t j = 1024;
    int32_t spid0 = fij2spid(f,i,j);
    SkyPatch sp0 = SkyPatch_from_id(spid0);

    // Read off neighbors of first row
    cout << format("Starting SkyPatch:\n%s") % sp0.str();
    cout << format("Neighbors of this SkyPatch:\n");
    // Offset into table for sp0
    int32_t idx0 = spid0*9;
    for (int j=0; j<9; j++)
    {
        // The jth neighbor
        int32_t spid1 = spn[idx0+j];
        // Only process *real* neighbors with non-negative spids
        if (spid1 >= 0)
        {
            SkyPatch sp1 = SkyPatch_from_id(spid1);
            // cout << format("spid: %d\n") % sp1.id();            
            cout << sp1.str();
        }
    }

    // Return the assembled SkyPatch neighbor table for use in the next test
    return spn;
}

// *****************************************************************************
void test_sky_patch_neighbor_distance(spn_type spn)
{
    // Build the SkyPatchNeighborDistance table
    print_newline();
    print_stars();
    cout << format("Building SkyPatch neighbor distance for N = %d...\n") % N;
    spnd_type spnd = make_sky_patch_neighbor_dist_table(spn);

    // Initialize arrays with summary statistics of neighbors by column j
    int dist_count[9];
    // int dist_count_dummy[9];
    double dist_sum[9];
    double dist_mean[9];
    double dist_min[9];
    double dist_max[9];
    // Initialize the 9 columns
    for (int j=0; j<9; j++)
    {
        dist_count[j] = 0;
        dist_sum[j] = 0.0;
        dist_min[j] = 2.0;
        dist_max[j] = 0.0;
    }
    // Iterate through all the rows, accumulating the summary stats
    for (int i=0; i<N_sp; i++)
    {
        // Base index for this row
        int idx0 = i*9;
        // Iterate through the 9 columns
        for (int j=0; j<9; j++)
        {
            // Index of this entry
            int idx = idx0+j;
            // Skip this entry if it's not a real neighbor
            if (spn[idx]<0) {continue;}
            // Get the distance
            double x = spnd[idx];
            // Accumulate the count
            dist_count[j]++;
            // Accumulate the total
            dist_sum[j] += x;
            // Accumulate the min
            if (x<dist_min[j]) {dist_min[j]=x;}
            // Accumulate the max
            if (x>dist_max[j]) {dist_max[j]=x;}
        }   // for j
    }   // for i
    // Calculate the mean from the total and count
    for (int j=0; j<9; j++)
    {
        dist_mean[j] = dist_sum[j] / dist_count[j];
    }

    // Report the summary statistics
    cout << "di  dj   MEAN   MIN    MAX    COUNT    HOLES\n";
    int k=0;    
    for (int di=-1; di<=1; di++)
    {
        for (int dj=-1; dj<=1; dj++)
        {
            // The number of records and holes
            int count = dist_count[k];
            int holes = N_sp - count;
            // Convert distances to arc seconds
            double mean = dist2sec(dist_mean[k]);
            double min = dist2sec(dist_min[k]);
            double max = dist2sec(dist_max[k]);
            cout << format("%+d  %+d   %6.1f %6.1f %6.1f %d %d\n") 
                % di % dj % mean % min % max % count % holes;
            // Increment the column counter k (looping on di and dj but need k for array column index)
            k++;
        }
    }

    // Global summary statistics
    int holes=0;
    double mean=0.0;
    double min=2.0;
    double max=0.0;
    for (int k=0; k<9;k++)
    {
        // Skip the case of k=4, since there is no move
        if (k==0) {continue;}
        holes += (N_sp - dist_count[k]);
        mean += dist_mean[k]/9.0;
        if (dist_min[k] < min) {min=dist_max[k];}
        if (dist_max[k] > max) {max=dist_max[k];}
    }
    // Convert from Cartesian distance to arc seconds
    double mean_sec = dist2sec(mean);
    double min_sec = dist2sec(min);
    double max_sec = dist2sec(max);

    // Report global results
    cout << format("Summary statistics over all non-trivial neighbor interactions.\n");
    cout << format("\nMean distance: %6.1f arc seconds\n") % mean_sec;
    cout << format("Max  distance: %6.1f arc seconds\n") % max_sec;
    cout << format("Number of holes: %d\n") % holes;

    // Test that number of holes is 24
    bool is_ok_holes = (holes == 24);
    report_test("\nSkyPatchNeighbor: total number of holes is 24?", is_ok_holes);

    // Test that max distance between neighbors is not too large
    double thresh_sec = 300.0;
    bool is_ok_max = (max_sec < thresh_sec);
    report_test((format("\nSkyPatchNeighbor: max distance < %d arc seconds?")%thresh_sec).str(), is_ok_max);

    // Test symmetry
    double thresh_sym = 0.01;

    // double min_rook = min_element({dist_mean[1], dist_mean[3], dist_mean[5], dist_mean[7]});
    // double max_rook = max_element({dist_mean[1], dist_mean[3], dist_mean[5], dist_mean[7]});
    // double err_rook = max_rook - min_rook;
    // bool is_ok_sym_rook = err_rook < thresh_sym;

    // double min_diag = min_element({dist_mean[0], dist_mean[2], dist_mean[6], dist_mean[8]});
    // double max_diag = max_element({dist_mean[0], dist_mean[2], dist_mean[6], dist_mean[8]});
    // double err_diag = max_diag - min_diag;
    // bool is_ok_sym_diag = err_diag < thresh_sym;
    // bool is_ok_sym = is_ok_sym_rook && is_ok_sym_diag;
    // report_test("SkyPatchNeighbor: distance symmetric across 'rook' and 'diagonal' style moves?", is_ok_holes);

    bool is_ok = is_ok_holes && is_ok_max;
    report_test("SkyPatchNeighbor: overall test results", is_ok);
}

// *****************************************************************************
int main()
{
    // Test CubeFace class
    test_cube_face();

    // Test SkyPatch class
    // test_sky_patch();

    // Test SkyPatch neighbor
    // spn_type spn = test_sky_patch_neighbor();

    // Test SkyPatch neighbor distance
    // test_sky_patch_neighbor_distance(spn);

    // DEBUG
    // std::cout << "\nBoost version:" << BOOST_LIB_VERSION << '\n';

    // Normal program exit
    return 0;
}
