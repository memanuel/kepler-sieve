/*****************************************************************************
 * Test harness for SkyPatch and sky_patch.
 * 
 * Michael S. Emanuel
 * 2021-06-24
 * ****************************************************************************/

// *****************************************************************************
// Included libraries
#include <cmath>
#include <fmt/format.h>

// Local dependencies
#include "utils.h"
#include "astro_utils.h"
#include "SkyPatchNeighbor.h"

// *****************************************************************************
// Names used
using std::min_element;
using std::max_element;
using fmt::print;
using fmt::format;
using ks::CubeFace;
using ks::SkyPatch;
using ks::sky_patch::N;
using ks::sky_patch::N_sp;
using ks::sky_patch::N_spn;
using ks::SkyPatch_from_id;
using ks::fij2spid;
using ks::SkyPatchNeighbor;
// using ks::sky_patch::spn_type;
// using ks::sky_patch::spnd_type;
// using ks::make_sky_patch_neighbor_table;
// using ks::make_sky_patch_neighbor_dist_table;
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
    print("CubeFace:\n");
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
    print("SkyPatch:\n{}", sp.str());
    // Calculate sky patch integer attributes
    print("\nInteger attributes\n");
    print("f:   {}\n", sp.f.id);
    print("i:   {}\n", sp.i);
    print("j:   {}\n", sp.j);
    print("id   {}\n", sp.id());
    // Test that recovered ID matches
    bool is_ok_id = (sky_patch_id == sp.id());
    string test_name = format("SkyPatch: recovered ID matches calculated ID ({})", sky_patch_id);
    report_test(test_name, is_ok_id);
    
    // Calculate local sky patch coordinates (u, v, w) on unit sphere
    double u = sp.u();
    double v = sp.v();
    double w = sp.w();
    // Display (u, v, w)
    print("\nLocal coordinates of midpoint on unit sphere\n");
    print("u:   {:+10.8f}\n", u);
    print("v:   {:+10.8f}\n", v);
    print("w:   {:+10.8f}\n", w);

    // Calculate global sky patch coordinates (x, y, z) on unit sphere
    double x = sp.x();
    double y = sp.y();
    double z = sp.z();
    // Display (x, y, z)
    print("Global coordinates of midpoint on unit sphere\n");
    print("x:   {:+10.8f}\n", x);
    print("y:   {:+10.8f}\n", y);
    print("z:   {:+10.8f}\n", z);

    // Check that (u, v, w) is really on the sphere
    double r_sph = sqrt(sqr(u) + sqr(v) + sqr(w));
    print("\nRecovered radius of point (u, v, w) on sphere:\n");
    print("r    : {:8.6f}\n", r_sph);
    double err = fabs(1.0 - r_sph);
    print("Error: {:8.6e}\n", err);
    bool is_ok_sphere = (err < 1.0E-15);
    report_test("SkyPatch: point(u, v, w) on unit sphere", is_ok_sphere);

    // Create a second SkyPatch from the integer ID and check that it's the same as the first
    SkyPatch sp2 = SkyPatch_from_id(sky_patch_id);
    print("\nSkyPatch built using this ID {}):\n", sky_patch_id);
    print("({:+8.6f}, {:+8.6f}, {:+8.6f})\n", sp2.u(), sp2.v(), sp2.w());
    print("id:   {}\n", sp2.id());
    bool is_ok_from_id = (sp2.id() == sky_patch_id);
    report_test("SkyPatch: instance built from ID matches input sky_patch_id", is_ok_from_id);

    // Combined test results
    bool is_ok = (is_ok_id && is_ok_sphere && is_ok_from_id);
    report_test("\nSkyPatch: overall test results", is_ok);
}

// *****************************************************************************
void test_sky_patch_neighbor()
{
    // Build the SkyPatchNeighbor table
    print_newline();
    print_stars();
    print("Building SkyPatch neighbors for N = {}...\n", N);
    SkyPatchNeighbor spn = SkyPatchNeighbor();
    print("Completed SkyPatch neighbor table spn.\n");

    // Array pointing to the 9 neighbors of a sky patch
    int32_t *neighbors;

    // Count the number of nonzero neighbors
    int neighbor_count = 0;
    int missing_count = 0;
    for (int spid=0; spid<N_spn; spid++)
    {
        // Get the 9 neibhbors of this sky ppatch
        neighbors = spn[spid];
        // Iterate through the neibhbors, counting the real entries
        for (int j=0; j<9; j++){
            if (spn[j]>= 0) 
            {
                neighbor_count++;
            }
            else
            {
                missing_count++;
            }
        }
    }
    // Report number of nonzero neighbors
    print("SkyPatchNeighbor table has {} valid entries and {} holes.\n", neighbor_count, missing_count);

    // Initialize a starting SkyPatch
    int8_t f = 0;
    int16_t i = 0;
    int16_t j = 1024;
    int32_t spid0 = fij2spid(f,i,j);
    SkyPatch sp0 = SkyPatch_from_id(spid0);

    // Read off neighbors of the selected sky patch
    print("Starting SkyPatch:\n{}", sp0.str());
    print("Neighbors of this SkyPatch:\n");
    // Offset into table for sp0
    neighbors = spn[spid0];
    for (int k=0; k<9; k++)
    {
        // The jth neighbor
        int32_t spid1 = neighbors[k];
        // Only process *real* neighbors with non-negative spids
        if (spid1 >= 0)
        {
            SkyPatch sp1 = SkyPatch_from_id(spid1);
            // print("spid: {}\n", sp1.id());
            print(sp1.str());
        }
    }
}

// *****************************************************************************
void test_sky_patch_neighbor_distance()
{

    // Build the SkyPatchNeighbor table
    print_newline();
    print_stars();
    print("Building SkyPatch neighbors for N = {}...\n", N);
    SkyPatchNeighbor spn = SkyPatchNeighbor();
    print("Completed SkyPatch neighbor table spn.\n");
    // Populate the distance table
    spn.build_neighbor_distance();

    // Array pointing to the 9 neighbors of a sky patch
    int32_t *neighbors;
    // Array pointing to the 9 neighbor distances of a sky patch
    double *dist;

    // Initialize arrays with summary statistics of neighbors by column j
    int dist_count[9];
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

    // Iterate through all the sky patches, accumulating the summary stats
    for (int spid=0; spid<N_sp; spid++)
    {
        // The IDs of the 9 neighbors
        neighbors = spn[spid];
        // The distances to the 9 neighbors of this sky patch
        dist = spn.neighbor_distance(spid);
        // Iterate through the 9 columns
        for (int j=0; j<9; j++)
        {
            // Skip this entry if it's not a real neighbor
            if (neighbors[j]<0) {continue;}
            // Get the distance
            double x = dist[j];
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
    print("di  dj    MEAN    MIN   MAX   COUNT    HOLES\n");
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
            print("{:+d}  {:+d}   {:6.1f} {:6.1f} {:6.1f} {:d} {:d}\n", 
                di, dj, mean, min, max, count, holes);
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
        if (k==4) {continue;}
        holes += (N_sp - dist_count[k]);
        mean += dist_mean[k]/9.0;
        if (dist_min[k] < min) {min=dist_min[k];}
        if (dist_max[k] > max) {max=dist_max[k];}
    }
    // Convert from Cartesian distance to arc seconds
    double mean_sec = dist2sec(mean);
    double min_sec = dist2sec(min);
    double max_sec = dist2sec(max);

    // Report global results
    print("Summary statistics over all non-trivial neighbor interactions.\n");
    print("\nMean distance: {:6.1f} arc seconds\n", mean_sec);
    print("Min  distance: {:6.1f} arc seconds\n", min_sec);
    print("Max  distance: {:6.1f} arc seconds\n", max_sec);
    print("Number of holes: {:d}\n", holes);

    // Test that number of holes is 24
    bool is_ok_holes = (holes == 24);
    report_test("\nSkyPatchNeighbor: total number of holes is 24?", is_ok_holes);

    // Test that max distance between neighbors is not too large
    double thresh_sec = 300.0;
    bool is_ok_max = (max_sec < thresh_sec);
    report_test((format("\nSkyPatchNeighbor: max distance < {:4.0f} arc seconds?", thresh_sec)), is_ok_max);

    // Test symmetry
    // double thresh_sym = 0.01;

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
    test_sky_patch();

    // Test SkyPatch neighbor
    // test_sky_patch_neighbor();

    // Test SkyPatch neighbor distance
    test_sky_patch_neighbor_distance();

    // Normal program exit
    return 0;
}
