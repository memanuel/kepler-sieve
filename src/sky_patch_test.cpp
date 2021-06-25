// *****************************************************************************
// Included libraries
#include <cmath>
#include <iostream>
#include <boost/format.hpp>

// Local dependencies
#include "utils.h"
#include "sky_patch.h"

// *****************************************************************************
// Names used
using std::cout;
using boost::format;
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

// *****************************************************************************
void test_cube_face()
{
    // Initialize a CubeFace
    CubeFace cf = CubeFace(0);
    // Report results
    print_stars();
    cout << "CubeFace:\n";
    // String attributes
    cout << format("description: %s\n") % cf.description();
    cout << format("code : %s\n") % cf.str();
    cout << format("alpha: %c\n") % cf.alpha();
    cout << format("beta : %c\n") % cf.beta();
    cout << format("gamma: %c\n") % cf.gamma();
    // Integer and float attributes
    cout << format("id   : %d\n") % static_cast<int>(cf.id);
    cout << format("k1   : %d\n") % cf.k1();
    cout << format("k2   : %d\n") % cf.k2();
    cout << format("k3   : %d\n") % cf.k3();
    cout << format("c    :%+3.1f\n") % cf.c();
}

// *****************************************************************************
void test_sky_patch()
{
    // Set SkyPatch integer values for test
    int8_t f_ = 0;
    int16_t i = 512;
    int16_t j = 1024+512;
    // int32_t sky_patch_id = (M2*f_ + M*i + j);
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
    cout << "j  MEAN    MIN     MAX     \n";
    for (int j=0; j<9; j++)
        {cout << format("%d %8.6f %8.6f %8.6f\n") % j % dist_mean[j] % dist_min[j] % dist_max[j];}
}

// *****************************************************************************
int main()
{
    // Test CubeFace class
    // test_cube_face();

    // Test SkyPatch class
    // test_sky_patch();

    // Test SkyPatch neighbor
    spn_type spn = test_sky_patch_neighbor();

    // Test SkyPatch neighbor distance
    test_sky_patch_neighbor_distance(spn);

    // DEBUG
    // SkyPatch sp0 = SkyPatch(0, 0, 0);
    // int32_t spid0 = sp0.id();
    // cout << format("\nSkyPatch sp0 built from spid=%d.\n") % spid0;
    // int di=-1;
    // int dj=0;
    // SkyPatch sp1 = sp0.shift(di, dj);
    // int32_t spid1 = sp1.id();
    // cout << format("\nSkyPatch sp1 is sp0 shifted (%d, %d);  spid1=%d.\n") % di % dj % spid1;
    // cout << format("%s\n") % sp1.str(); 
    // SkyPatch sp1 = SkyPatch_from_id(spid1);
    // cout << format("sp1:\n%s\n") % sp1.str();

    // Normal program exit
    return 0;
}
