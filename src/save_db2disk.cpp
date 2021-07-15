/** @file save_db2disk.cpp
 *  @brief Batch program to save output of DB queries of fast cache to disk.
 *
 *  @author Michael S. Emanuel
 *  @date 2021-07-08
 * 
 * Example calls:
 * ./save_db2disk.x --DetectionTime
 * ./save_db2disk.x --DetectionCandidate
 * ./save_db2disk.x --Detection
 * ./save_db2disk.x --BodyVector
 * ./save_db2disk.x --PlanetVector
 * ./save_db2disk.x --PlanetElement
 * ****************************************************************************/

// *****************************************************************************
// Library dependencies
#include <filesystem>
    using std::filesystem::exists;
#include <fmt/format.h>
    using fmt::print;
#include <boost/program_options.hpp>
    namespace po = boost::program_options;

// Local dependencies
#include "constants.hpp"
    using ks::cs::mpd;
#include "utils.hpp"
    using ks::print_stars;
    using ks::report_test;
#include "db_utils.hpp"
    using ks::db_conn_type;
    using ks::get_db_conn;
    using ks::sp_run;
#include "astro_utils.hpp"
    using ks::SolarSystemBody_bv;
#include "MassiveBody.hpp"
    using ks::MassiveBody;
    using ks::MassiveBodyTable;
#include "DetectionTime.hpp"
    using ks::DetectionTime;
    using ks::DetectionTimeTable;
#include "DetectionCandidate.hpp"
    using ks::DetectionCandidate;
    using ks::DetectionCandidateTable;
#include "Detection.hpp"
    using ks::Detection;
    using ks::DetectionTable;
#include "BodyVector.hpp"
    using ks::BodyVector;
#include "PlanetVector.hpp"
    using ks::PlanetVector;
#include "PlanetElement.hpp"
    using ks::PlanetElement;

// *****************************************************************************
// Declare functions
void save_MassiveBody(db_conn_type& conn);
void save_DetectionTime(db_conn_type& conn);
void save_DetectionCandidate(db_conn_type& conn);
void save_Detection(db_conn_type& conn);
void save_BodyVector(db_conn_type& conn);
void save_PlanetVector(db_conn_type& conn);
void save_PlanetElement(db_conn_type& conn);

// *****************************************************************************
int main(int argc, char* argv[])
{
    // *****************************************************************************
    // Process commandline arguments.
    // *****************************************************************************

    // Flags from commandline arguments
    bool run_MassiveBody = false;
    bool run_DetectionTime = false;
    bool run_Detection = false;
    bool run_DetectionCandidate = false;
    bool run_BodyVector = false;
    bool run_PlanetVector = false;
    bool run_PlanetElement = false;

    // Set up parser for named commandline arguments ("options")
    po::options_description desc("Find detections near asteroids");
    desc.add_options()
        ("help,h", "Produce help message")
        ("MassiveBody", po::bool_switch(&run_MassiveBody), 
            "Save contents of stored procedure KS.GetMassiveBody")
        ("DetectionTime", po::bool_switch(&run_DetectionTime), 
            "Save contents of stored procedure KS.GetDetectionTimes")
        ("Detection", po::bool_switch(&run_Detection), 
            "Save contents of stored procedure KS.GetDetections on full range")
        ("DetectionCandidate", po::bool_switch(&run_DetectionCandidate), 
            "Save contents of stored procedure KS.GetCandidateDetections on full range")
        ("BodyVector", po::bool_switch(&run_BodyVector), 
            "Save contents of stored procedure KS.GetVectors_Sun and KS.GetVectors_Earth on full range")
        ("PlanetVector", po::bool_switch(&run_PlanetVector), 
            "Save contents of stored procedure KS.GetVectors_Planets on full range")
        ("PlanetElement", po::bool_switch(&run_PlanetElement), 
            "Save contents of stored procedure KS.GetElements_Planets on full range")
    ;
    po::variables_map vm;

    // Parse commandline arguments including positional arguments
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    // *****************************************************************************
    // Program body
    // *****************************************************************************

    // Establish DB connection
    db_conn_type conn = get_db_conn();

    // Save MassiveBody if requested
    if (run_MassiveBody) {save_MassiveBody(conn);}

    // Save DetectionTime if requested
    if (run_DetectionTime) {save_DetectionTime(conn);}

    // Save Detection if requested
    if (run_Detection) {save_Detection(conn);}

    // Save DetectionCandidate if requested
    if (run_DetectionCandidate) {save_DetectionCandidate(conn);}

    // Save BodyVector if requested
    if (run_BodyVector) {save_BodyVector(conn);}

    // Save PlanetVector if requested
    if (run_PlanetVector) {save_PlanetVector(conn);}

    // Save PlanetElement if requested
    if (run_PlanetElement) {save_PlanetElement(conn);}

    // Close DB connection
    conn->close();

    // Normal program exit
    return 0;
}

// *****************************************************************************
// Save a file to disk with cached data from database.
// *****************************************************************************
// *****************************************************************************
void save_MassiveBody(db_conn_type& conn)
{
    // Instantiate a MassiveBodyTable
    MassiveBodyTable mbt(false);
    // Load it from DB
    mbt.load(conn);
    // Save it
    mbt.save();
    print("Saved MassiveBodyTable to disk.\n");
}

// *****************************************************************************
void save_DetectionTime(db_conn_type& conn)
{
    // Load the detection time table from database
    DetectionTimeTable dtt = DetectionTimeTable(conn);
    // Save detection time table to disk
    dtt.save();
    print("Saved DetectionTime to disk.\n");
}

// *****************************************************************************
void save_DetectionCandidate(db_conn_type& conn)
{
    // Load the detection candidate table from database
    bool progbar = true;
    DetectionCandidateTable dt = DetectionCandidateTable(conn, progbar);
    // Save detection candidate table to disk
    dt.save();
}

// *****************************************************************************
void save_Detection(db_conn_type& conn)
{
    // Load the detection candidate table from database
    bool progbar = true;
    DetectionTable dt = DetectionTable(conn, progbar);
    // Save detection table to disk
    dt.save();
}

// *****************************************************************************
void save_BodyVector(db_conn_type& conn)
{
    // Suported bodies: Sun, Earth
    constexpr std::array<SolarSystemBody_bv, 2> body_enums = 
        {SolarSystemBody_bv::sun, SolarSystemBody_bv::earth};
    
    // Iterate over supported bodies
    for (SolarSystemBody_bv body: body_enums)
    {
        // The body name
        string body_name = get_body_name(body);
        // Initialize BodyVector for this body using database
        print("Loading BodyVector for {:s} from DB...\n", body_name);
        BodyVector bv(body, conn);
        // Save to disk
        bv.save();
        print("Saved BodyVector to disk for {:s}.\n", body_name);
    }
}

// *****************************************************************************
void save_PlanetVector(db_conn_type& conn)
{
    // Build PlanetVector
    print("Building PlanetVector...\n");
    PlanetVector pv(conn);
    // Status update
    double mjd0 = pv.mjd0;
    double mjd1 = pv.mjd1;
    int dt_min = (mjd1-mjd0)*mpd/(pv.N_t-1);
    print("\nBuilt PlanetVector object from mjd0 {:8.4f} to mjd1 {:8.4f} with time step {:d} minutes.\n", 
            mjd0, mjd1, dt_min);
    // Save to disk
    pv.save();
    print("Saved PlanetVector to disk\n");
}

// *****************************************************************************
void save_PlanetElement(db_conn_type& conn)
{
    // Build PlanetElement
    print("Building PlanetElement...\n");
    PlanetElement pe(conn);
    // Status update
    double mjd0 = pe.get_mjd()[0];
    double mjd1 = pe.get_mjd()[pe.N_t-1];
    int dt_min = (mjd1-mjd0)*mpd/(pe.N_t-1);
    print("\nBuilt PlanetElement object from mjd0 {:8.4f} to mjd1 {:8.4f} with time step {:d} minutes.\n", 
            mjd0, mjd1, dt_min);
    // Save to disk
    pe.save();
    print("Saved PlanetElement to disk\n");
}
