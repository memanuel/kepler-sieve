/** @file Detection.cpp
 *  @brief Implmentation of Detection class.
 *
 *  @author Michael S. Emanuel
 *  @date 2021-06-28
 */

// *****************************************************************************
// Included files
#include "Detection.hpp"

// *****************************************************************************
// Local names used
using ks::Detection;
using ks::DetectionTable;

// Set batch size to one million
constexpr int batch_size = 1000000;

// Location of file with serialized data
const string file_name = "data/cache/DetectionTable.bin";

// *****************************************************************************
// Class DetectionTable
// *****************************************************************************

// *****************************************************************************
DetectionTable::DetectionTable(): 
    d0(0),
    d1(0),
    dt(vector<Detection>(0)),
    dtsp(vector<vector<int32_t>>(N_sp))
{
    // Load data from disk
    load();
}

// *****************************************************************************
// This constructor is used for loading data in chunks from the database in conjunction with load.
DetectionTable::DetectionTable(int d0, int d1): 
    d0(d0),
    d1(d1),
    // Initialize dt to a vector with sz entries, one for each possible detection in the interval
    dt(vector<Detection>(d1-d0)),
    // Initialize dtsp to a vector with N_sp entries, one for each SkyPatch (whether occupied or not)
    dtsp(vector<vector<int32_t>>(N_sp))
    {}

// *****************************************************************************
DetectionTable::DetectionTable(db_conn_type &conn, bool progbar): 
    // Delegate to range constructor, using DB to compute d1
    DetectionTable(0, sp_run_int(conn, "KS.GetMaxDetectionID")) 
{
    // Load data from database
    load(conn, progbar);
}

// *****************************************************************************
void DetectionTable::load(db_conn_type& conn, bool progbar)
{
    // Write 0 into DetectionID field so we can later ignore any holes (e.g. DetectionID=0)
    int sz = d1-d0;
    for (int i=0; i<sz; i++)
    {
        dt[i].detection_id=0;
    }

    // Calculate number of batches
    int batch_count = std::max((sz+1)/batch_size, 1);
    // Status update
    if (progbar) 
    {
        print("Processing {:d} detections from {:d} to {:d} in {:d} batches of {:d}...\n",
                sz, d0, d1, batch_count, batch_size);
    }

    // Timer for processing detections from DB
    Timer t;
    t.tick();

    // Iterate over the batches
    for (int i0=d0; i0<d1; i0+=batch_size)
    {
        // Upper limit for this batch
        int i1 = std::min(i0+batch_size, d1);
        // Process SQL data in this batch
        process_rows(conn, i0, i1);
        // Progress bar
        if (progbar) 
        {
            print("."); 
            flush_console();
        }
    }
    if (progbar) 
    {
        print("\nLoaded detection candidates table with {:d} entries.\n", dt.size());
        t.tock_msg();
    }
}

// *****************************************************************************
/// Helper function: Process a batch of rows, writing data from stored procedure output to vector of detections.
void DetectionTable::process_rows(db_conn_type& conn, int i0, int i1)
{
    // Run the stored procedure to get detections including the observatory position
    // string sp_name = "KS.GetDetectionsObs";
    string sp_name = "KS.GetDetections";
    vector<string> params = {to_string(i0), to_string(i1)};
    ResultSet* rs = sp_run(conn, sp_name, params);

    // Loop through resultset
    while (rs->next()) 
    {
        // Unpack the fields in the resultset
        // 3 key fields
        int32_t detection_id = rs->getInt("DetectionID");
        int32_t sky_patch_id = rs->getInt("SkyPatchID");
        int32_t time_id = rs->getInt("TimeID");
        int32_t detection_time_id = rs->getInt("DetectionTimeID");
        // The detection time
        double mjd = rs->getDouble("mjd");
        // 3 direction components
        double ux = rs->getDouble("ux");
        double uy = rs->getDouble("uy");
        double uz = rs->getDouble("uz");
        // Magnitude
        double mag = rs->getDouble("mag");

        // Initialize the Detection at this location in dt
        int idx = detection_id - d0;
        dt[idx] = 
        {
            .detection_id = detection_id,
            .sky_patch_id = sky_patch_id,
            .time_id = time_id,
            .detection_time_id = detection_time_id,
            .mjd = mjd,
            .ux = ux,
            .uy = uy,
            .uz = uz,
            .mag = mag,
        };
        // Write this DetectionID to the vector keyed by this SkyPatchID
        (dtsp[sky_patch_id]).push_back(detection_id);
    }   // while rs
    // Close the resultset and free memory
    rs->close();
    delete rs;
}

// *****************************************************************************
///Default destructor is OK here
DetectionTable::~DetectionTable() {}

// *****************************************************************************
const int DetectionTable::size() const
{
    return (d1-d0);
}

// *****************************************************************************
Detection DetectionTable::operator[](int32_t id) const
{
    return dt[id];
}

// *****************************************************************************
vector<int32_t> DetectionTable::get_skypatch(int32_t spid) const
{
    return dtsp[spid];
} // function

// *****************************************************************************
void DetectionTable::save()
{
    // Open output filestream in binary output; truncate file contents
    std::ofstream fs;
    std::ios_base::openmode file_mode = (std::ios::out | std::ios::binary | std::ios::trunc);
    fs.open(file_name, file_mode);

    // Write the number of rows in binary as long int
    long sz = dt.size();
    fs.write((char*) &sz, sizeof(sz));

    // Write the rows to the file in binary
    for (Detection d : dt)
    {
        fs.write((char*) &d, sizeof(d));
    }

    // Status
    // print("Wrote data to {:s}\n", file_name);

    // Close output filestream
    fs.close();
}

// *****************************************************************************
void DetectionTable::load()
{
    // Open input filestream in binary mode
    std::ifstream fs;
    std::ios_base::openmode file_mode = (std::ios::in | std::ios::binary);
    fs.open(file_name, file_mode);

    // Read the number of rows in the file
    long sz=-1;
    fs.read( (char*) &sz, sizeof(sz));

    // Two modes of operation.
    // (1) auto mode: d1==0, load the whole file
    // (2) manual mode: d1>0 was previously set; load only this chunk    
    // Update paramater d1; if placeholder of 0, update it to sz; otherwise it is min(d1, sz)
    d1 = (d1>0) ? min(d1, (int) sz) : sz;

    // Allocate space for the detection table
    dt.resize(sz);

    // Read the rows from the file in binary
    Detection d;
    for (int i=0; i<sz; i++)
    {
        // Read this row into d
        fs.read( (char*) &d, sizeof(d));
        // Get the detection_id
        int32_t detection_id = d.detection_id;
        // If the detection_id not in the requested range [d0, d1), skip it
        if (detection_id < d0) {continue;}
        if (detection_id >= d1) {break;}
        // Save this element to dt
        dt[detection_id-d0] = d;
        // Write this DetectionID to the vector keyed by this SkyPatchID
        (dtsp[d.sky_patch_id]).push_back(detection_id);
    }

    // Close input filestream
    fs.close();
}
