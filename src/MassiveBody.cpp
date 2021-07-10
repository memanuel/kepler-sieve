/** @file MassiveBody.cpp
 *  @brief Implmentation of MassiveBody class.
 *  
 *  @author Michael S. Emanuel
 *  @date 2021-07-08
 */

// *****************************************************************************
// Included files
#include "MassiveBody.hpp"

// *****************************************************************************
// Local names used
using ks::MassiveBody;
using ks::MassiveBodyTable;

// *****************************************************************************
// Constants used in this module

// Location of file with serialized data
const string file_name = "data/cache/MassiveBody.bin";

// *****************************************************************************
// Constructor and destructor
// *****************************************************************************

// *****************************************************************************
// Explicit constructor
MassiveBodyTable::MassiveBodyTable(bool load_from_disk):
    body_id(vector<int32_t>(0)),
    M(vector<double>(0)),
    GM(vector<double>(0)),
    row_map(map<int32_t, int32_t>{})
{
    // Load data if requested
    if (load_from_disk) {load();}
}

// *****************************************************************************
// Default constructor delegates to explicit constructor with load = true
MassiveBodyTable::MassiveBodyTable():
    MassiveBodyTable(true)
    {};

// *****************************************************************************
// Default destructor
MassiveBodyTable::~MassiveBodyTable()
{}

// *****************************************************************************
// Access data about one body
// *****************************************************************************

// *****************************************************************************
const double MassiveBodyTable::get_M(int32_t body_id) const
{
    const int32_t row_id = get_row_id(body_id);
    return M[row_id];
}

// *****************************************************************************
const double MassiveBodyTable::get_GM(int32_t body_id) const
{
    const int32_t row_id = get_row_id(body_id);
    return GM[row_id];
}

// *****************************************************************************
const MassiveBody MassiveBodyTable::operator[](int32_t body_id) const
{
    // Get the mass and field strength
    double M = get_M(body_id);
    double GM = get_GM(body_id);
    // Wrap up the answer into a MassiveBody instance
    return MassiveBody {.body_id=body_id, .M= M, .GM=GM};
}

// *****************************************************************************
// Access arrays of data
// *****************************************************************************

// *****************************************************************************
const int MassiveBodyTable::size() const 
    {return body_id.size();}

// *****************************************************************************
const vector<int32_t>& MassiveBodyTable::get_body_id() const 
    {return body_id;}

// *****************************************************************************
const vector<double>& MassiveBodyTable::get_M() const 
    {return M;}

// *****************************************************************************
const vector<double>& MassiveBodyTable::get_GM() const 
    {return GM;}

// *****************************************************************************
// Implementation functions - work with row IDs
// *****************************************************************************

// *****************************************************************************
const int32_t MassiveBodyTable::get_row_id(int32_t body_id) const
    {return row_map.at(body_id);}

// *****************************************************************************
const int32_t MassiveBodyTable::current_row_id() const 
    {return body_id.size();}

// *****************************************************************************
// Load from database
// *****************************************************************************

// *****************************************************************************
void MassiveBodyTable::load(db_conn_type& conn)
{
    // Run the stored procedure to get detections including the observatory position
    string sp_name = "KS.GetMassiveBodies";
    vector<string> params = {};
    ResultSet* rs = sp_run(conn, sp_name, params);

    // Loop through resultset
    while (rs->next())
    {
        // Get the body_id
        int32_t body_id = rs->getInt("BodyID");
        // Get the mass
        double M = rs->getDouble("M");
        // Get the grivational field strength
        double GM = rs->getDouble("GM");

        // Get the current row_id
        int row_id = current_row_id();
        // Add this row_id to the map keyed by body_id
        row_map[body_id] = row_id;

        // Write the components into the vectors
        this->body_id.push_back(body_id);
        this->M.push_back(M);
        this->GM.push_back(GM);
    }
    // Close the resultset and free memory
    rs->close();
    delete rs;
}

// *****************************************************************************
// Load from database
// *****************************************************************************

// *****************************************************************************
void MassiveBodyTable::save()
{
    // Open output filestream in binary output; truncate file contents
    std::ofstream fs;
    std::ios_base::openmode file_mode = (std::ios::out | std::ios::binary | std::ios::trunc);
    fs.open(file_name, file_mode);

    // Write the number of rows in binary as long int
    long sz = body_id.size();
    fs.write((char*) &sz, sizeof(sz));

    // Write the rows to the file in binary
    for (int32_t body_id : this->body_id)
    {
        MassiveBody mb = (*this)[body_id];
        fs.write((char*) &mb, sizeof(mb));
    }

    // Close output filestream
    fs.close();
}

// *****************************************************************************
const int MassiveBodyTable::file_length() const
{
    // Open input filestream in binary mode
    std::ifstream fs;
    std::ios_base::openmode file_mode = (std::ios::in | std::ios::binary);
    fs.open(file_name, file_mode);
    // Read the number of rows in the file
    long sz=-1;
    fs.read( (char*) &sz, sizeof(sz));
    // Close input filestream
    fs.close();
    // Return file as an int, not a long
    return (int) sz;
}

// *****************************************************************************
void MassiveBodyTable::load()
{
    // Open input filestream in binary mode
    std::ifstream fs;
    std::ios_base::openmode file_mode = (std::ios::in | std::ios::binary);
    fs.open(file_name, file_mode);

    // Read the number of rows in the file
    long sz=-1;
    fs.read( (char*) &sz, sizeof(sz));
    // Status
    // print("Opened file {:s} with {:d} rows of DetectionTime data.\n", file_name, sz);

    // Read the rows from the file in binary
    MassiveBody mb;
    for (int i=0; i<sz; i++)
    {
        // Read this row into mb
        fs.read( (char*) &mb, sizeof(mb));
        // Get the current row_id
        int row_id = current_row_id();
        // Add this row_id to the map keyed by body_id
        row_map[mb.body_id] = row_id;

        // Write the components into the vectors
        this->body_id.push_back(mb.body_id);
        this->M.push_back(mb.M);
        this->GM.push_back(mb.GM);
    }
    // Close input filestream
    fs.close();
}
