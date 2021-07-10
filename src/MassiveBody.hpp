/** @file   MassiveBody.hpp
 *  @brief  Class to load all massive bodies in the DE-435 integration.
 *  See DB table KS.MassiveBody and stored procedure KS.GetMassiveBodies.
 *  
 *  @author Michael S. Emanuel
 *  @date   2021-07-08
 */

// *****************************************************************************
#pragma once

// *****************************************************************************
// Included libraries
#include <cstdint>
#include <string>
    using std::string;
#include <fstream>
    using std::ifstream;
    using std::ofstream;
#include <vector>
    using std::vector;
#include <map>
    using std::map;


// *****************************************************************************
#include "db_utils.hpp"
    using ks::db_conn_type;
    using ks::get_db_conn;
    using ks::sp_run;

// *****************************************************************************
namespace ks {

// *****************************************************************************
struct MassiveBody
{
    // The integer ID of this body; FK to KS.Body
    int32_t body_id;
    // The mass of this body in solar units
    double M;
    // The graviational field strength in units of au^3 / day^2
    double GM;
};

// *****************************************************************************
class MassiveBodyTable
{
public:
    /// Default constructor initializes MassiveBodyTable from cached file on disk
    MassiveBodyTable();
    /// Explicit constructor initializes table with option to load or leave it empty
    MassiveBodyTable(bool load_from_disk);
    /// Default destructor is fine
    ~MassiveBodyTable();

    /// Get the mass of a body in Msun given its body_id
    const double get_M(int32_t body_id) const;
    /// Get the gravitational field strength of a body in AU^2 / day^3 given its body_id
    const double get_GM(int32_t body_id) const;
    /// Get a MassiveBody structure for a body given its body_id
    const MassiveBody operator[](int32_t body_id) const;

    /// Get size of this table (number of bodies)
    const int size() const;
    /// Get array of body_id
    const vector<int32_t>& get_body_id() const;
    /// Get array of masses
    const vector<double>& get_M() const;
    /// Get array of graviational field strength
    const vector<double>& get_GM() const;

    /// Load with data from database
    void load(db_conn_type& conn);
    /// Save to disk
    void save();
    /// Load with data from disk
    void load();

private:
    /// Vector of body_ids
    vector<int32_t> body_id;
    /// Vector of body masses; position in this array is row_id
    vector<double> M;
    /// Vector of body gravitational field strength; position in this array is row_id
    vector<double> GM;
    /// Map with key = body_id, value = row_id
    map<int32_t, int32_t> row_map;
    /// Get the row_id of a body given its body_id
    const int32_t get_row_id(int32_t body_id) const;
    /// The current row_id
    const int32_t current_row_id() const;
    /// Number of rows in data file
    const int file_length() const;
};

// *****************************************************************************
} // namespace ks
