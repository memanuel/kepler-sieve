// *****************************************************************************
// Library dependencies
#include <fmt/format.h>
    using fmt::print;

// Local dependencies
#include "db_utils.hpp"
    using ks::db_conn_type;
    using ks::get_db_conn;
    using ks::sp_run;

#include "DetectionCandidate.hpp"
    using ks::DetectionCandidate;
    using ks::DetectionCandidateTable;

// *****************************************************************************
int main()
{
    // Establish DB connection
    db_conn_type conn = get_db_conn();

    // Load the detection candidate table
    int d0 = 0;
    int d1 = 1000000;
    bool progbar = true;
    // DetectionCandidateTable dt = DetectionCandidateTable(conn, d0, d1, progbar);

    // Save detection table to disk
    // dt.serialize();

    // Load detection table from disk
    DetectionCandidateTable dt;
    dt.load();

}