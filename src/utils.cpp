// *****************************************************************************
// Included files
#include <iostream>
#include <boost/format.hpp>
#include "utils.h"

// *****************************************************************************
// Names used
using std::cout;
using boost::format;

// *****************************************************************************
// Put functions into the namespace ks (for Kepler Sieve)
namespace ks {

// *****************************************************************************
void print_stars()
{
    cout << "********************************************************************************\n";
}

// *****************************************************************************
void print_newline()
{
    cout << "\n";
}

// *****************************************************************************
const string test_message(bool is_ok)
{
    if (is_ok)
        return string("PASS");
    else
        return string("FAIL");
}

// *****************************************************************************
void report_test(const string test_name, bool is_ok)
{
    cout << format("%s:\n") % test_name;
    cout << format("**** %s ****\n") % test_message(is_ok);
}

// *****************************************************************************
}; // namespace
