"""
Convert text file with masses used in DE 435 integration to a CSV suitable for database import.
Michael S. Emanuel
2020-08-04
"""

# ************************************************************************************************    
def convert_line(l: str):
    """Convert one line to CSV format"""
    tokens = l.split()
    param_name = tokens[0]
    mass_str = tokens[1].replace('D', 'E')
    return f'{param_name},{mass_str}\n'
    
# ************************************************************************************************    
def main():
    """Write out masses.csv"""
    # Names of input and output files
    fname_in = 'masses.txt'
    fname_out = 'masses.csv'

    # Read the input file and read the lines
    with open(fname_in) as fh_in:
        lines = fh_in.readlines()

    # The first line has headers; convert the remaining lines
    with open(fname_out, 'w') as fh_out :
        # Write out the headers for the processed file
        fh_out.write('ParameterName,GM_AU3_per_Day2\n')
        # Write out the data rows
        for line_in in lines[1:]:
            line_out = convert_line(line_in)
            fh_out.write(line_out)
    
# ************************************************************************************************    
if __name__ == '__main__':\
    main()
