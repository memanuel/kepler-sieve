"""
Harvard IACS Masters Thesis
Assemble orbital elements for known asteroids Pandas DataFrame.
Main interface used by consumers is load_ast_elt()

Michael S. Emanuel
Fri Aug 23 16:13:28 2019
"""

# Library imports
import numpy as np
import pandas as pd
import rebound

# Utility
from datetime import datetime

# Types
from typing import List, Tuple, Dict, Optional

# ********************************************************************************************************************* 
def load_data_numbered() -> pd.DataFrame:
    """Load the asteroid data for numbered asteroids into a Pandas DataFrame"""
    # The source for this file is at https://ssd.jpl.nasa.gov/?sb_elem
    fname: str = '../data/jpl/orbital_elements/asteroid_numbered.txt'

    # The field names in the JPL file and their column positions
    names: List[str] = ['Num', 'Name', 'Epoch', 'a', 'e', 'i', 'w', 'Node', 'M', 'H', 'G', 'Ref']
    colspec_tbl: Dict[str, Tuple[int, int]] = {
        'Num': (0,6), 
        'Name': (7, 25), 
        'Epoch': (25, 30), 
        'a': (31, 41), 
        'e': (42, 52), 
        'i': (54, 62), 
        'w': (63, 72),
        'Node': (73, 82),
        'M': (83, 94),
        'H': (95, 100),
        'G': (101, 105),
        'Ref': (106, 113),
    }
    
    # Other arguments for Pandas file import
    colspecs: List[Tuple[int, int]] = [colspec_tbl[nm] for nm in names]
    header: int = 0
    skiprows: List[int] = [1]
    dtype: Dict[str, int] = {
        'Num': int,
        'Name': str,
        'Epoch': float,
        'a': float,
        'e': float,
        'i': float,
        'w': float,
        'Node': float,
        'M': float,
        'H': float,
        'G': float,
        'Ref': str,
    }

    # Read the DataFrame
    df: pd.DataFrame = pd.read_fwf(fname, colspecs=colspecs, header=header, names=names, skiprows=skiprows, dtype=dtype)
    # Set the asteroid number field to be the index
    df.set_index(keys=['Num'], drop=False, inplace=True)
    return df

# ********************************************************************************************************************* 
def load_data_unnumbered() -> pd.DataFrame:
    """Load the asteroid data for unnumbered asteroids into a Pandas DataFrame"""

    fname: str = '../data/jpl/orbital_elements/asteroid_unnumbered.txt'

    # The field names in the JPL file and their column positions
    names: List[str] = ['Name', 'Epoch', 'a', 'e', 'i', 'w', 'Node', 'M', 'H', 'G', 'Ref']
    colspec_tbl: Dict[str, Tuple[int, int]] = {
        'Name': (0, 12), 
        'Epoch': (12, 17), 
        'a': (18, 28), 
        'e': (29, 39), 
        'i': (40, 50), 
        'w': (51, 60),
        'Node': (61, 70),
        'M': (71, 82),
        'H': (83, 88),
        'G': (89, 93),
        'Ref': (94, 104),
    }

    # Other arguments for Pandas file import
    colspecs: List[Tuple[int, int]] = [colspec_tbl[nm] for nm in names]
    header: int = 0
    skiprows: List[int] = [1]
    dtype: Dict[str, int] = {
        'Name': str,
        'Epoch': float,
        'a': float,
        'e': float,
        'i': float,
        'w': float,
        'Node': float,
        'M': float,
        'H': float,
        'G': float,
        'Ref': str,
    }

    # Read the DataFrame
    df: pd.DataFrame = pd.read_fwf(fname, colspecs=colspecs, header=header, names=names, skiprows=skiprows, dtype=dtype)
    # Populate the asteroid_num field and add it to DataFrame
    ast_num_offset = 1000001
    ast_num = df.index.values + ast_num_offset
    df.insert(loc=0, column='Num', value=ast_num)
    # Set the asteroid number field to be the index
    df.set_index(keys=['Num'], drop=False, inplace=True)
    return df

# ********************************************************************************************************************* 
def load_data_impl() -> pd.DataFrame:
    """Load the combined asteroid data into a Pandas DataFrame"""
    # Load main data file
    df1 = load_data_numbered()
    # Load auxiliary data file
    df2 = load_data_unnumbered()
    # Return combined DataFrame
    df = pd.concat([df1, df2])
    return df

# ********************************************************************************************************************* 
def convert_data(df_in: pd.DataFrame, epoch: Optional[float]=None) -> pd.DataFrame:
    """
    Convert data from the JPL format to be friendly to rebound integrator and matching selected epoch
    INPUTS:
        df_in: DataFrame with orbital elements in JPL format
        epoch: Optional epoch as an MJD; used to filter for only matching epochs.
               Defaults to the mode, e.g. 58600.0
    """
    # Apply the default value of epoch if it was not input
    if epoch is None:
        epoch = pd.Series.mode(df_in.Epoch)[0]

    # Create a mask with only the matching rows
    mask = (df_in.Epoch == epoch)
    
    # Initialize Dataframe with asteroid numbers
    df = pd.DataFrame(data=df_in.Num[mask])

    # Add fields one at a time
    df['Name'] = df_in.Name[mask]
    df['epoch'] = df_in.Epoch[mask]
    df['a'] = df_in.a[mask]
    df['e'] = df_in.e[mask]
    df['inc'] = np.radians(df_in.i[mask])
    df['Omega'] = np.radians(df_in.Node[mask])
    df['omega'] = np.radians(df_in.w[mask])
    df['M'] = np.radians(df_in.M[mask])
    df['H'] = df_in.H[mask]
    df['G'] = df_in.G[mask]
    df['Ref'] = df_in.Ref[mask]
    
    # Set the asteroid number field to be the index
    df.set_index(keys=['Num'], drop=False, inplace=True)

    # Return the newly assembled DataFrame
    return df

# ********************************************************************************************************************* 
def load_ast_elt() -> pd.DataFrame:
    """Load the asteroid orbital elements data into a Pandas Dataframe"""
    # The name for the saved DataFrame
    fname: str = '../data/jpl/orb_elements_asteroid.h5'
    
    # Try to load from disk if available
    ast_elt: pd.DataFrame
    try:
        ast_elt = pd.read_hdf(fname, key='ast_elt')
    except:
        # Load data from JPL asteroids file
        df_in = load_data_impl()
        # Convert data to rebound format
        ast_elt = convert_data(df_in=df_in)
        # Add the calculated orbital elements
        ast_elt = ast_data_add_calc_elements(ast_elt)
        # Add the row number field
        ast_elt['row_num'] = np.arange(ast_elt.shape[0], dtype=np.int32)
        # Save it to h5
        ast_elt.to_hdf(fname, key='ast_elt', mode='w')
    
    return ast_elt
