"""
Astronomy: Transformation between a direction in the sky and a SkyPatch.

Michael S. Emanuel
2021-03-29
"""

# Core
import numpy as np
import pandas as pd

# ********************************************************************************************************************
# Load one copy of the Cube DataFrame into memory
# cf = sp2df('KS.GetCubeFace')
cf = pd.DataFrame(
    {'f': np.arange(6),
    'j1': np.array([1,3,2,2,3,1], dtype=np.int32),
    'j2': np.array([2,1,3,3,1,2], dtype=np.int32),
    'i' : np.array([3,2,1,1,2,3], dtype=np.int32),
    'ci': np.array([1,1,1,-1,-1,-1])    
    })

# ********************************************************************************************************************
def dir2SkyPatchID(dir: np.array, N: int=1024):
    """
    Compute a vector of SkyPatchIDs from a 3 vectors of direction components
    INPUTS:
        dir: An array [ux, uy, uz] of directions on the unit sphere in the ecliptic frame; shape N_row x 3
        N: Grid size for SkyPatch
    OUTPUTS:
        SkyPatchID: An array of N_row SkyPatchIDs.
    """

    # Unpack components of directions
    x = dir[:,0]
    y = dir[:,1]
    z = dir[:,2]

    # The search array for k is ordered z, y, x, -x, -y, -z
    sa = np.stack([z, y, x, -x, -y, -z]).T

    # The index of of the largest element; zero based.  CubeFaceID = f+1
    f = np.int32(np.argmax(sa, axis=1))

    # Row indexer for 2D array indexing
    N_row: int = sa.shape[0]
    row_idx = np.arange(N_row)

    # Indices for axes corresponding to u and v and w
    idx_j1 = cf.j1[f].values-1
    idx_j2 = cf.j2[f].values-1
    idx_i = cf.i[f].values-1

    # Calculate three components (u, v, w) on unit sphere
    u = dir[row_idx, idx_j1]
    v = dir[row_idx, idx_j2]
    w = dir[row_idx, idx_i]

    # The factor t required to dilate (u, v, w) to the nearest cube face
    t = 1.0 / np.abs(w)

    # The projection of the direction onto the cube
    a = t*u
    b = t*v
    # c = t*w

    # Compute the integer entries (i, j)
    i = np.int32(np.floor(N * (1.0 + a)))
    j = np.int32(np.floor(N * (1.0 + b)))

    # Caclulate SkyPatchID
    M: int = np.int32(2*N)
    M2: int = np.int32(M*M)
    SkyPatchID = f*M2 + i*M + j

    return SkyPatchID

# ********************************************************************************************************************
def SkyPatchID2dir(SkyPatchID: np.array, N: int):
    """
    Compute a vector of SkyPatchIDs from a 3 vectors of direction components
    INPUTS:
        SkyPatchID: An array of N_row SkyPatchIDs.
        N: Grid size for SkyPatch
    OUTPUTS:
        dir_mid: An array [ux, uy, uz] of directions on the unit sphere in the ecliptic frame; shape N_row x 3
                 These directions are taken at the midpoint of the SkyPatch.
    """

    # Number of rows
    N_row: int = SkyPatchID.shape[0]

    # Calculate M and M2 from N
    M: int = np.int32(2*N)
    M2: int = np.int32(M*M)

    # Extract face number f and grid coordinates i, j
    f, x = np.divmod(SkyPatchID, M2)
    i, j = np.divmod(x, M)

    # Compute a, b, c from f, i, j
    a = -1.0 + (i + 0.5) / N
    b = -1.0 + (j + 0.5) / N
    c = cf.ci.loc[f].values

    # The distance r at the cube face
    r = np.sqrt(np.square(a) + np.square(b) + 1.0)

    # Compute (u, v, w) from (a, b, c)
    u = a / r
    v = b / r
    w = c / r

    # Index positions of (u, v, w) in (x, y, z)
    idx = cf.loc[f, ['j1', 'j2', 'i']]

    # Compute direction of the SkyPatch midpoint
    # dir_mid = (x, y, z) from (u, v, w)
    dir_mid = np.zeros((N_row, 3))
    dir_mid[:,0] = (idx.j1==1)*u + (idx.j2==1)*v + (idx.i==1)*w
    dir_mid[:,1] = (idx.j1==2)*u + (idx.j2==2)*v + (idx.i==2)*w
    dir_mid[:,2] = (idx.j1==3)*u + (idx.j2==3)*v + (idx.i==3)*w
    return dir_mid
