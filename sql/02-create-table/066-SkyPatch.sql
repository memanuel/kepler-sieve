-- ************************************************************************************************
-- SkyPatchGridCell describes the 4N^2 square grid cells on one major face
CREATE OR REPLACE TABLE KS.SkyPatchGrid(
    i SMALLINT NOT NULL
        COMMENT "Counter for the first Cartesian coordinate on this face",
    j SMALLINT NOT NULL
        COMMENT "Counter for the second Cartesian coordinate on this face",
    k INT NOT NULL
    	COMMENT "Counter for the entire grid",
    -- Center point
    u DOUBLE NOT NULL COMMENT "Midpoint of u (coordinate corresponding to i)",
    v DOUBLE NOT NULL COMMENT "Midpoint of v (coordinate corresponding to j)",
    w DOUBLE NOT NULL COMMENT "Midpoint of w (third coordinate, not i or j)",
	-- Lower left corner (u00, v00, w00)
    u00 DOUBLE NOT NULL COMMENT "Lower left corner, u",
    v00 DOUBLE NOT NULL COMMENT "Lower left corner, v",
    w00 DOUBLE NOT NULL COMMENT "Lower left corner, w",
	-- Upper left corner (u01, v01, w01)
    u01 DOUBLE NOT NULL COMMENT "Upper left corner, u",
    v01 DOUBLE NOT NULL COMMENT "Upper left corner, v",
    w01 DOUBLE NOT NULL COMMENT "Upper left corner, w",
	-- Lower right corner (u10, v10, w10)
    u10 DOUBLE NOT NULL COMMENT "Lower right corner, u",
    v10 DOUBLE NOT NULL COMMENT "Lower right corner, v",
    w10 DOUBLE NOT NULL COMMENT "Lower right corner, w",
	-- Upper right corner (u11, v11, w11)
    u11 DOUBLE NOT NULL COMMENT "Upper right corner, u",
    v11 DOUBLE NOT NULL COMMENT "Upper right corner, v",
    w11 DOUBLE NOT NULL COMMENT "Upper right corner, w",
    -- Unique key
    PRIMARY KEY (i,j)
    	COMMENT "The pair (i,j) uniquely determines one grid cell of a major face",
    UNIQUE KEY UNQ_SkyPatchGrid_k (k)
)
ENGINE='Aria' TRANSACTIONAL=0
COMMENT = "SkyPatchGrid describes the 4N^2 square grid cells on one major face";

-- ************************************************************************************************
-- The whole SkyPatch table; six major faces
CREATE OR REPLACE TABLE KS.SkyPatch(
	SkyPatchID INT NOT NULL PRIMARY KEY
        COMMENT "Integer ID for this patch of sky; assigned by formula (2N)^2*(f-1) + (2N)i + j",
    CubeFaceID TINYINT NOT NULL
        COMMENT "The major face of the cube on which the 2Nx2N grid is inscribed; foreign key to CubeFace",
    i SMALLINT NOT NULL
        COMMENT "Counter for the first Cartesian coordinate on this CubeFace, which is named alpha on the CubeFace table",
    j SMALLINT NOT NULL
        COMMENT "Counter for the second Cartesian coordinate on this CubeFace, which is named beta on the CubeFace table",
    -- Center point (x, y, z)
    x DOUBLE NOT NULL COMMENT "SkyPatch midpoint: x",
    y DOUBLE NOT NULL COMMENT "SkyPatch midpoint: y",
    z DOUBLE NOT NULL COMMENT "SkyPatch midpoint: z",
    -- Coordinates of lower left grid cell (x00, y00, z00)
    x00 DOUBLE NOT NULL COMMENT "SkyPatch lower left corner (00): x",
    y00 DOUBLE NOT NULL COMMENT "SkyPatch lower left corner (00): y",
    z00 DOUBLE NOT NULL COMMENT "SkyPatch lower left corner (00): z",
    -- Coordinates of upper left grid cell (x01, y01, z01)
    x01 DOUBLE NOT NULL COMMENT "SkyPatch upper left corner (01): x",
    y01 DOUBLE NOT NULL COMMENT "SkyPatch upper left corner (01): y",
    z01 DOUBLE NOT NULL COMMENT "SkyPatch upper left corner (01): z",
    -- Coordinates of lower right grid cell (x10, y10, z10)
    x10 DOUBLE NOT NULL COMMENT "SkyPatch lower right corner (10): x",
    y10 DOUBLE NOT NULL COMMENT "SkyPatch lower right corner (10): y",
    z10 DOUBLE NOT NULL COMMENT "SkyPatch lower right corner (10): z",
    -- Coordinates of lower upper right cell (x11, y11, z11)
    x11 DOUBLE NOT NULL COMMENT "SkyPatch upper right corner (11): x",
    y11 DOUBLE NOT NULL COMMENT "SkyPatch upper right corner (11): y",
    z11 DOUBLE NOT NULL COMMENT "SkyPatch upper right corner (11): z",
    -- Unique key
    UNIQUE KEY UNQ_SkyPatch_CubeFaceID_i_j(CubeFaceID, i, j)
    	COMMENT "The pair (i,j) determines one small patch on a major face; the trio (f,i,j) is unique",
    UNIQUE KEY UNQ_SkyPatch_CubeFaceID_j_i(CubeFaceID, j, i)
    	COMMENT "also support searching when f and j are known but not i",
    -- Foreign key
    CONSTRAINT FK_SkyPatch_CubeFaceID FOREIGN KEY (CubeFaceID) REFERENCES CubeFace(CubeFaceID)
)
ENGINE='Aria' TRANSACTIONAL=0
COMMENT "Collection of discrete patches of the sky corresponding to a cube in which the unit sphere is inscribed.";

-- ************************************************************************************************
-- Neighbor distance on the SkyPatchGrid table
CREATE OR REPLACE TABLE KS.SkyPatchGridDistance(
	-- Keys for the two grid cells being connected
	i1 INT NOT NULL
        COMMENT "i coordinate of the first SkyPatchGrid cell",
	j1 INT NOT NULL
        COMMENT "j coordinate of the first SkyPatchGrid cell",
	i2 INT NOT NULL
        COMMENT "i coordinate of the second SkyPatchGrid cell",
	j2 INT NOT NULL
        COMMENT "j coordinate of the second SkyPatchGrid cell",
    -- Distances: midpoint and minimum
	dr_mid DOUBLE NOT NULL
		COMMENT "Distance from center point to center point",
	dr_min DOUBLE NOT NULL
		COMMENT "Estimated minimum distance based on getting closer by up to half the width on each end",
	-- Keys
	PRIMARY KEY (i1, j1, i2, j2),
	INDEX IDX_i1_j1_dr_min (i1, j1, dr_min)
)
ENGINE='Aria' TRANSACTIONAL=0
COMMENT "Distance bewteen two SkyPatchGrid cells; only cataloged for neighbors that are reasonably close.";


-- ************************************************************************************************
-- Neighbor distance on the SkyPatch table
CREATE OR REPLACE TABLE KS.SkyPatchDistance(
	SkyPatchID_1 INT NOT NULL
        COMMENT "The first SkyPatchID",
	SkyPatchID_2 INT NOT NULL
        COMMENT "The second SkyPatchID",
	dr_mid DOUBLE NOT NULL
		COMMENT "Distance from center point to center point",
	dr_min DOUBLE NOT NULL
		COMMENT "Estimated minimum distance based on getting closer by up to half the width on each end",
	IsCrossFace BOOL,
	-- Keys
	PRIMARY KEY (SkyPatchID_1, SkyPatchID_2),
	INDEX (SkyPatchID_1, dr_min)
		COMMENT "Support searching for all SkyPatch cells within a distance of SkyPatchID_1",
	INDEX IDX_SkyPatchDistance_IsCrossFace (IsCrossFace)	
)
COMMENT "Distance bewteen two SkyPatch cells; only cataloged for neighbors that are reasonably close.";
