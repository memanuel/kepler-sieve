-- ************************************************************************************************
-- Set the variable N
-- SET @N = POW(2,9);
SET @N = 4;

-- The size of each grid is 2Nx2N with 4N^2 entries
SET @M = 4 * @N * @N;

-- Width for populating the SkyPatchGridNeighbor table
SET @grid_width = 5;

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
    	COMMENT "The pair (i,j) determines one small patch on a major face; the trio (f,i,j) is unique"
    -- Foreign key
    -- CONSTRAINT FK_SkyPatch_CubeFaceID FOREIGN KEY (CubeFaceID) REFERENCES CubeFace(CubeFaceID)
)
ENGINE='Aria' TRANSACTIONAL=0
COMMENT "Collection of discrete patches of the sky corresponding to a cube in which the unit sphere is inscribed.";

-- SELECT * FROM KS.SkyPatch;
INSERT INTO KS.SkyPatch
(SkyPatchID, CubeFaceID, i, j, x, y, z, x00, y00, z00, x01, y01, z01, x10, y10, z10, x11, y11, z11)
SELECT
	-- Integer IDs
	(cf.CubeFaceID-1)*@M + 2*gr.i*@N + j AS SkyPatchID,
	cf.CubeFaceID,
	-- Local grid coordinates (i, j) on the major face
	gr.i,
	gr.j,
	-- Coordinates of midpoint (x, y, z)
	IF(cf.alpha='X', gr.u, 0.0) + IF(cf.beta='X', gr.v, 0.0) + IF(cf.gamma='X', gr.w*cf.ci, 0.0) AS x,
	IF(cf.alpha='Y', gr.u, 0.0) + IF(cf.beta='Y', gr.v, 0.0) + IF(cf.gamma='Y', gr.w*cf.ci, 0.0) AS y,
	IF(cf.alpha='Z', gr.u, 0.0) + IF(cf.beta='Z', gr.v, 0.0) + IF(cf.gamma='Z', gr.w*cf.ci, 0.0) AS z,
	-- Coordinates of lower left corner (x00, y00, z00)
	IF(cf.alpha='X', gr.u00, 0.0) + IF(cf.beta='X', gr.v00, 0.0) + IF(cf.gamma='X', gr.w00*cf.ci, 0.0) AS x00,
	IF(cf.alpha='Y', gr.u00, 0.0) + IF(cf.beta='Y', gr.v00, 0.0) + IF(cf.gamma='Y', gr.w00*cf.ci, 0.0) AS y00,
	IF(cf.alpha='Z', gr.u00, 0.0) + IF(cf.beta='Z', gr.v00, 0.0) + IF(cf.gamma='Z', gr.w00*cf.ci, 0.0) AS z00,
	-- Coordinates of upper left corner (x01, y01, z01)
	IF(cf.alpha='X', gr.u01, 0.0) + IF(cf.beta='X', gr.v01, 0.0) + IF(cf.gamma='X', gr.w01*cf.ci, 0.0) AS x01,
	IF(cf.alpha='Y', gr.u01, 0.0) + IF(cf.beta='Y', gr.v01, 0.0) + IF(cf.gamma='Y', gr.w01*cf.ci, 0.0) AS y01,
	IF(cf.alpha='Z', gr.u01, 0.0) + IF(cf.beta='Z', gr.v01, 0.0) + IF(cf.gamma='Z', gr.w01*cf.ci, 0.0) AS z01,
	-- Coordinates of lower right corner (x10, y10, z10)
	IF(cf.alpha='X', gr.u10, 0.0) + IF(cf.beta='X', gr.v10, 0.0) + IF(cf.gamma='X', gr.w10*cf.ci, 0.0) AS x10,
	IF(cf.alpha='Y', gr.u10, 0.0) + IF(cf.beta='Y', gr.v10, 0.0) + IF(cf.gamma='Y', gr.w10*cf.ci, 0.0) AS y10,
	IF(cf.alpha='Z', gr.u10, 0.0) + IF(cf.beta='Z', gr.v10, 0.0) + IF(cf.gamma='Z', gr.w10*cf.ci, 0.0) AS z10,
	-- Coordinates of upper right corner (x11, y11, z11)
	IF(cf.alpha='X', gr.u11, 0.0) + IF(cf.beta='X', gr.v11, 0.0) + IF(cf.gamma='X', gr.w11*cf.ci, 0.0) AS x11,
	IF(cf.alpha='Y', gr.u11, 0.0) + IF(cf.beta='Y', gr.v11, 0.0) + IF(cf.gamma='Y', gr.w11*cf.ci, 0.0) AS y11,
	IF(cf.alpha='Z', gr.u11, 0.0) + IF(cf.beta='Z', gr.v11, 0.0) + IF(cf.gamma='Z', gr.w11*cf.ci, 0.0) AS z11
FROM
	KS.CubeFace AS cf,
	KS.SkyPatchGrid AS gr
ORDER BY SkyPatchID;

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
	PRIMARY KEY (i1, j1, i2, j2)
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
	PRIMARY KEY (SkyPatchID_1, SkyPatchID_2)
)
COMMENT "Distance bewteen two SkyPatch cells; only cataloged for neighbors that are reasonably close.";
