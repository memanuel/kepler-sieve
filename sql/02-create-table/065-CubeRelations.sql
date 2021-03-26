-- ************************************************************************************************
-- Relate Cube vertex to edge
CREATE OR REPLACE TABLE KS.CubeEdgeVertex(
    CubeEdgeID TINYINT NOT NULL,
    CubeVertexID TINYINT NOT NULL,
    PRIMARY KEY (CubeEdgeID, CubeVertexID)
)
ENGINE='Aria' TRANSACTIONAL=1
COMMENT = "Relate a cube edge to a cube vertex if the vertex is one end of the edge.";

INSERT INTO KS.CubeEdgeVertex
(CubeEdgeID, CubeVertexID)
SELECT
	e.CubeEdgeID,
	cv.CubeVertexID
FROM
	KS.CubeEdge AS e
	INNER JOIN KS.CubeVertex AS cv ON cv.CubeVertexID = e.CubeVertexID_1
UNION
SELECT
	e.CubeEdgeID,
	cv.CubeVertexID
FROM
	KS.CubeEdge AS e
	INNER JOIN KS.CubeVertex AS cv ON cv.CubeVertexID = e.CubeVertexID_2;

-- ************************************************************************************************
-- Relate Cube face to edge
CREATE OR REPLACE TABLE KS.CubeFaceEdge(
	CubeFaceID TINYINT NOT NULL,
    CubeEdgeID TINYINT NOT NULL,
    -- CubeVertexID_1 TINYINT NOT NULL,
    -- CubeVertexID_2 TINYINT NOT NULL,
    PRIMARY KEY (CubeFaceID, CubeEdgeID),
    UNIQUE KEY UNQ_CubeFaceEdge_EdgeID_FaceID (CubeEdgeID, CubeFaceID)
)
ENGINE='Aria' TRANSACTIONAL=1
COMMENT = "Relate a cube face to a cube edge if the edge is on the perimeter of the face.";

INSERT INTO KS.CubeFaceEdge
(CubeFaceID, CubeEdgeID)
WITH t1 AS (
SELECT
	f.CubeFaceID,
	f.CubeFaceCD,
	f.i,
	f.j1,
	f.j2,
	f.ci,
	s1._ AS s1,
	s2._ AS s2,
	IF(f.i=1, f.ci, 0) + IF(f.j1=1, s1._, 0) + IF(f.j2=1, s2._, 0) AS x,
	IF(f.i=2, f.ci, 0) + IF(f.j1=2, s1._, 0) + IF(f.j2=2, s2._, 0) AS y,
	IF(f.i=3, f.ci, 0) + IF(f.j1=3, s1._, 0) + IF(f.j2=3, s2._, 0) AS z
FROM
	KS.CubeFace AS f
	INNER JOIN KS.CounterSigned AS s1 ON s1._ IN (-1, 1) 
	INNER JOIN KS.CounterSigned AS s2 ON s2._ IN (-1, 1)
),
t2a AS (
SELECT
	t1.CubeFaceID,
	v.CubeVertexID
FROM
	t1
	INNER JOIN KS.CubeVertex AS v ON v.x=t1.x AND v.y=t1.y AND v.z=t1.z
),
t2b AS (
SELECT
	t1.CubeFaceID,
	v.CubeVertexID
FROM
	t1
	INNER JOIN KS.CubeVertex AS v ON v.x=t1.x AND v.y=t1.y AND v.z=t1.z
)
SELECT
	f.CubeFaceID,
	e.CubeEdgeID
	-- e.CubeVertexID_1,
	-- e.CubeVertexID_2
FROM
	KS.CubeFace AS f
	INNER JOIN t2a AS v1 ON v1.CubeFaceID=f.CubeFaceID
	INNER JOIN t2b AS v2 ON v2.CubeFaceID=f.CubeFaceID
	INNER JOIN KS.CubeEdge AS e ON
		e.CubeVertexID_1 = v1.CubeVertexID AND
		e.CubeVertexID_2 = v2.CubeVertexID
WHERE
	v1.CubeVertexID <> v2.CubeVertexID;

-- ************************************************************************************************
-- Relate Cube face to its four neighbors by orientation
CREATE OR REPLACE TABLE KS.CubeFaceNeighbor(
	CubeFaceID TINYINT NOT NULL PRIMARY KEY
		COMMENT "The face whose neighbors are being enumerated",
	CubeFaceID_i0 TINYINT NOT NULL
		COMMENT "Neighbor on grid when i=0; this is where axis j1 is equal to -1",
	CubeFaceID_i1 TINYINT NOT NULL
		COMMENT "Neighbor on grid when i=2N-1; this is where axis j1 is equal to +1",
	CubeFaceID_j0 TINYINT NOT NULL
		COMMENT "Neighbor on grid when j=0; this is where axis j2 is equal to -1",
	CubeFaceID_j1 TINYINT NOT NULL
		COMMENT "Neighbor on grid when j=2N-1; this is where axis j3 is equal to +1"
)
ENGINE='Aria' TRANSACTIONAL=1
COMMENT = "Relate a cube face to its four neighbors based on the direction of approach.";

INSERT INTO KS.CubeFaceNeighbor
(CubeFaceID, CubeFaceID_i0, CubeFaceID_i1, CubeFaceID_j0, CubeFaceID_j1)
SELECT
	f.CubeFaceID,
	-- f.CubeFaceCD,
	f_j1_lo.CubeFaceID AS CubeFaceID_i0,
	f_j1_hi.CubeFaceID AS CubeFaceID_i1,
	f_j2_lo.CubeFaceID AS CubeFaceID_j0,
	f_j2_hi.CubeFaceID AS CubeFaceID_j1
FROM
	KS.CubeFace AS f
	INNER JOIN KS.CubeFace AS f_j1_lo ON f_j1_lo.i=f.j1 AND f_j1_lo.ci=-1
	INNER JOIN KS.CubeFace AS f_j1_hi ON f_j1_hi.i=f.j1 AND f_j1_hi.ci=1
	INNER JOIN KS.CubeFace AS f_j2_lo ON f_j2_lo.i=f.j2 AND f_j2_lo.ci=-1
	INNER JOIN KS.CubeFace AS f_j2_hi ON f_j2_hi.i=f.j2 AND f_j2_hi.ci=1;
