-- ************************************************************************************************
-- The eight vertices (corners) of a cube spanning [-1,+1]^3
CREATE OR REPLACE TABLE KS.CubeVertex(
    CubeVertexID TINYINT NOT NULL PRIMARY KEY,
    x DOUBLE NOT NULL,
    y DOUBLE NOT NULL,
    z DOUBLE NOT NULL,
    UNIQUE KEY UNQ_CubeVertex_x_y_z (x, y, z)
)
ENGINE='Aria' TRANSACTIONAL=1
COMMENT = "The eight vertices of a cube where each coordinate is at +/- 1.";

INSERT INTO KS.CubeVertex
(CubeVertexID, x, y, z)
WITH s AS (
SELECT
	CAST(2*x._-1.0 AS DOUBLE) AS t
FROM
	KS.Counter AS x
WHERE x._ < 2
)
SELECT
	row_number() OVER (ORDER BY s1.t, s2.t, s3.t) AS CubeVertexID,
	s1.t AS x,
	s2.t AS y,
	s3.t AS z
FROM 
	s AS s1,
	s AS s2,
	s AS s3
ORDER BY s1.t, s2.t, s3.t;

-- ************************************************************************************************
CREATE OR REPLACE TABLE KS.CubeEdge(
	CubeEdgeID TINYINT NOT NULL PRIMARY KEY,
	CubeVertexID_1 TINYINT NOT NULL
		COMMENT "Vertex of this edge with smaller VertexID",
	CubeVertexID_2 TINYINT NOT NULL
		COMMENT "Vertex of this edge with larger VertexID",
    i TINYINT NOT NULL
    	COMMENT "Index of the axis that varying in this edge, e.g. i=1 means that x covers [-1,+1]",
    j1 TINYINT NOT NULL
    	COMMENT "Index of first axis that is constant on this edge; j1=i+1 mod 3.",
    j2 TINYINT NOT NULL
    	COMMENT "Index of second axis that is constant on this edge; j2=i+2 mod 3.",
    cj1 DOUBLE NOT NULL
    	COMMENT "Constant value of axis j1 on this edge",
    cj2 DOUBLE NOT NULL
    	COMMENT "Constant value of axis j2 on this edge",
	UNIQUE KEY UNQ_CubeEdge_CubeVertexID_Pair (CubeVertexID_1, CubeVertexID_2)
)
ENGINE='Aria' TRANSACTIONAL=1
COMMENT = "The 12 edges of a cube characterized by the pair of vertices sorted in ascending order.";

INSERT INTO KS.CubeEdge
(CubeEdgeID, CubeVertexID_1, CubeVertexID_2, i, j1, j2, cj1, cj2)
WITH t1 AS(
SELECT
	(i._+1) AS i,
	((i._+1) MOD 3) + 1 AS j1,
	((i._+2) MOD 3) + 1 AS j2,
	2*s1._-1 AS c1,
	2*s2._-1 AS c2
FROM
	KS.Counter AS i,
	KS.Counter AS s1,
	KS.Counter AS s2
WHERE
	i._ < 3 AND s1._<2 AND s2._<2
), t2 AS(
SELECT
	t1.i,
	t1.j1,
	t1.j2,
	t1.c1,
	t1.c2,
	-- Coordinates of the first vertex (x1, y1, z1)
	IF(t1.i=1,-1,0) + IF(t1.j1=1,t1.c1,0) + IF(t1.j2=1,t1.c2,0) AS x1,
	IF(t1.i=2,-1,0) + IF(t1.j1=2,t1.c1,0) + IF(t1.j2=2,t1.c2,0) AS y1,
	IF(t1.i=3,-1,0) + IF(t1.j1=3,t1.c1,0) + IF(t1.j2=3,t1.c2,0) AS z1,
	-- Coordinates of the second vertex (x2, y2, z2)
	IF(t1.i=1,1,0) + IF(t1.j1=1,t1.c1,0) + IF(t1.j2=1,t1.c2,0) AS x2,
	IF(t1.i=2,1,0) + IF(t1.j1=2,t1.c1,0) + IF(t1.j2=2,t1.c2,0) AS y2,
	IF(t1.i=3,1,0) + IF(t1.j1=3,t1.c1,0) + IF(t1.j2=3,t1.c2,0) AS z2
FROM 
	t1
)
SELECT
	row_number() OVER (ORDER BY i, j1, j2, cj1, cj2) AS CubeEdgeID,
	cv1.CubeVertexID AS CubeVertexID_1,
	cv2.CubeVertexID AS CubeVertexID_2,
	t2.i,
	t2.j1,
	t2.j2,
	t2.c1 AS cj1,
	t2.c2 AS cj2
FROM
	t2
	INNER JOIN KS.CubeVertex AS cv1 ON cv1.x = t2.x1 AND cv1.y= t2.y1 AND cv1.z=t2.z1
	INNER JOIN KS.CubeVertex AS cv2 ON cv2.x = t2.x2 AND cv2.y= t2.y2 AND cv2.z=t2.z2
ORDER BY i, j1, j2, cj1, cj2;

-- ************************************************************************************************
-- The six faces of a cube
CREATE OR REPLACE TABLE KS.CubeFace(
    CubeFaceID TINYINT NOT NULL PRIMARY KEY
    	COMMENT "CubeFaceID is chosen to be consistent with a six sided die; orient it with 1 on top, 5 facing you, 3 on right.",
    CubeFaceCD CHAR(2) NOT NULL UNIQUE,
    i TINYINT NOT NULL
    	COMMENT "Index of the axis that is constant on this face, e.g. i=3 is for the face Z+",
    j1 TINYINT NOT NULL
    	COMMENT "Index of the first axis that is varying on this face, e.g. j1=1 for x on Z+",
    j2 TINYINT NOT NULL
    	COMMENT "Index of the second axis that is varying on this face, e.g. j2=2 for y on Z+",
    ci DOUBLE NOT NULL
    	COMMENT "Constant value of axis i on this face; one of -1 and +1; e.g. ci=+1 on Z+",
    alpha CHAR(1) NOT NULL
    	COMMENT "The first axis name that varies on this face; corresponds to j1",
    beta CHAR(1) NOT NULL
    	COMMENT "The second axis name that varies on this face; corresponds to j2",
    gamma CHAR(1) NOT NULL
    	COMMENT "The axis name that is +/- 1 on this face; corresponds to i"    	
)
ENGINE='Aria' TRANSACTIONAL=1
COMMENT = "The six faces of a cube described as the index of the constant axis and its value.";

INSERT INTO KS.CubeFace
(CubeFaceID, CubeFaceCD, i, j1, j2, ci, alpha, beta, gamma)
WITH t1 AS(
SELECT
	(i._+1) AS i,
	((i._+1) MOD 3) + 1 AS j1,
	((i._+2) MOD 3) + 1 AS j2,	
	CAST(2*s._ - 1 AS DOUBLE) AS c
FROM
	KS.Counter AS i
	INNER JOIN KS.Counter AS s ON s._ < 2
WHERE
	i._ < 3
)
SELECT
	-- Integer and string ID
	row_number() OVER (ORDER BY t1.i*t1.c DESC) AS CubeFaceID,
	CONCAT(CHAR(ASCII('X')+t1.i-1), IF(c>0, '+', '-')) AS CubeFaceCD,
	-- Index numbers for the axes and constant value on face
	t1.i,
	t1.j1,
	t1.j2,
	t1.c AS ci,
	-- Names of the three axes
	CHAR(ASCII('W')+t1.j1) AS alpha,
	CHAR(ASCII('W')+t1.j2) AS beta,
	CHAR(ASCII('W')+t1.i) AS gamma	
FROM
	t1
ORDER BY (t1.i*t1.c) DESC;
