SET @N = 1024;
SET @M = 2*@N;
SET @dr_max = 0.25/@N;

WITH t1 AS (
SELECT
	db.DetectionBatchID,
	db.DetectionTimeID,
	db.DetectionID,
	cf.CubeFaceID,
	db.ux AS x,
	db.uy AS y,
	db.uz AS z,
	-- (CASE cf.i WHEN 1 THEN db.ux WHEN 2 THEN db.uy WHEN 3 THEN db.uz END)*cf.ci AS ww,
	-- w is the absolute value of the largest component
	greatest(db.ux, db.uy, db.uz, -db.ux, -db.uy, -db.uz) AS w,
	1.0 / greatest(db.ux, db.uy, db.uz, -db.ux, -db.uy, -db.uz) AS t,
	-- calculate u and v
	CASE cf.alpha WHEN 'X' THEN db.ux WHEN 'Y' THEN db.uy WHEN 'Z' THEN db.uz END AS u,
	CASE cf.beta  WHEN 'X' THEN db.ux WHEN 'Y' THEN db.uy WHEN 'Z' THEN db.uz END AS v,
	db.mag
FROM
	KS.DetectionBatch AS db
	INNER JOIN KS.CubeFace AS cf ON
		(CASE cf.i WHEN 1 THEN db.ux WHEN 2 THEN db.uy WHEN 3 THEN db.uz END)*cf.ci = 
		greatest(db.ux, db.uy, db.uz, -db.ux, -db.uy, -db.uz)
-- GROUP BY db.DetectionBatchID		
)
SELECT
	t1.DetectionTimeID,
	t1.DetectionID,
	t1.CubeFaceID,
	t1.x,
	t1.y,
	t1.z,
	t1.mag,
	t1.t,
	t1.u,
	t1.v,
	t1.w,
	t1.t * t1.u AS a,
	t1.t * t1.v AS b,
	t1.t * t1.w AS c	
FROM
	t1;
