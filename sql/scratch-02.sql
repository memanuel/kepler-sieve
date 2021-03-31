SELECT * FROM KS.Tracklet;
SELECT * FROM KS.TrackletBatch;

CREATE OR REPLACE TEMPORARY TABLE KS.TrackletBatch LIKE KS.Tracklet;

INSERT INTO KS.TrackletBatch
(DetectionID_1, DetectionID_2, SkyPatchID, mjd, ux, uy, uz, ux_dot, uy_dot, uz_dot)

SELECT
	-- The two DetectionIDs	
	d1.DetectionID AS DetectionID_1,
	d2.DetectionID AS DetectionID_2,
	-- The SkyPatchID of the mean
	sp.SkyPatchID,
	-- The mean date
	dtp.mjd,
	-- Average direction
	(d1.ux + d2.ux)/2 AS ux,
	(d1.uy + d2.uy)/2 AS uy,
	(d1.uz + d2.uz)/2 AS uz,
	-- Derivative of direction
	(d2.ux - d2.ux)/dtp.dt AS ux_dot,
	(d2.uy - d2.uy)/dtp.dt AS uy_dot,
	(d2.uz - d2.uz)/dtp.dt AS uz_dot
FROM
	KS.Detection AS d1
	INNER JOIN KS.DetectionTimePair AS dtp ON dtp.DetectionTimeID_1 = d1.DetectionTimeID
	INNER JOIN KS.Detection AS d2 ON 
		d2.DetectionTimeID = dtp.DetectionTimeID_2 AND 
		d2.SkyPatchID = d1.SkyPatchID
	-- Possible values of the SkyPatch
	INNER JOIN KS.SkyPatch AS sp ON sp.SkyPatchID IN (d1.SkyPatchID, d2.SkyPatchID)
WHERE
	d1.SkyPatchID=d2.SkyPatchID
LIMIT 100;
