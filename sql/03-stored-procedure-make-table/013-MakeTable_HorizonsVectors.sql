DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE JPL.MakeTable_HorizonsVectors()
COMMENT "Populate the HorizonsVectors table from HorizonsImport"
BEGIN 

INSERT IGNORE INTO JPL.HorizonsVectors 
(TimeID, HorizonsBodyID, mjd, qx, qy, qz, vx, vy, vz)
SELECT 
	ht.TimeID,
	hb.HorizonsBodyID,
	ht.mjd,
	hi.qx,
	hi.qy,
	hi.qz,
	hi.vx,
	hi.vy,
	hi.vz
FROM
	JPL.HorizonsImport AS hi
	INNER JOIN KS.BodyType AS bt ON
		bt.BodyTypeCD = hi.BodyTypeCD
	INNER JOIN JPL.HorizonsBody AS hb ON
		hb.BodyTypeID = bt.BodyTypeID AND
		hb.HorizonsBodyNumber = hi.BodyNumber
	INNER JOIN JPL.HorizonsTime AS ht ON
		ht.JD = hi.JD;

END
$$

DELIMITER ;
