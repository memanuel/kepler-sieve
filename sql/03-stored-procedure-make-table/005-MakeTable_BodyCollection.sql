DELIMITER $$

CREATE OR REPLACE 
DEFINER = kepler
PROCEDURE JPL.MakeTable_MassiveBody()
COMMENT "Populate the MassiveBody table from MassiveBodyImport"
BEGIN 

-- Empty BodyCollectionEntry and BodyCollection
TRUNCATE TABLE BodyCollectionEntry;
TRUNCATE TABLE BodyCollection;
	
-- Populate BodyCollection
INSERT INTO KS.BodyCollection
(BodyCollectionID, BodyCollectionCD, BodyCollectionName)
VALUES
(1, 'P', 'Planets'),
(2, 'D', 'DE435'),
(3, 'D4', 'DE435-Top-16'),
(101, 'PS1', 'Mercury System'),
(102, 'PS2', 'Venus System'),
(103, 'PS3', 'Earth System'),
(104, 'PS4', 'Mars System'),
(105, 'PS5', 'Jupiter System'),
(106, 'PS6', 'Saturn System'),
(107, 'PS7', 'Uranus System'),
(108, 'PS8', 'Neptune System'),
(109, 'PS9', 'Pluto System');

-- Temporary table with ranks of massive bodies sorted by mass
CREATE OR REPLACE TEMPORARY TABLE JPL.MassRank
SELECT
	b.BodyName,
	b.BodyID,
	(mb.M / @SolarSystemMass) AS M_ss,
	(mb.M / @AsteroidMass) AS M_ast,
	row_number() OVER (ORDER BY mb.M DESC) AS MassRank
FROM
	JPL.MassiveBody AS mb
	INNER JOIN JPL.HorizonsBody AS hb ON hb.HorizonsBodyID = mb.HorizonsBodyID
	INNER JOIN KS.Body AS b ON b.BodyID = hb.BodyID
-- Avoid the EMB because we want to treat the Earth and Moon separately	
WHERE b.BodyName != 'Earth-Moon Barycenter';

-- Populate KS.BodyCollection for the planets
INSERT INTO KS.BodyCollectionEntry
(BodyCollectionID, BodyID, BodyNumber)
SELECT
	bc.BodyCollectionID,
	mr.BodyID,
	row_number() OVER (ORDER BY b.SortOrder) AS BodyNumber
FROM 
	JPL.MassRank AS mr
	INNER JOIN KS.Body AS b ON b.BodyID = mr.BodyID
	INNER JOIN KS.BodyType AS bt ON bt.BodyTypeID = b.BodyTypeID
	INNER JOIN KS.BodyCollection AS bc ON bc.BodyCollectionName = 'Planets'
WHERE bt.IsLargeBody_JPL;

-- Populate KS.BodyCollection for the DE435 
INSERT INTO KS.BodyCollectionEntry
(BodyCollectionID, BodyID, BodyNumber)
SELECT
	bc.BodyCollectionID,
	mr.BodyID,
	row_number() OVER (ORDER BY b.SortOrder) AS BodyNumber
FROM 
	JPL.MassRank AS mr
	INNER JOIN KS.Body AS b ON b.BodyID = mr.BodyID
	INNER JOIN KS.BodyCollection AS bc ON bc.BodyCollectionName = 'DE435'
;

-- Populate KS.BodyCollection for the DE435-Top-16
INSERT INTO KS.BodyCollectionEntry
(BodyCollectionID, BodyID, BodyNumber)
SELECT
	bc.BodyCollectionID,
	mr.BodyID,
	row_number() OVER (ORDER BY b.SortOrder) AS BodyNumber
FROM 
	JPL.MassRank AS mr
	INNER JOIN KS.Body AS b ON b.BodyID = mr.BodyID
	INNER JOIN KS.BodyCollection AS bc ON bc.BodyCollectionName = 'DE435-Top-16'
WHERE mr.MassRank <= 16;

-- Populate planetary system collections
INSERT INTO KS.BodyCollectionEntry
(BodyCollectionID, BodyID, BodyNumber)
SELECT
	bc.BodyCollectionID,
	b.BodyID,
	-- Use mod 99 trick so the planet with k=99 appears with body number 1, then the moons in order
	(k._ % 99) + 1 AS BodyNumber
FROM
	-- Counter from 1 to 9
	KS.Counter AS n
	-- The body number in the system according to horizons; moons start at 1, planet is 99
	INNER JOIN KS.Counter AS k ON k._ <= 99
	-- Bodies in this planetary system
	INNER JOIN KS.Body AS b ON b.BodyID = n._*100 + k._
	-- Planetary system containing this planet and its satellites
	INNER JOIN KS.BodyCollection AS bc ON bc.BodyCollectionCD = CONCAT('PS', n._)
WHERE n._ <= 9
ORDER BY BodyCollectionID, BodyNumber;

-- Apply the total mass to the BodyCollection table
-- First handle sum of bodies whose mass is listed on MassiveBody
CREATE TEMPORARY TABLE KS.bctm AS
SELECT
	bc.BodyCollectionID,
	count(bce.BodyID) AS BodyCount,
	sum(mb.M) AS TotalMass,
	sum(mb.GM) AS TotalGM
FROM
	KS.BodyCollection AS bc
	INNER JOIN KS.BodyCollectionEntry AS bce ON bce.BodyCollectionID = bc.BodyCollectionID
	INNER JOIN KS.Body AS b ON b.BodyID = bce.BodyID
	INNER JOIN KS.MassiveBody AS mb ON mb.BodyID = b.BodyID
GROUP BY bc.BodyCollectionID;

UPDATE 
	KS.BodyCollection AS bc
	INNER JOIN KS.bctm ON bctm.BodyCollectionID = bc.BodyCollectionID
SET
	bc.BodyCount = bctm.BodyCount,
	bc.TotalMass = bctm.TotalMass;

-- Next handle the planetary systems other than Earth
UPDATE
	Counter AS n
	INNER JOIN KS.BodyCollection AS bc ON bc.BodyCollectionCD = CONCAT('PS', n._)
	INNER JOIN KS.MassiveBody AS mb ON mb.BodyID = n._
SET	
	bc.TotalMass = mb.M	
WHERE n._ BETWEEN 1 AND 9 AND n._ <> 3;	

-- Populate BarycenterWeight table
TRUNCATE TABLE KS.BarycenterWeight;

INSERT INTO KS.BarycenterWeight
(BodyCollectionID, BodyID, M, Weight)
SELECT
	bc.BodyCollectionID,
	b.BodyID,
	COALESCE(mb.M, 0.0) AS M,
	COALESCE(mb.M, 0.0) / bc.TotalMass AS Weight
FROM
	KS.BodyCollection AS bc
	INNER JOIN KS.BodyCollectionEntry AS bce ON	bce.BodyCollectionID = bc.BodyCollectionID
	INNER JOIN KS.Body AS b ON b.BodyID= bce.BodyID
	LEFT JOIN KS.MassiveBody AS mb ON mb.BodyID = b.BodyID
;

-- Clean up
DROP TEMPORARY TABLE JPL.MassRank;
DROP TEMPORARY TABLE KS.bctm;
	
END
$$

DELIMITER ;
