-- Drop tables in reverse dependency order
DROP TABLE IF EXISTS KS.BodyCollectionEntry;
DROP TABLE IF EXISTS KS.BodyCollection;

-- BodyCollection
CREATE OR REPLACE TABLE KS.BodyCollection(
	BodyCollectionID SMALLINT NOT NULL PRIMARY KEY,
	BodyCollectionCD VARCHAR(4) NOT NULL,
	BodyCollectionName VARCHAR(32) NOT NULL,
	BodyCount SMALLINT NOT NULL DEFAULT 0,
	TotalMass DOUBLE NOT NULL DEFAULT 0.0,
	UNIQUE KEY UNQ_BodyCollectionCD (BodyCollectionCD),
	UNIQUE KEY UNQ_BodyCollectionName (BodyCollectionName)
)
COMMENT "Collections of bodies used in solar system integrations.";

-- Populate BodyCollection
INSERT INTO KS.BodyCollection
(BodyCollectionID, BodyCollectionCD, BodyCollectionName)
VALUES
(1, 'P', 'Planets'),
(2, 'D', 'DE435'),
(3, 'D4', 'DE435-Top-16'),
(4, 'D5', 'DE435-Top-32'),
(5, 'D6', 'DE435-Top-64'),
(6, 'D7', 'DE435-Top-128'),
(7, 'D8', 'DE435-Top-256');

-- BodyCollectionEntry
CREATE OR REPLACE TABLE KS.BodyCollectionEntry(
	BodyCollectionID SMALLINT NOT NULL,
	BodyID INT NOT NULL,
	BodyNumber SMALLINT NOT NULL,
	PRIMARY KEY (BodyCollectionID, BodyID),
	UNIQUE KEY UNQ_BodyCollectionEntry_BodyCollectionID_BodyNumber (BodyCollectionID, BodyNumber),
	CONSTRAINT FK_BodyCollectionEntry_BodyCollectionID
		FOREIGN KEY (BodyCollectionID) REFERENCES KS.BodyCollection(BodyCollectionID),
	CONSTRAINT FK_BodyCollectionEntry_BodyID
		FOREIGN KEY (BodyID) REFERENCES KS.Body(BodyID)
)
COMMENT "Members of body collections.";

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

-- Populate KS.BodyCollection for the DE-435-Top-32
INSERT INTO KS.BodyCollectionEntry
(BodyCollectionID, BodyID, BodyNumber)
SELECT
	bc.BodyCollectionID,
	mr.BodyID,
	row_number() OVER (ORDER BY b.SortOrder) AS BodyNumber
FROM 
	JPL.MassRank AS mr
	INNER JOIN KS.Body AS b ON b.BodyID = mr.BodyID
	INNER JOIN KS.BodyCollection AS bc ON bc.BodyCollectionName = 'DE435-Top-32'
WHERE mr.MassRank <= 32;

-- Populate KS.BodyCollection for the DE-435-Top-64
INSERT INTO KS.BodyCollectionEntry
(BodyCollectionID, BodyID, BodyNumber)
SELECT
	bc.BodyCollectionID,
	mr.BodyID,
	row_number() OVER (ORDER BY b.SortOrder) AS BodyNumber
FROM 
	JPL.MassRank AS mr
	INNER JOIN KS.Body AS b ON b.BodyID = mr.BodyID
	INNER JOIN KS.BodyCollection AS bc ON bc.BodyCollectionName = 'DE435-Top-64'
WHERE mr.MassRank <= 64;

-- Populate KS.BodyCollection for the DE-435-Top-128
INSERT INTO KS.BodyCollectionEntry
(BodyCollectionID, BodyID, BodyNumber)
SELECT
	bc.BodyCollectionID,
	mr.BodyID,
	row_number() OVER (ORDER BY b.SortOrder) AS BodyNumber
FROM 
	JPL.MassRank AS mr
	INNER JOIN KS.Body AS b ON b.BodyID = mr.BodyID
	INNER JOIN KS.BodyCollection AS bc ON bc.BodyCollectionName = 'DE435-Top-128'
WHERE mr.MassRank <= 128;

-- Populate KS.BodyCollection for the DE-435-Top-256
INSERT INTO KS.BodyCollectionEntry
(BodyCollectionID, BodyID, BodyNumber)
SELECT
	bc.BodyCollectionID,
	mr.BodyID,
	row_number() OVER (ORDER BY b.SortOrder) AS BodyNumber
FROM 
	JPL.MassRank AS mr
	INNER JOIN KS.Body AS b ON b.BodyID = mr.BodyID
	INNER JOIN KS.BodyCollection AS bc ON bc.BodyCollectionName = 'DE435-Top-256'
WHERE mr.MassRank <= 256;

-- Populate total mass of collections
CREATE OR REPLACE TEMPORARY TABLE KS.BodyCollectionTotalMass
SELECT
	bce.BodyCollectionID,
	bce.BodyID,
	count(bce.BodyID) AS BodyCount,
	sum(mb.M) AS TotalMass
FROM
	KS.BodyCollectionEntry AS bce
	INNER JOIN KS.Body AS b ON b.BodyID = bce.BodyID
	INNER JOIN JPL.HorizonsBody AS hb ON hb.BodyID = bce.BodyID
	INNER JOIN JPL.MassiveBody AS mb ON mb.HorizonsBodyID = hb.HorizonsBodyID
GROUP BY bce.BodyCollectionID;

-- Apply the total mass to the BodyCollection table
UPDATE 
	KS.BodyCollection AS bc
	INNER JOIN KS.BodyCollectionTotalMass AS bctm ON bctm.BodyCollectionID = bc.BodyCollectionID
SET
	bc.BodyCount = bctm.BodyCount,
	bc.TotalMass = bctm.TotalMass;

-- Clean up
DROP TEMPORARY TABLE JPL.MassRank;
DROP TEMPORARY TABLE KS.BodyCollectionTotalMass;