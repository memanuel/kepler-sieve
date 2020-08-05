-- Drop tables in reverse dependency order
DROP TABLE IF EXISTS KS.BodyCollectionEntry;
DROP TABLE IF EXISTS KS.BodyCollection;

-- BodyCollection
CREATE OR REPLACE TABLE KS.BodyCollection(
	BodyCollectionID SMALLINT NOT NULL PRIMARY KEY,
	BodyCollectionName VARCHAR(32) NOT NULL,
	UNIQUE KEY UNQ_BodyCollectionName (BodyCollectionName)
)
COMMENT "Collections of bodies used in solar system integrations.";

-- Populate BodyCollection
INSERT INTO KS.BodyCollection
(BodyCollectionID, BodyCollectionName)
VALUES
(1, 'Planets'),
(2, 'DE-435');

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

CREATE OR REPLACE TEMPORARY TABLE KS.BodyCollectionEntryInput(
	BodyCollectionName VARCHAR(32) NOT NULL,
	BodyName VARCHAR(32) NOT NULL,
	BodyNumber SMALLINT NOT NULL,
	PRIMARY KEY (BodyCollectionName, BodyName)
);

# Populate BodyCollectionEntryInput
-- INSERT INTO KS.BodyCollectionEntryInput
-- (BodyCollectionName, BodyName, BodyNumber)
-- VALUES
-- ('Planets', 'LB.Sun', 1),
-- ('Planets', 'LB.Mercury Barycenter', 2),
-- ('Planets', 'LB.Venus Barycenter', 3),
-- ('Planets', 'LB.Earth', 4),
-- ('Planets', 'LB.Moon', 5),
-- ('Planets', 'LB.Mars Barycenter', 6),
-- ('Planets', 'LB.Jupiter Barycenter', 7),
-- ('Planets', 'LB.Saturn Barycenter', 8),
-- ('Planets', 'LB.Uranus Barycenter', 9),
-- ('Planets', 'LB.Neptune Barycenter', 10),
-- ('Planets', 'LB.Pluto Barycenter', 11);

# Populate KS.BodyCollectionEntry from the input table with string IDs
-- INSERT INTO KS.BodyCollectionEntry
-- (BodyCollectionID, BodyID, BodyNumber)
-- SELECT
-- 	bc.BodyCollectionID,
-- 	b.BodyID,
-- 	bcei.BodyNumber
-- FROM 
-- 	KS.BodyCollectionEntryInput AS bcei
-- 	LEFT JOIN KS.BodyCollection AS bc ON bc.BodyCollectionName = bcei.BodyCollectionName
-- 	LEFT JOIN KS.Body AS b ON b.BodyName = bcei.BodyName
-- ORDER BY bcei.BodyCollectionName, BodyNumber;

# Populate KS.BodyCollection for Planets collection
INSERT INTO KS.BodyCollectionEntry
(BodyCollectionID, BodyID, BodyNumber)
SELECT
	bc.BodyCollectionID,
	b.BodyID,
	row_number() OVER (ORDER BY b.BodyID) AS BodyNumber
FROM
	KS.Body AS b
	INNER JOIN KS.BodyType AS bt ON bt.BodyTypeID = b.BodyTypeID 
	INNER JOIN KS.BodyCollection AS bc ON bc.BodyCollectionName = 'Planets'
WHERE 
(
bt.BodyTypeCD IN ('S', 'PS') OR
(b.BodyName IN ('LB.Earth', 'LB.Moon'))
);

# Populate KS.BodyCollection for the DE-435 
INSERT INTO KS.BodyCollectionEntry
(BodyCollectionID, BodyID, BodyNumber)
SELECT
	bc.BodyCollectionID,
	b.BodyID,
	row_number() OVER (ORDER BY b.BodyID) AS BodyNumber
FROM
	JPL.MassiveBody AS mb
	INNER JOIN JPL.HorizonsBody AS hb ON hb.HorizonsBodyID = mb.HorizonsBodyID
	INNER JOIN KS.Body AS b ON b.BodyID = hb.BodyID
	INNER JOIN KS.BodyCollection AS bc ON bc.BodyCollectionName = 'DE-435'