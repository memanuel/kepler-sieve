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
