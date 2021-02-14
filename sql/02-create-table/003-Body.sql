-- Drop tables to clear dependencies
DROP TABLE IF EXISTS JPL.HorizonsBody;
DROP TABLE IF EXISTS JPL.LargeBody;
DROP TABLE IF EXISTS JPL.SmallBody;

-- Body
CREATE OR REPLACE TABLE KS.Body(
	BodyID INT NOT NULL PRIMARY KEY,
	BodyName VARCHAR(32) NOT NULL,
	BodyTypeID TINYINT NOT NULL,
    SortOrder INT NOT NULL DEFAULT 0,
	CONSTRAINT FK_Body_BodyTypeID FOREIGN KEY (BodyTypeID) REFERENCES BodyType(BodyTypeID),
	UNIQUE KEY UNQ_BodyTypeID_BodyName (BodyTypeID, BodyName)
)
ENGINE='Aria' TRANSACTIONAL=1
COMMENT "Solar System bodies used in the Kepler Sieve application.";

-- Large Body
CREATE OR REPLACE TABLE JPL.LargeBody(
	LargeBodyID SMALLINT NOT NULL PRIMARY KEY,
	LargeBodyName VARCHAR(32) NOT NULL UNIQUE,
	BodyTypeID TINYINT NOT NULL,
	BodyID INT NOT NULL UNIQUE,
	CONSTRAINT FK_LargeBody_BodyTypeID FOREIGN KEY (BodyTypeID) REFERENCES KS.BodyType(BodyTypeID),
	CONSTRAINT FK_LargeBody_BodyID FOREIGN KEY (BodyID) REFERENCES KS.Body(BodyID)
)
ENGINE='Aria' TRANSACTIONAL=1
COMMENT "Large Body as defined by JPL / Horizons system.  Includes stars, planets, moons.";

-- Small Body
CREATE OR REPLACE TABLE JPL.SmallBody(
	SmallBodyID INT NOT NULL PRIMARY KEY,
	SmallBodyName VARCHAR(32) NOT NULL UNIQUE,
	BodyTypeID TINYINT NOT NULL,
	BodyID INT NOT NULL UNIQUE,
	CONSTRAINT FK_SmallBody_BodyTypeID FOREIGN KEY (BodyTypeID) REFERENCES KS.BodyType(BodyTypeID),
	CONSTRAINT FK_SmallBody_BodyID FOREIGN KEY (BodyID) REFERENCES KS.Body(BodyID)
)
ENGINE='Aria' TRANSACTIONAL=1
COMMENT "Small Body as defined by JPL / Horizons system.  Includes asteroids.";

-- HorizonsBody
CREATE OR REPLACE TABLE JPL.HorizonsBody(
	HorizonsBodyID INT NOT NULL PRIMARY KEY,
	HorizonsBodyNumber INT NOT NULL,
	HorizonsBodyName VARCHAR(32) NOT NULL UNIQUE,
	BodyTypeID TINYINT NOT NULL,
	BodyID INT NOT NULL UNIQUE,
	LargeBodyID SMALLINT NULL,
	SmallBodyID INT NULL,
	UNIQUE KEY UNQ_BodyTypeID_BodyNumber (BodyTypeID, HorizonsBodyNumber),
	UNIQUE KEY UNQ_BodyTypeID_BodyName (BodyTypeID, HorizonsBodyName),
	CONSTRAINT FK_HorizonsBody_BodyTypeID FOREIGN KEY (BodyTypeID) REFERENCES KS.BodyType(BodyTypeID),
	CONSTRAINT FK_HorizonsBody_BodyID FOREIGN KEY (BodyID) REFERENCES KS.Body(BodyID),
	CONSTRAINT FK_HorizonsBody_LargeBodyID FOREIGN KEY (LargeBodyID) REFERENCES JPL.LargeBody(LargeBodyID),
	CONSTRAINT FK_HorizonsBody_SmallBodyID FOREIGN KEY (SmallBodyID) REFERENCES JPL.SmallBody(SmallBodyID)
)
ENGINE='Aria' TRANSACTIONAL=1
COMMENT "Solar System bodies with data downloaded from Horizons.";
