-- Drop tables to clear dependencies
DROP TABLE IF EXISTS JPL.HorizonsBody;
DROP TABLE IF EXISTS JPL.LargeBody;
DROP TABLE IF EXISTS JPL.SmallBody;

-- Body
CREATE OR REPLACE TABLE KeplerDB.Body(
	BodyID INT NOT NULL PRIMARY KEY,
	BodyName VARCHAR(32) NOT NULL,
	BodyTypeID TINYINT NOT NULL,
	CONSTRAINT FK_Body_BodyTypeID FOREIGN KEY (BodyTypeID) REFERENCES BodyType(BodyTypeID),
	UNIQUE KEY UNQ_BodyTypeID_BodyName (BodyTypeID, BodyName)
)
	COMMENT "Solar System bodies used in the Kepler Sieve application.";

-- Large Body
CREATE OR REPLACE TABLE JPL.LargeBody(
	LargeBodyID SMALLINT NOT NULL PRIMARY KEY,
	LargeBodyName VARCHAR(32) NOT NULL UNIQUE,
	BodyTypeID TINYINT NOT NULL,
	BodyID INT NOT NULL UNIQUE,
	CONSTRAINT FK_LargeBody_BodyTypeID FOREIGN KEY (BodyTypeID) REFERENCES KeplerDB.BodyType(BodyTypeID),
	CONSTRAINT FK_LargeBody_BodyID FOREIGN KEY (BodyID) REFERENCES KeplerDB.Body(BodyID)
)
	COMMENT "Large Body as defined by JPL / Horizons system.  Includes stars, planets, moons.";

-- Small Body
CREATE OR REPLACE TABLE JPL.SmallBody(
	SmallBodyID INT NOT NULL PRIMARY KEY,
	SmallBodyName VARCHAR(32) NOT NULL UNIQUE,
	BodyTypeID TINYINT NOT NULL,
	BodyID INT NOT NULL UNIQUE,
	CONSTRAINT FK_SmallBody_BodyTypeID FOREIGN KEY (BodyTypeID) REFERENCES KeplerDB.BodyType(BodyTypeID),
	CONSTRAINT FK_SmallBody_BodyID FOREIGN KEY (BodyID) REFERENCES KeplerDB.Body(BodyID)
)
	COMMENT "Small Body as defined by JPL / Horizons system.  Includes asteroids.";

-- HorizonsBody
CREATE OR REPLACE TABLE JPL.HorizonsBody(
	HorizonsBodyID INT NOT NULL PRIMARY KEY,
	HorizonsBodyName VARCHAR(32) NOT NULL UNIQUE,
	BodyTypeID TINYINT NOT NULL,
	BodyID INT NOT NULL UNIQUE,
	LargeBodyID SMALLINT NULL,
	SmallBodyID INT NULL,
	UNIQUE KEY UNQ_BodyTypeID_BodyName (BodyTypeID, HorizonsBodyName),
	CONSTRAINT FK_Body_BodyTypeID FOREIGN KEY (BodyTypeID) REFERENCES KeplerDB.BodyType(BodyTypeID),
	CONSTRAINT FK_HorizonsBody_BodyID FOREIGN KEY (BodyID) REFERENCES KeplerDB.Body(BodyID),
	CONSTRAINT FK_HorizonsBody_LargeBodyID FOREIGN KEY (LargeBodyID) REFERENCES JPL.LargeBody(LargeBodyID),
	CONSTRAINT FK_HorizonsBody_SmallBodyID FOREIGN KEY (SmallBodyID) REFERENCES JPL.SmallBody(SmallBodyID)
)
	COMMENT "Solar System bodies with data downloaded from Horizons.";
