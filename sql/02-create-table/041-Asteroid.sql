CREATE OR REPLACE TABLE KS.Asteroid(
	AsteroidID INT NOT NULL PRIMARY KEY
		COMMENT "Internal integer ID for Asteroids; in practice this is the same as AsteroidNumber field",
	AsteroidNumber INT NOT NULL UNIQUE
		COMMENT "The IAU asteroid number when it exists; otherwise a sequential counter starting at 1,000,000",
	AsteroidName VARCHAR(32) NOT NULL UNIQUE
		COMMENT "The IAU asteroid name where applicable or asteroid designation",
    AsteroidIndex INT NOT NULL AUTO_INCREMENT UNIQUE
        COMMENT "Internal sequence counter of all asteroids with no skips, used for batching jobs",
	BodyID INT NOT NULL UNIQUE
		COMMENT "Foreign key to KS.Body table",
    IsNumberedAsteroid BOOL NOT NULL
        COMMENT "Flag indicating whether or not this is an IAU numbered asteroid",
	H DOUBLE NOT NULL
		COMMENT "The H brightness parameter",
	G DOUBLE NOT NULL
		COMMENT "The G brightness parameter",
	CONSTRAINT FK_Asteroid_BodyID FOREIGN KEY (BodyID) REFERENCES KS.Body(BodyID)
)
ENGINE='Aria' TRANSACTIONAL=1
COMMENT "Census of all known asteroids including reference to KS.Body table.  Orbital Elements stored separately.";
