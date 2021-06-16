-- ****************************************************************************************************
CREATE OR REPLACE TABLE KS.AsteroidSkyPatch(
	-- The asteroid
	AsteroidID INT NOT NULL
		COMMENT "The asteroid whose direction SkyPatch is being saved; foreign key to Asteroid table",
	Segment INT NOT NULL
		COMMENT "The segment number; table describes segments where asteroid is one one patch of the sky",
	-- The SkyPatch of the direction
	SkyPatchID INT NOT NULL
		COMMENT "The SkyPatch where this asteroid could be found during this segment; foreign key to SkyPatch table",
	-- The time when the direction is calculated
	TimeID_0 INT NOT NULL
		COMMENT "The TimeID when this segment starts",
	TimeID_1 INT NOT NULL
		COMMENT "The TimeID when this segment ends",
	-- Primary key
	PRIMARY KEY (AsteroidID, Segment),
	-- Also key the table by SkyPatchID and TimeID to facilitate searching for a SkyPatchID and TimeID to match a detection
	INDEX IDX_SkyPatchID (SkyPatchID),
	-- Foreign keys
	CONSTRAINT FK_AsteroidSkyPatch_AsteroidID FOREIGN KEY (AsteroidID) REFERENCES KS.Asteroid(AsteroidID),
	CONSTRAINT FK_AsteroidSkyPatch_SkyPatchID FOREIGN KEY (SkyPatchID) REFERENCES KS.SkyPatch(SkyPatchID)
    -- Don't build these FKs because don't want to bloat table with superfluous indices
	-- CONSTRAINT FK_AsteroidSkyPatch_TimeID_0 FOREIGN KEY (TimeID_0) REFERENCES KS.HiResTime(TimeID),
	-- CONSTRAINT FK_AsteroidSkyPatch_TimeID_1 FOREIGN KEY (TimeID_1) REFERENCES KS.HiResTime(TimeID)	
)
ENGINE='Aria' TRANSACTIONAL=0
COMMENT = "For all known asteroids, at every IntegrationTimeSlice, compute the direction of the asteroid and its SkyPatch.";

-- ****************************************************************************************************
CREATE OR REPLACE TABLE KS.AsteroidSkyPatch_Stage(
	-- The asteroid
	AsteroidID INT NOT NULL
		COMMENT "The asteroid whose direction SkyPatch is being saved; foreign key to Asteroid table",
	Segment INT NOT NULL
		COMMENT "The segment number; table describes segments where asteroid is one one patch of the sky",
	-- The SkyPatch of the direction
	SkyPatchID INT NOT NULL
		COMMENT "The SkyPatch where this asteroid could be found during this segment; foreign key to SkyPatch table",
	-- The time when the direction is calculated
	TimeID_0 INT NOT NULL
		COMMENT "The TimeID when this segment starts",
	TimeID_1 INT NOT NULL
		COMMENT "The TimeID when this segment ends",
	-- Primary key
	PRIMARY KEY (AsteroidID, Segment)
)
ENGINE='Aria' TRANSACTIONAL=0
COMMENT = "Staging table for AsteroidSkypatch";
