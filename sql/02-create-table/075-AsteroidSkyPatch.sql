CREATE OR REPLACE TABLE KS.AsteroidSkyPatch(
	-- The time when the direction is calculated
	IntegrationTimeID INT NOT NULL
		COMMENT "The TimeID when the direction is calculated; sampled from TimeSliceID; foreign key to IntegrationTime table",
	-- The asteroid
	AsteroidID INT NOT NULL
		COMMENT "The asteroid whose direction SkyPatch is being saved; foreign key to Asteroid table",
	-- The SkyPatch of the direction
	SkyPatchID INT NOT NULL
		COMMENT "The SkyPatch where this asteroid could be found at the given IntegrationTime; foreign key to SkyPatch table",
	-- Primary key
	PRIMARY KEY (IntegrationTimeID, AsteroidID),
	-- Also key the table by time and direction for searching later (?)
	-- Foreign keys
	CONSTRAINT FK_AsteroidSkyPatch_IntegrationTimeID FOREIGN KEY (IntegrationTimeID) REFERENCES KS.IntegrationTime(IntegrationTimeID),
	CONSTRAINT FK_AsteroidSkyPatch_AsteroidID FOREIGN KEY (AsteroidID) REFERENCES KS.Asteroid(AsteroidID),
	CONSTRAINT FK_AsteroidSkyPatch_SkyPatchID FOREIGN KEY (SkyPAtchID) REFERENCES KS.SkyPatch(SkyPatchID)
)
ENGINE='Aria' TRANSACTIONAL=0
COMMENT = "For all known asteroids, at every IntegrationTimeSlice, compute the direction of the asteroid and its SkyPatch.";