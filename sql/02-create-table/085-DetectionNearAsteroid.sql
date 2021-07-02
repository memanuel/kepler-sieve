-- ********************************************************************************
-- Interactions where a detection is near the direction from a known asteroid to earth
CREATE OR REPLACE TABLE KS.DetectionNearAsteroid(
	-- Key fields
    DetectionID INT NOT NULL
        COMMENT "The detection that is near a known asteroid",
	AsteroidID INT NOT NULL
		COMMENT "The known asteroid that is close in the sky to this detection",
    -- Data payload
	s DOUBLE NOT NULL
		COMMENT "The Cartesian distance between the detection and the known asteroid",
    LightTime DOUBLE NOT NULL
        COMMENT "Splined light time from asteroid to Earth; used for high precision direction calculation",
	-- Keys and constraints
	PRIMARY KEY (DetectionID, AsteroidID)
		COMMENT "Optimize PK to find the nearest asteroid to a detection",
	UNIQUE KEY (AsteroidID, DetectionID)
		COMMENT "Support querying by asteroid to find all nearby detections"
	-- Don't actually build these FKs, don't want to bloat the table 
	-- CONSTRAINT FK_DetectionNearAsteroid_DetectionID FOREIGN KEY (DetectionID) REFERENCES KS.Detection(DetectionID),
	-- CONSTRAINT FK_AsteroidDirections_AsteroidID	FOREIGN KEY (AsteroidID) REFERENCES KS.Asteroid(AsteroidID)
)
ENGINE='Aria' TRANSACTIONAL=0
COMMENT "Interactions between detections and asteroids where a known asteroid is close in the sky to a detection.";
