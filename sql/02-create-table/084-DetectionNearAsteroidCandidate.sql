-- ********************************************************************************
-- Interactions where a detection is near the direction from a known asteroid to earth
CREATE OR REPLACE TABLE KS.DetectionNearAsteroidCandidate(
	-- Key fields
    DetectionID INT NOT NULL
        COMMENT "The detection that is near a known asteroid",
	AsteroidID INT NOT NULL
		COMMENT "The known asteroid that is close in the sky to this detection",
	-- Keys and constraints
	PRIMARY KEY (DetectionID, AsteroidID)
		COMMENT "Optimize PK to find the nearest asteroid to a detection"
)
ENGINE='Aria' TRANSACTIONAL=0
COMMENT "Interactions between detections and asteroids where a known asteroid is close in the sky to a detection; not yet enriched with precise distance and light time.";
