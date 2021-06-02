-- Interactions where a detection is near the direction from a known asteroid to earth
CREATE OR REPLACE TABLE KS.DetectionNearAsteroid(
    DetectionID INT NOT NULL
        COMMENT "The asteroid detection that is near a known asteroid",
	AsteroidID INT NOT NULL
		COMMENT "The known asteroid that is close in the sky to this detection",
	s DOUBLE NOT NULL
		COMMENT "The Cartesian distance between the detection and the known asteroid",
    LightTime DOUBLE NOT NULL
        COMMENT "Splined light time from asteroid to Earth; used for high precision direction calculation",
	-- Keys and constraints
	PRIMARY KEY (DetectionID, AsteroidID),
	CONSTRAINT FK_DetectionNearAsteroid_DetectionID FOREIGN KEY (DetectionID) REFERENCES KS.Detection(DetectionID),
	CONSTRAINT FK_AsteroidDirections_AsteroidID	FOREIGN KEY (AsteroidID) REFERENCES KS.Asteroid(AsteroidID)
)
ENGINE='Aria' TRANSACTIONAL=0
COMMENT "Interactions between detections and asteroids where a known asteroid is close in the sky to a ."
