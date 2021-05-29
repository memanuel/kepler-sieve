-- Times when orbital elements of planets wind around
CREATE OR REPLACE TABLE KS.AsteroidWindingNumber(
	AsteroidID INT NOT NULL,
	TimeID INT NOT NULL,
    WindingNumber INT NOT NULL,
	PRIMARY KEY (AsteroidID, TimeID),
	UNIQUE KEY UNQ_BodyID_WindingNumber(AsteroidID, WindingNumber)
)
ENGINE='Aria' TRANSACTIONAL=0
COMMENT "Times where the mean anomaly M of asteroid orbital elements windw from positive to negative values.";
