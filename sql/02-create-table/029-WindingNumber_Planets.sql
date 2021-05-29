-- Times when orbital elements of planets wind around
CREATE OR REPLACE TABLE KS.WindingNumber_Planets(
	BodyID INT NOT NULL
		COMMENT "The Body whose orbital elements are described; FK to KS.Body",
	TimeID INT NOT NULL
		COMMENT "Integer ID for the timestamp of these orbital elements; FK to KS.IntegrationTime",
    WindingNumber INT NOT NULL,
	PRIMARY KEY (BodyID, TimeID),
	UNIQUE KEY UNQ_BodyID_WindingNumber(BodyID, WindingNumber),
	CONSTRAINT FK_WindingTime_Planets_TimeID
		FOREIGN KEY (TimeID) REFERENCES KS.IntegrationTime(TimeID),
	CONSTRAINT FK_WindingTime_Planets_BodyID
		FOREIGN KEY (BodyID) REFERENCES KS.Body(BodyID)
)
ENGINE='Aria' TRANSACTIONAL=0
COMMENT "Times where the mean anomaly M of planetary orbital elements windw from positive to negative values.";
