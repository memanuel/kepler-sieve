CREATE OR REPLACE TABLE KS.IntegrationDiff(
	BodyCollectionID SMALLINT NOT NULL COMMENT "Collection of bodies that were integrated in Rebound; FK to KS.BodyCollection",
	TimeID INT NOT NULL COMMENT "Time when state vectors are compared between Rebound and JPL; FK to KS.IntegrationTime",
	BodyID INT NOT NULL COMMENT "Body whose state vectors are compared between Rebound and JPL; FK to KS.Body",
	MJD DOUBLE NOT NULL COMMENT "The Modified Julian Date as of which the state vectors are compared",
	dq DOUBLE NOT NULL COMMENT "Difference in position between Rebound and JPL, in AU",
	dv DOUBLE NOT NULL COMMENT "Difference in velocity between Rebound and JPL, in AU/day",
	PRIMARY KEY (BodyCollectionID, TimeID, BodyID) 
		COMMENT "One integration difference is defined by the collection that was integrated, the time stamp, and the body whose difference is measured.",
	CONSTRAINT FK_IntegrationDiff_BodyCollectionID FOREIGN KEY (BodyCollectionID) REFERENCES KS.BodyCollection(BodyCollectionID),
	CONSTRAINT FK_IntegrationDiff_TimeID FOREIGN KEY (TimeID) REFERENCES KS.IntegrationTime(TimeID),
	CONSTRAINT FK_IntegrationDiff_BodyID FOREIGN KEY (BodyID) REFERENCES KS.Body(BodyID)	
)
COMMENT "Compare the state vectors (position and velocity) for Solar System bodies computed in Rebound using various collections of bodies to JPL results.";
