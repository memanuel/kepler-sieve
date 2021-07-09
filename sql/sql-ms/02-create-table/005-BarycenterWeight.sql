DROP TABLE IF EXISTS KS.BarycenterWeight;

-- BarycenterWeight
CREATE OR REPLACE TABLE KS.BarycenterWeight(
	BodyCollectionID SMALLINT NOT NULL,
	BodyID INT NOT NULL,
    M DOUBLE not null,
    Weight DOUBLE not null,
	PRIMARY KEY (BodyCollectionID, BodyID),
	CONSTRAINT FK_BarycenterWeight_BodyCollectionID
		FOREIGN KEY (BodyCollectionID) REFERENCES KS.BodyCollection(BodyCollectionID),
	CONSTRAINT FK_BarycenterWeight_BodyID
		FOREIGN KEY (BodyID) REFERENCES KS.Body(BodyID)
)
ENGINE='Aria' TRANSACTIONAL=1
COMMENT "Weighting factors for bodies in collections.";
