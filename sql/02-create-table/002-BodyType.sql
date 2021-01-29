-- Drop tables to clear dependencies
DROP TABLE IF EXISTS JPL.LargeBody;
DROP TABLE IF EXISTS JPL.SmallBody;
DROP TABLE IF EXISTS KS.Body;

CREATE OR REPLACE TABLE KS.BodyType(
	BodyTypeID TINYINT NOT NULL PRIMARY KEY,
	BodyTypeCD VARCHAR(2) NOT NULL UNIQUE
		COMMENT "Short code",
	BodyTypeName VARCHAR(32)
		COMMENT "Full name",
	IsLargeBody_JPL BOOL NOT NULL
		COMMENT "Whether this type of body is considered 'large' or 'small' by JPL.",
	SortOrder TINYINT NOT NULL
)
	COMMENT "Types of Solar System bodies.";

INSERT INTO KS.BodyType
(BodyTypeID, BodyTypeCD, BodyTypeName, IsLargeBody_JPL, SortOrder)
VALUES
(0, 'SS', 'Solar System Barycenter', True, 1),
(1, 'S', 'Star', True, 1),
(2, 'PS', 'Planetary System Barycenter', True, 2),
(3, 'PB', 'Planet Single Body', True, 3),
(4, 'M', 'Moon', True, 4),
(5, 'A', 'Asteroid', False, 5);
