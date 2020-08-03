CREATE OR REPLACE TABLE JPL.AsteroidElement_Numbered(
	AsteroidNumber INT NOT NULL PRIMARY KEY
		COMMENT "The IAU asteroid number",
	AsteroidName VARCHAR(32) NOT NULL UNIQUE
		COMMENT "The IAU asteroid name",
	Epoch DOUBLE NOT NULL
		COMMENT "The epoch as of which the orbital elements are computed; expressed as an MJD",
	a DOUBLE NOT NULL
		COMMENT "The semimajor axis in AU",
	e DOUBLE NOT NULL
		COMMENT "The eccentricity; dimensionless",
	i DOUBLE NOT NULL
		COMMENT "The inclination in degrees",
	w DOUBLE NOT NULL
		COMMENT "The argument of perhelion omega in degrees",
	Node DOUBLE NOT NULL
		COMMENT "The longitude of the ascending node in degrees",
	M DOUBLE NOT NULL
		COMMENT "The mean anomaly M in degrees",
	H DOUBLE NOT NULL
		COMMENT "The H brightness parameter",
	G DOUBLE NOT NULL
		COMMENT "The G brightness parameter",
	`Ref` VARCHAR(32)
		COMMENT "JPL description of the integration used to calculate this ephemeris"
)
COMMENT "Import JPL file with orbital elements of numbered asteroids. Primary is the Sun, NOT Solar System Barycenter!";
