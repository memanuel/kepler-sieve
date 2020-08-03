CREATE OR REPLACE TABLE JPL.AsteroidElement(
	AsteroidNumber INT NOT NULL PRIMARY KEY
		COMMENT "The IAU asteroid number when it exists; otherwise a sequential counter starting at 1,000,000",
	AsteroidName VARCHAR(32) NOT NULL UNIQUE
		COMMENT "The IAU asteroid name where applicable or asteroid designation",
    IsNumberedAsteroid BOOL NOT NULL
        COMMENT "Flag indicating whether or not this is an IAU numbered asteroid",
	epoch DOUBLE NOT NULL
		COMMENT "The epoch as of which the orbital elements are computed; expressed as an MJD",
	a DOUBLE NOT NULL
		COMMENT "The semimajor axis in AU",
	e DOUBLE NOT NULL
		COMMENT "The eccentricity; dimensionless",
	inc DOUBLE NOT NULL
		COMMENT "The inclination in radians",
	Omega_node DOUBLE NOT NULL
		COMMENT "The longitude of the ascending node in radians",
	omega_peri DOUBLE NOT NULL
		COMMENT "The argument of perhelion omega in radians",
    f DOUBLE NOT NULL
        COMMENT "The true anomaly f in radians",
	M DOUBLE NOT NULL
		COMMENT "The mean anomaly M in radians",
    eccentric_anomaly DOUBLE NOT NULL
        COMMENT "The eccentric anomaly E in radians",
    period DOUBLE NOT NULL
        COMMENT "The orbital period in days",
    mean_motion DOUBLE NOT NULL
        COMMENT "The mean motion in radians per day",
	H DOUBLE NOT NULL
		COMMENT "The H brightness parameter",
	G DOUBLE NOT NULL
		COMMENT "The G brightness parameter",
	`Ref` VARCHAR(32)
		COMMENT "JPL description of the integration used to calculate this ephemeris",
	`row_num` int NOT NULL
)
COMMENT "Import JPL files with orbital elements of both numbered and unnumbered asteroids. Primary is the Sun, NOT Solar System Barycenter! Performs calculations on quoted elements to convert them to radians and add true anomaly f from mean anomly M.";
