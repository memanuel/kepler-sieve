CREATE OR REPLACE TABLE KS.AsteroidElement_Ref(
	AsteroidID INT NOT NULL 
		COMMENT "Foreign key to Asteroid table",
	TimeID INT NOT NULL
		COMMENT "Foreign key to IntegrationTime table",
	epoch DOUBLE NOT NULL
		COMMENT "The epoch as of which the orbital elements are computed; expressed as an MJD",
	-- The six main orbital elements plus mean anomaly for interpolation
	a DOUBLE NOT NULL
		COMMENT "The semimajor axis, a, in AU",
	e DOUBLE NOT NULL
		COMMENT "The eccentricity, e; dimensionless",
	inc DOUBLE NOT NULL
		COMMENT "The inclination, inc, in radians",
	Omega_node DOUBLE NOT NULL
		COMMENT "The longitude of the ascending node, Omega, in radians",
	omega_peri DOUBLE NOT NULL
		COMMENT "The argument of perhelion, omega, in radians",
    f DOUBLE NOT NULL
        COMMENT "The true anomaly, f, in radians",
	M DOUBLE NOT NULL
		COMMENT "The mean anomaly, M, in radians",
    -- Physics calculations
	d DOUBLE NOT NULL
		COMMENT "The radial distance from the primary in AU",
	v DOUBLE NOT NULL
		COMMENT "The speed relative to the primary in AU/day",
	h DOUBLE NOT NULL
		COMMENT "The specific angular momentum",
	-- Period and motion
    period DOUBLE NOT NULL
        COMMENT "The orbital period, P, in days",
    mean_motion DOUBLE NOT NULL
        COMMENT "The mean motion, n, in radians per day",
    T_peri DOUBLE NOT NULL
    	COMMENT "The time of pericenter passage",
	-- Additional angles
	pomega DOUBLE NOT NULL
		COMMENT "longitude of pericenter in radians",
	EA DOUBLE AS (MOD(2.0*ATAN(SQRT((1.0-e)/(1.0+e))*TAN(0.5*f))+2.0*PI(), 2.0*PI())) PERSISTENT
		COMMENT "The eccentric anomaly; derived from the true anomaly",		
    mean_longitude DOUBLE AS (MOD(Omega_node + omega_peri + M, 2.0*PI())) PERSISTENT
        COMMENT "Mean longitude, l, in radians = Omega + omega + M",        
    true_longitude DOUBLE AS (MOD(Omega_node + omega_peri + f, 2.0*PI())) PERSISTENT
        COMMENT "The angle, theta, in radians = Omega + omega + f",
	-- Keys and constraints
    PRIMARY KEY (AsteroidID, TimeID)
        COMMENT "Reference elements uniquely defined by the asteroid and the epoch",
    UNIQUE KEY (TimeID, AsteroidID)
        COMMENT "Support fast searching by time ID",
    CONSTRAINT FK_AsteroidElement_Ref_AsteroidID
        FOREIGN KEY (AsteroidID) REFERENCES KS.Asteroid(AsteroidID),
    CONSTRAINT FK_AsteroidElement_Ref_TimeID
        FOREIGN KEY (TimeID) REFERENCES KS.IntegrationTime(TimeID)
)
ENGINE='Aria' TRANSACTIONAL=0
COMMENT "Orbital elements of asteroids as of reference dates; used to initialize integrations. Primary for these elements is the Sun, NOT the Solar System Barycenter!  See https://rebound.readthedocs.io/en/latest/python_api.html for details on rebound orbital elements";
