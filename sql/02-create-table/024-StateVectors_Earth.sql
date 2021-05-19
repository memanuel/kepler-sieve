-- State vectors from MSE integration using Rebound
CREATE OR REPLACE TABLE KS.StateVectors_Earth(
	TimeID INT NOT NULL
		COMMENT "Integer ID for the timestamp of these state vectors; FK to KS.HiResTime",
	mjd DOUBLE NOT NULL
		COMMENT "The Modified Julian Date in the TDB (barycentric dynamical time) frame; derivable from IntegrationTimeID but included for performance.",
	-- Position of Earth q = [qx, qy, qz]
	qx DOUBLE NOT NULL
		COMMENT "Position of Earth (x coordinate) in AU in the barcycentric mean ecliptic frame",
	qy DOUBLE NOT NULL
		COMMENT "Position of Earth (y coordinate) in AU in the barcycentric mean ecliptic frame",
	qz DOUBLE NOT NULL
		COMMENT "Position of Earth (z coordinate) in AU in the barcycentric mean ecliptic frame",
	-- Velocity of Earth v = [vx, vy, vz]
	vx DOUBLE NOT NULL
		COMMENT "Velocity of Earth (x coordinate) in AU/day in the barcycentric mean ecliptic frame",
	vy DOUBLE NOT NULL
		COMMENT "Velocity of Earth (y coordinate) in AU/day in the barcycentric mean ecliptic frame",
	vz DOUBLE NOT NULL
		COMMENT "Velocity of Earth (z coordinate) in AU/day in the barcycentric mean ecliptic frame",
	-- Keys and constraints
	PRIMARY KEY (TimeID)
		COMMENT "An Earth state vector is identified by the time stamp; use integer time ID for performance.",
	CONSTRAINT FK_StateVectors_Earth_TimeID
		FOREIGN KEY (TimeID) REFERENCES KS.HiResTime(TimeID)
)
ENGINE='Aria' TRANSACTIONAL=0
COMMENT "State vectors (position and velocity) for Earth computed in rebound using the planets as massive bodies and initial conditions from DE435 at MJD 59000."
