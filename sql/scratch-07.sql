CREATE OR REPLACE TABLE KS.HiResTime(
	TimeID INT NOT NULL PRIMARY KEY
		COMMENT "MJD as integer number of minutes, e.g. floor(mjd*24*60)",
	mjd DOUBLE NOT NULL UNIQUE
		COMMENT "The Modified Julian Date in the TDB (barycentric dynamical time) frame",
	JD DOUBLE AS (mjd + 2400000.5)
		COMMENT "The Julian Date in the TDB (barycentric dynamical time) frame"
)
ENGINE='Aria' TRANSACTIONAL=1
COMMENT "Distinct time stamps at which MSE integrated positions of solar system bodies are available (every 5 minutes).";
