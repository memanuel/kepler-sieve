CREATE OR REPLACE TABLE KS.IntegrationTime(
	TimeID INT NOT NULL PRIMARY KEY
		COMMENT "MJD as integer number of minutes, e.g. floor(mjd*24*60)",
	mjd DOUBLE NOT NULL UNIQUE
		COMMENT "The Modified Julian Date in the TDB (barycentric dynamical time) frame",
	JD DOUBLE AS (mjd + 2400000.5)
		COMMENT "The Julian Date in the TDB (barycentric dynamical time) frame",
	CalendarDate DATE NOT NULL
		COMMENT "The date and time on the Gregorian calendar in the TDB frame",
	CalendarDateTime DATETIME(6) NOT NULL UNIQUE
		COMMENT "The date and time on the Gregorian calendar in the TDB frame",
	delta_T DOUBLE NOT NULL
		COMMENT "The difference between the TDB and terrestial atomic time frames"
)
ENGINE='Aria' TRANSACTIONAL=1
COMMENT "Distinct time stamps at which MSE integrated positions of solar system bodies are available.";
