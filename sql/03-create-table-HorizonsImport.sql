-- Create tables for Horizons data import
CREATE OR REPLACE TABLE JPL.HorizonsImport(
	HorizonsImportID INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
	BodyTypeCD VARCHAR(2) NOT NULL
		COMMENT "The type of this body.  One of 'S' (star), 'PS' (planet system barycenter), 'PB', (planet single body) 'M' (moon), 'A' (asteroid)",
	BodyNumber INT NOT NULL
		COMMENT "JPL number of this body; separate series for small and large bodies",
	BodyName VARCHAR(32) NOT NULL
		COMMENT "JPL name of this body",
	IntegrationSource VARCHAR(32) NOT NULL
		COMMENT "JPL description of the integration used to calculate this ephemeris",
	JD DOUBLE NOT NULL
		COMMENT "The Julian Date in the TDB (barycentric dynamical time) frame",
	CalendarDateTime DATETIME(6) NOT NULL
		COMMENT "The date and time on the Gregorian calendar in the TDB frame",
	delta_T DOUBLE NOT NULL
		COMMENT "The difference between the TDB and terrestial atomic time frames",
	qx DOUBLE NOT NULL
		COMMENT "Position of body (x coordinate) in AU in the barcycentric mean ecliptic frame",
	qy DOUBLE NOT NULL
		COMMENT "Position of body (y coordinate) in AU in the barcycentric mean ecliptic frame",
	qz DOUBLE NOT NULL
		COMMENT "Position of body (z coordinate) in AU in the barcycentric mean ecliptic frame",
	vx DOUBLE NOT NULL
		COMMENT "Velocity of body (x coordinate) in AU/day in the barcycentric mean ecliptic frame",
	vy DOUBLE NOT NULL
		COMMENT "Velocity of body (y coordinate) in AU/day in the barcycentric mean ecliptic frame",
	vz DOUBLE NOT NULL
		COMMENT "Velocity of body (z coordinate) in AU/day in the barcycentric mean ecliptic frame",
	KEY BodyTypeCD_BodyNumber (BodyTypeCD, BodyNumber, JD)
		COMMENT "Don't enforce uniqueness here to allow loading e.g. daily and weekly data files without a collision."
)
	COMMENT "Staging table to import data from files downloaded from Horizons API."
