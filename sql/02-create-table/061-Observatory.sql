CREATE OR REPLACE TABLE KS.Observatory(
	ObservatoryID TINYINT NOT NULL PRIMARY KEY,
    ObservatoryCD CHAR(3) NOT NULL UNIQUE
        COMMENT "Minor Planet Center (MPC) code for this observatory",
	ObservatoryShortName VARCHAR(16) NOT NULL UNIQUE
		COMMENT "Name of this observatory",
	ObservatoryName VARCHAR(32) NOT NULL UNIQUE
		COMMENT "Long name of this observatory",
    Latitude DOUBLE NOT NULL
        COMMENT "Earth latitude as a double; positive is North, negative is South",
    Longitude DOUBLE NOT NULL
        COMMENT "Earth longitude as a double; positive is East, negative is West ",
    Height DOUBLE NOT NULL
        COMMENT "Elevation of this site from mean sea level in meters",
    qx DOUBLE NOT NULL
        COMMENT "x coordinate according to Skyfield on IERS2010 Topos model",
    qy DOUBLE NOT NULL
        COMMENT "y coordinate according to Skyfield on IERS2010 Topos model",
    qz DOUBLE NOT NULL
        COMMENT "z coordinate according to Skyfield on IERS2010 Topos model",
    SkyfieldSiteName VARCHAR(16) NOT NULL
        COMMENT "Name of this site in Skyfield library",
	SortOrder TINYINT NOT NULL
)
ENGINE='Aria' TRANSACTIONAL=1
COMMENT "Astronomical observatories and their position on earth.";

-- See https://en.wikipedia.org/wiki/List_of_observatory_codes for Minor Planet Center
-- list of observatory codes

-- See https://en.wikipedia.org/wiki/Palomar_Observatory
-- https://geohack.toolforge.org/geohack.php?pagename=Palomar_Observatory&params=33.3564_N_116.865_W_
INSERT INTO KS.Observatory
(ObservatoryID, ObservatoryCD, ObservatoryShortName, ObservatoryName, 
 Latitude, Longitude, Height, qx, qy, qz, SkyfieldSiteName, SortOrder)
VALUES
(0, '___', 'Geocenter', 'Earth Geocenter', 
 90.0, 0.0, -6356752.31424518, 0.0, 0.0, 0.0, 'geocenter', 0),
(1, 'I41', 'ZTF', 'Zwicky Transient Facility',
 33.3564, -116.865, 1712.0, -2410503.9835643, -4758565.10396658, 3487983.06777811, 'palomar', 1);
