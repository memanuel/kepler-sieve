CREATE OR REPLACE TABLE ZTF.Detection(
	DetectionID BIGINT NOT NULL PRIMARY KEY
		COMMENT "Integer ID for this asteroid detection",
    ObjectCD VARCHAR(12) NOT NULL
    	COMMENT "String ID for the object this detection was associated with by Alerce",
    mjd DOUBLE NOT NULL
    	COMMENT "Time (mjd) of this detection",
   	ra DOUBLE NOT NULL
	   	COMMENT "Right Ascension (RA) of this detection",
   	`dec` DOUBLE NOT NULL
	   	COMMENT "Declination (DEC) of this detection",
	MagPSF DOUBLE NOT NULL
		COMMENT "Magnitude of the PSF of this detection",
	MagApp DOUBLE NOT NULL
		COMMENT "Apparent Magnitude of this detection",
	MagNR DOUBLE NOT NULL
		COMMENT "Some other magnitude of this detection, not sure exactly what",
	Sigma_RA DOUBLE NOT NULL
		COMMENT "Standard devation of RA estimate for this detection",
	Sigma_DEC DOUBLE NOT NULL
		COMMENT "Standard deviation of DEC estimate for this detection"
   	-- Constraints
   	-- CONSTRAINT FK_Detection_ObjectCD FOREIGN KEY (ObjectCD) REFERENCES ZTF.Object(ObjectCD)
)
ENGINE='Aria' TRANSACTIONAL=0
COMMENT "Asteroid detections in ZTF2 data according to Alerce.";
