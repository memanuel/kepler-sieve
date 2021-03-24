CREATE OR REPLACE TABLE ZTF.Object(
    ObjectID INT NOT NULL AUTO_INCREMENT PRIMARY KEY
    	COMMENT "Integer ID for this object; assigned by MSE",
    ObjectCD VARCHAR(12) NOT NULL
    	COMMENT "String ID for this object, assigned by Alerce",
    ObservationCount INT NOT NULL
    	COMMENT "Number of observations associated with this object",
    mjd0 DOUBLE NOT NULL
    	COMMENT "Time (mjd) this object was first detected",
    mjd1 DOUBLE NOT NULL
    	COMMENT "Time (mjd) this object was last detected",
	MeanRA DOUBLE NOT NULL
		COMMENT "Mean right ascension (RA) of the detections of this object",
	MeanDEC DOUBLE NOT NULL
		COMMENT "Mean declination (DEC) of the detections of this object",
	MeanMag_g DOUBLE NOT NULL
		COMMENT "Mean magnitude in green band of the dections of this object",
	MeanMag_r DOUBLE NOT NULL
		COMMENT "Mean magnitude in red band of the dections of this object",
	MeanMagPSF_g DOUBLE NOT NULL
		COMMENT "Mean PSF (point spread function) magnitude in green band of the dections of this object",
	MeanMagPSF_r DOUBLE NOT NULL
		COMMENT "Mean PSF (point spread function) magnitude in red band of the dections of this object",
	ObjectClassID TINYINT NOT NULL
		COMMENT "The ObectClassID that Alerce classified this object as, e.g. 21 for Asteroid",
	ClassificationProb DOUBLE NOT NULL
   		COMMENT "Alerce classification probability that this object was in the relevant class",
-- Keys
UNIQUE KEY UNQ_ObjectCD (ObjectCD)
	COMMENT "The string code for each object is unique"
)
ENGINE='Aria' TRANSACTIONAL=0
COMMENT "Objects detected by Alerce in ZTF2 data.";
