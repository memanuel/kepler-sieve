CREATE OR REPLACE TABLE ZTF.DetectionTime(
	DetectionTimeID INT NOT NULL AUTO_INCREMENT PRIMARY KEY
		COMMENT "Integer ID for this detection time",
    mjd DOUBLE NOT NULL UNIQUE
    	COMMENT "Time (mjd) of this detection"
)
ENGINE='Aria' TRANSACTIONAL=0
COMMENT "Distinct times at which ZTF detections are available.";
