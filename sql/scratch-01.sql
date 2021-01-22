# truncate TABLE KS.Integration_Planets;

LOAD DATA LOCAL INFILE 
'D:/Harvard/kepler-sieve/data/rebound/csv/Integration_Planets.csv'
IGNORE
INTO TABLE KS.Integration_Planets
FIELDS TERMINATED BY ','
LINES TERMINATED BY '\n'
IGNORE 1 LINES
(TimeID, BodyID, MJD, qx, qy, qz, vx, vy, vz);
