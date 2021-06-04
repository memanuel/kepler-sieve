SELECT * FROM JPL.HorizonsVectors WHERE HorizonsBodyID = 1000001;

SELECT 
min(mjd) AS mjd0,
max(mjd) AS mjd1
FROM JPL.HorizonsVectors 
WHERE HorizonsBodyID = 1000001;