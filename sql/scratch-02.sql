SELECT * FROM KS.Integration_Planets_chunk_000;

INSERT IGNORE INTO KS.Integration_Planets
(TimeID, BodyID, MJD, qz, qy, qx, vx, vy, vz)
SELECT
TimeID, BodyID, MJD, qz, qy, qx, vx, vy, vz
FROM KS.Integration_Planets_chunk_000;