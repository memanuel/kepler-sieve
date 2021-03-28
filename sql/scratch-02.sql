SELECT * FROM KS.SkyPatchDistance_cand;
SELECT * FROM KS.SkyPatchDistance_cand ORDER BY dr_mid desc;

SELECT max(dr_mid) AS dr_max FROM KS.SkyPatchDistance_cand;
SELECT count(*) AS dr_max FROM KS.SkyPatchDistance_cand AS cand_count;

SELECT max (dr_mid) AS dr_max FROM KS.SkyPatchGridDistance;

SELECT 0.0041431949*360.0 AS dr_max_deg;

DELETE FROM KS.SkyPatchGridDistance WHERE dr_min > (1.0 / 360.0);