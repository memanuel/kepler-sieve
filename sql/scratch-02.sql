SET @N = 16;
SET @dr_max = GREATEST(1.0/@N, 0.5 / 360.0);

CALL KS.MakeTable_SkyPatchGrid(@N, @dr_max);

SELECT * FROM KS.SkyPatchGridDistance WHERE dr_mid < dr_min;
