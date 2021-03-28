-- ************************************************************************************************
-- Global variables used at bottom to build the tables
-- Set the variable N
SET @N = 1024;

-- Maximum distance for SkyPatchGridDistance
SET @dr_max = (0.5 / 360.0);

-- ************************************************************************************************
-- Build all the tables
-- CALL KS.MakeTable_SkyPatchDistance(@dr_max);
CALL KS.MakeTable_SkyPatchDistance_face(@dr_max);
-- CALL KS.MakeTable_SkyPatchDistance_extend(@dr_max);
-- CALL KS.MakeTable_SkyPatchDistance_min();