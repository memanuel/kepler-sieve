-- ************************************************************************************************
-- Global variables to build the tables
-- Set the grid size N
SET @N = 1024;

-- Maximum distance for SkyPatchGridDistance and SkyPatchDistance
SET @dr_max = (0.5 / 360.0);

-- ************************************************************************************************
-- Build all the tables in the SkyPatch distance family
CALL KS.MakeTable_SkyPatchGrid(N, dr_max);
CALL KS.MakeTable_SkyPatch(@dr_max);

-- Individual stages of building SkyPatchDistance for testing / re-running part of job
-- CALL KS.MakeTable_SkyPatchDistance(@dr_max);
-- CALL KS.MakeTable_SkyPatchDistance_extend(@dr_max);
-- CALL KS.MakeTable_SkyPatchDistance_min();
