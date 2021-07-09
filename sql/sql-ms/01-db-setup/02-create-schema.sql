USE KeplerDB
GO

-- Get rid of default schema dbo
DROP SCHEMA dbo;
GO

-- Main schema
CREATE SCHEMA KS
-- AUTHORIZATION kepler
GO

-- Data imported from JPL
CREATE SCHEMA JPL
-- AUTHORIZATION kepler
GO

