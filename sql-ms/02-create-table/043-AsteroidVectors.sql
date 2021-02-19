Use KeplerDB;

-- Create tables for integrated vectors for asteroids
DROP TABLE IF EXISTS KS.AsteroidVectors;

CREATE TABLE KS.AsteroidVectors(
	TimeID INT NOT NULL, 
	AsteroidID INT NOT NULL, 
	MJD DOUBLE PRECISION NOT NULL, 
	-- Position q = [qx, qy, qz]
	qx DOUBLE PRECISION NOT NULL,
	qy DOUBLE PRECISION NOT NULL,
	qz DOUBLE PRECISION NOT NULL,
	-- Velocity v = [vx, vy, vz]
	vx DOUBLE PRECISION NOT NULL,
	vy DOUBLE PRECISION NOT NULL,
	vz DOUBLE PRECISION NOT NULL,
	-- Primary key and unique constraints
	CONSTRAINT PK_AsteroidVectors PRIMARY KEY (TimeID, AsteroidID),
	CONSTRAINT UNQ_AsteroidID_TimeID UNIQUE (AsteroidID, TimeID)
	-- Foreign keys
	-- CONSTRAINT FK_AsteroidVectors_TimeID FOREIGN KEY (TimeID) REFERENCES KS.DailyTime(TimeID)
	-- CONSTRAINT FK_AsteroidVectors_BodyID	FOREIGN KEY (AsteroidID) REFERENCES KS.Asteroid(AsteroidID)
)

-- Table comments
exec sys.sp_addextendedproperty 
@name=MS_Description,
@level0type = 'Schema', @level0name='KS',
@level1type = 'Table',  @level1name='AsteroidVectors',
@value=
"State vectors (position and velocity) for asteroids computed in Rebound using the planets as massive bodies";

exec sys.sp_addextendedproperty 
@name=MS_Description,
@level0type = 'Schema',  @level0name='KS',
@level1type = 'Table',   @level1name='AsteroidVectors',
@level2type = 'Column',  @level2name='TimeID',
@value=
"Integer ID for the timestamp of these state vectors; FK to KS.IntegrationTime";

exec sys.sp_addextendedproperty 
@name=MS_Description ,
@level0type = 'Schema',  @level0name='KS',
@level1type = 'Table',   @level1name='AsteroidVectors',
@level2type = 'Column',  @level2name='AsteroidID',
@value=
"The asteroid whose state vectors are described; FK to KS.Asteroid";

exec sys.sp_addextendedproperty 
@name=MS_Description ,
@level0type = 'Schema',  @level0name='KS',
@level1type = 'Table',   @level1name='AsteroidVectors',
@level2type = 'Column',  @level2name='MJD',
@value=
"The Modified Julian Date in the TDB frame; derivable from TimeID but included for performance.";

exec sys.sp_addextendedproperty 
@name=MS_Description ,
@level0type = 'Schema',  @level0name='KS',
@level1type = 'Table',   @level1name='AsteroidVectors',
@level2type = 'Column',  @level2name='qx',
@value=
"Position of body (x coordinate) in AU in the barcycentric mean ecliptic frame";

exec sys.sp_addextendedproperty 
@name=MS_Description ,
@level0type = 'Schema',  @level0name='KS',
@level1type = 'Table',   @level1name='AsteroidVectors',
@level2type = 'Column',  @level2name='qy',
@value=
"Position of body (y coordinate) in AU in the barcycentric mean ecliptic frame";

exec sys.sp_addextendedproperty 
@name=MS_Description ,
@level0type = 'Schema',  @level0name='KS',
@level1type = 'Table',   @level1name='AsteroidVectors',
@level2type = 'Column',  @level2name='qz',
@value=
"Position of body (z coordinate) in AU in the barcycentric mean ecliptic frame";

exec sys.sp_addextendedproperty 
@name=MS_Description ,
@level0type = 'Schema',  @level0name='KS',
@level1type = 'Table',   @level1name='AsteroidVectors',
@level2type = 'Column',  @level2name='vx',
@value=
"Velocity of body (x coordinate) in AU/day in the barcycentric mean ecliptic frame";

exec sys.sp_addextendedproperty 
@name=MS_Description ,
@level0type = 'Schema',  @level0name='KS',
@level1type = 'Table',   @level1name='AsteroidVectors',
@level2type = 'Column',  @level2name='vy',
@value=
"Velocity of body (y coordinate) in AU/day in the barcycentric mean ecliptic frame";

exec sys.sp_addextendedproperty 
@name=MS_Description ,
@level0type = 'Schema',  @level0name='KS',
@level1type = 'Table',   @level1name='AsteroidVectors',
@level2type = 'Column',  @level2name='vz',
@value=
"Velocity of body (z coordinate) in AU/day in the barcycentric mean ecliptic frame";

-- Comments on constraints
exec sys.sp_addextendedproperty 
@name=MS_Description ,
@level0type = 'Schema',		 @level0name='KS',
@level1type = 'Table',		 @level1name='AsteroidVectors',
@level2type = 'Constraint',  @level2name='PK_AsteroidVectors',
@value=
"A state vector is identified by the body and time stamp; use integer time ID for performance."

exec sys.sp_addextendedproperty 
@name=MS_Description,
@level0type = 'Schema',		 @level0name='KS',
@level1type = 'Table',		 @level1name='AsteroidVectors',
@level2type = 'Constraint',  @level2name='UNQ_AsteroidID_TimeID',
@value=
"Allow fast search keyed first by AsteroidID.";