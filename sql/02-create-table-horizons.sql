-- Create tables for Horizons data import
use JPL;

create or replace table HorizonsImport(
	HorizonsImportID int not null auto_increment primary key,
	BodyNumber int unsigned not null
		comment 'JPL number of this body; separate series for small and large bodies',
	BodyName varchar(32) not null
		comment 'JPL name of this body',
	IntegrationSource varchar(32) not null
		comment 'JPL description of the integration used to calculate this ephemeris',
	JD double not null
		comment 'The Julian Date in the TDB (barycentric dynamical time) frame',
	CalendarDateTime datetime not null
		comment 'The date and time on the Gregorian calendar in the TDB frame',
	delta_T double not null
		comment 'The difference between the TDB and terrestial atomic time frames',
	qx double not null
		comment 'Position of body (x coordinate) in AU in the barcycentric mean ecliptic frame',
	qy double not null
		comment 'Position of body (y coordinate) in AU in the barcycentric mean ecliptic frame',
	qz double not null
		comment 'Position of body (z coordinate) in AU in the barcycentric mean ecliptic frame',
	vx double not null
		comment 'Velocity of body (x coordinate) in AU/day in the barcycentric mean ecliptic frame',
	vy double not null
		comment 'Velocity of body (y coordinate) in AU/day in the barcycentric mean ecliptic frame',
	vz double not null
		comment 'Velocity of body (z coordinate) in AU/day in the barcycentric mean ecliptic frame'
)
