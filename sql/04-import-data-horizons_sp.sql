DELIMITER //

create or replace procedure JPL.HorizonsImport_load_csv()
BEGIN 
	# load data infile "/home/Harvard/kepler-sieve/data/jpl/horizons/planets/010_sun.csv"
	#	into table JPL.HorizonsImport 
	#	fields terminated by ','
	#	lines terminated by '\n'
	#	ignore 1 lines
	#	(BodyNumber, BodyName, IntegrationSource, JD, CalendarDateTime, delta_T, qx, qy, qz, vx, vy, vz)
	#	set HorizonsImportID = NULL;
	select 1 as result;
END

DELIMITER ;