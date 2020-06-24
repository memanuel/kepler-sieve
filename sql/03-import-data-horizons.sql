truncate table JPL.HorizonsImport;

# load data infile "/ssd1/tmp/mysql/test.csv"
load data infile "/home/Harvard/kepler-sieve/data/jpl/horizons/planets/010_sun.csv"
into table JPL.HorizonsImport 
fields terminated by ','
lines terminated by '\n'
ignore 1 lines
(BodyNumber, BodyName, IntegrationSource, JD, CalendarDateTime, delta_T, qx, qy, qz, vx, vy, vz)
set HorizonsImportID = NULL;