use JPL;

# load data local infile "d:/Harvard/kepler-sieve/data/jpl/horizons/planets/010_sun.csv" into table HorizonsImport;
# load data infile "/home/michael/Harvard/kepler-sieve/data/jpl/horizons/planets/010_sun.csv" into table HorizonsImport; 

load data infile "/ssd1/tmp/mysql/010_sun.csv" 
into table HorizonsImport 
fields terminated by ','
(BodyNumber,BodyName,IntegrationSource,JD,CalendarDateTime,delta_T,qx,qy,qz,vx,vy,vz)
# lines terminated by '\n'
ignore 1 lines;
# set HorizonsImportID = NULL
