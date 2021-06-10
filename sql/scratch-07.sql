truncate TABLE JPL.AsteroidDirectionImport;

-- LOAD DATA INFILE '/ssd1/tmp/mysql/jpl/horizons/asteroid_directions/ast_geocenter.csv'
-- INTO TABLE JPL.AsteroidDirectionImport
-- FIELDS TERMINATED BY ","
-- LINES TERMINATED BY "\r\n"
-- IGNORE 1 LINES
-- (AsteroidID, ObservatoryID, JD, RA_ast, DEC_ast, RA_app, DEC_app, Mag, Brightness, r, rDot, delta, deltaDot, LightTime);
-- 
-- LOAD DATA INFILE '/ssd1/tmp/mysql/jpl/horizons/asteroid_directions/ast_palomar.csv'
INTO TABLE JPL.AsteroidDirectionImport
FIELDS TERMINATED BY ","
LINES TERMINATED BY "\r\n"
IGNORE 1 LINES
(AsteroidID, ObservatoryID, JD, RA_ast, DEC_ast, RA_app, DEC_app, Mag, Brightness, r, rDot, delta, deltaDot, LightTime);
