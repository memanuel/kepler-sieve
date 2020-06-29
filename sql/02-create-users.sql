-- Create a user named kepler for the Kepler Sieve application
-- DROP USER IF EXISTS kepler@'%';
-- CREATE USER kepler@'%' IDENTIFIED BY 'kepler';
DROP USER IF EXISTS kepler;
CREATE USER kepler IDENTIFIED BY 'kepler';

-- Grant the user kepler global FILE privilege to load CSV files
GRANT FILE ON *.* TO 'kepler';

-- Grant the user Kepler on databases of the Kepler Sieve application
GRANT ALL PRIVILEGES ON KeplerDB.* TO kepler;
GRANT ALL PRIVILEGES ON JPL.* TO kepler;
GRANT DROP ON KeplerDB.* to kepler;
GRANT DROP ON JPL.* to kepler;
