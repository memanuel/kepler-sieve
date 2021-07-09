use master;
go

CREATE DATABASE KeplerDB
ON  PRIMARY 
(NAME = N'KeplerDB_dat', FILENAME = N'/ssd1/mssql/dbfiles/KeplerDB_dat.mdf', MAXSIZE = UNLIMITED),
( NAME = N'KeplerDB_dat2', FILENAME = N'/hdd1/mssql/dbfiles/KeplerDB_dat2.ndf' , MAXSIZE = UNLIMITED), 
FILEGROUP [MemoryOptimizedData] CONTAINS MEMORY_OPTIMIZED_DATA  DEFAULT
( NAME = N'MemOptData', FILENAME = N'/ssd1/mssql/dbfiles/MemOptData/KeplerDB' , MAXSIZE = UNLIMITED)
LOG ON 
( NAME = N'KeplerDB_log', FILENAME = N'/hdd1/mssql/tlogfiles/KeplerDB_log.ldf'), 
( NAME = N'KeplerDB_log2', FILENAME = N'/hdd1/mssql/tlogfiles/KeplerDB_log2.ldf')
WITH CATALOG_COLLATION = DATABASE_DEFAULT[dbo].[Color]
GO
