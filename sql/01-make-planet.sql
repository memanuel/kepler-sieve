use KeplerDB;
go;

drop table if exists Planet;
create table Planet(
PlanetID smallint not null,
PlanetName varchar(64) not null,
Mass_Solar double not null,
Mass_KG double as Mass_Solar * 
constraint PK_Planet_PlanetID primary key (PlanetID),
constraint UNQ_Planet_PlanetName unique (PlanetName),
);