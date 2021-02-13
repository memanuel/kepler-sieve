ALTER TABLE KS.OrbitalElements_Planets DROP CONSTRAINT FK_OrbitalElements_Planets_TimeID;
ALTER TABLE KS.OrbitalElements_Planets DROP CONSTRAINT FK_OrbitalElements_Planets_BodyID;
ALTER TABLE KS.OrbitalElements_Planets engine='aria' transactional=0;
ALTER TABLE KS.OrbitalElements_Planets ADD CONSTRAINT FK_OrbitalElements_Planets_TimeID FOREIGN KEY (TimeID) REFERENCES KS.IntegrationTime(TimeID);
ALTER TABLE KS.OrbitalElements_Planets ADD CONSTRAINT FK_OrbitalElements_Planets_BodyID FOREIGN KEY (BodyID) REFERENCES KS.Body(BodyID);