ALTER TABLE KS.OrbitalElements_DE435 DROP CONSTRAINT FK_OrbitalElements_DE435_TimeID;
ALTER TABLE KS.OrbitalElements_DE435 DROP CONSTRAINT FK_OrbitalElements_DE435_BodyID;
ALTER TABLE KS.OrbitalElements_DE435 engine='aria' transactional=0;
ALTER TABLE KS.OrbitalElements_DE435 ADD CONSTRAINT FK_OrbitalElements_DE435_TimeID FOREIGN KEY (TimeID) REFERENCES KS.IntegrationTime(TimeID);
ALTER TABLE KS.OrbitalElements_DE435 ADD CONSTRAINT FK_OrbitalElements_DE435_BodyID FOREIGN KEY (BodyID) REFERENCES KS.Body(BodyID);