-- ************************************************************************************************
-- Neighbor distance on the SkyPatchGrid table
CREATE OR REPLACE TABLE KS.SkyPatchGridNeighbor(
	i1 INT NOT NULL
        COMMENT "i coordinate of the first SkyPatchGrid cell",
	j1 INT NOT NULL
        COMMENT "j coordinate of the first SkyPatchGrid cell",
	i2 INT NOT NULL
        COMMENT "i coordinate of the second SkyPatchGrid cell",
	j2 INT NOT NULL
        COMMENT "ji coordinate of the second SkyPatchGrid cell",
	dr DOUBLE NOT NULL
		COMMENT "Distance from center point to center point",
	dr_min DOUBLE NOT NULL
		COMMENT "Estimated minimum distance based on getting closer by up to half the width on each end",
	PRIMARY KEY (i1, j1, i2, j2)
)
ENGINE='Aria' TRANSACTIONAL=0
COMMENT "Catalog neighbor interactions on SkyPatchGrid table.";


-- ************************************************************************************************
-- Neighbor distance on the SkyPatch table
CREATE OR REPLACE TABLE KS.SkyPatchNeighbor(
	SkyPatchID_1 INT NOT NULL
        COMMENT "The first SkyPatchID",
	SkyPatchID_2 INT NOT NULL
        COMMENT "The second SkyPatchID",
	dr DOUBLE NOT NULL
		COMMENT "Distance from center point to center point",
	dr_min DOUBLE NOT NULL
		COMMENT "Estimated minimum distance based on getting closer by up to half the width on each end",
	PRIMARY KEY (SkyPatchID_1, SkyPatchID_2)
)
ENGINE='Aria' TRANSACTIONAL=0
COMMENT "Catalog neighbor interactions on SkyPatch table.";
