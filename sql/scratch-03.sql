SELECT
	gd.dr_mid,
	gd.dr_min
FROM KS.SkyPatchGridDistance AS gd
WHERE gd.dr_mid < gr_min;