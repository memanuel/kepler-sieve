SELECT
	COUNT(spgd.dr_min) as num,
	MIN(spgd.dr_min*360) AS dr_min,
	AVG(spgd.dr_min*360) AS dr_avg,
	MAX(spgd.dr_min*360) AS dr_max
FROM KS.SkyPatchGridDistance AS spgd
where spgd.dr_min > 0;

SELECT
	COUNT(spd.dr_min) as num,
	MIN(spd.dr_min*360) AS dr_min,
	AVG(spd.dr_min*360) AS dr_avg,
	MAX(spd.dr_min*360) AS dr_max
FROM KS.SkyPatchDistance AS spd
where spd.dr_min > 0;

select * from KS.SkyPatchGridDistance;
select * from KS.SkyPatchDistance;


select * from KS.Corner;