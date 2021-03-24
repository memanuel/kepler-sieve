SELECT
	obj.oid AS ZtfObjectCD,
	obj.nobs AS ObservationCount,
	obj.mean_magap_g AS Mean_MagAp_g,
	obj.mean_magap_r AS Mean_MagAp_r,
	obj.mean_magpsf_g AS Mean_psf_g,
	obj.mean_magpsf_r AS Mean_psf_r,
	obj.meanra AS Mean_ra,
	obj.meandec AS Mean_dec,
	obj.firstmjd AS mjd0,
	obj.lastmjd AS mjd1
FROM
	public.objects AS obj
WHERE
	obj.classearly = 21
LIMIT 100;