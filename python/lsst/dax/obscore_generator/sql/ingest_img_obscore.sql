
insert into imgserv.obscore
(access_url,calib_level,core_id,dataproduct_subtype,dataproduct_type,ds_bg3_repo,facility_name,em_filter_name,
instrument_name,lastmodified,obs_collection,obs_id,obs_publisher_did,planeid,position_bounds_spoly,s_dec,
s_fov,s_pixel_scale,s_ra,s_region,s_xel1,s_xel2,t_exptime,t_max,t_min,target_name)
select access_url,CAST(calib_level AS integer),CAST(core_id AS uuid),dataproduct_subtype,dataproduct_type,
ds_bg3_repo,facility_name,em_filter_name,instrument_name,CAST(lastmodified AS timestamp),obs_collection,
obs_id,obs_publisher_did,CAST(planeid AS uuid),CAST(position_bounds_spoly AS spoly),CAST(s_dec AS float),
CAST(s_fov AS float),CAST(s_pixel_scale AS float),CAST(s_ra AS float),CAST(s_region AS varchar(512)),
CAST(s_xel1 AS integer),CAST(s_xel2 AS integer),CAST(t_exptime AS integer),CAST(t_max AS float),
CAST(t_min AS float),target_name
from imgserv.obscore_ci_hsc_raw

