
delete from imgserv.obscore;

copy imgserv.obscore(access_url,calib_level,core_id,dataproduct_subtype,dataproduct_type,ds_bg3_repo,em_filter_name,em_max,em_min,facility_name,instrument_name,lastmodified,obs_collection,obs_id,obs_publisher_did,planeID,position_bounds_spoly,s_dec,s_fov,s_pixel_scale,s_ra,s_region,s_xel1,s_xel2,t_exptime,t_max,t_min,target_name)
from 'obscore_out.csv' delimiter ',' CSV HEADER;
