
drop table if exists imgserv.obscore;

create table imgserv.obscore
(
    dataproduct_type varchar(64),
    dataproduct_subtype varchar(64),
    calib_level integer,
    obs_collection varchar(64) not null,
    facility_name varchar(64),
    instrument_name varchar(64),
    obs_id varchar(256) not null,
    obs_publisher_did varchar(512) not null,
    access_url varchar(512) not null,
    access_format varchar default 'image/fits',
    s_ra double precision not null,
    s_dec double precision not null,
    s_pixel_scale double precision not null,
    o_ucd varchar default 'phot.count',
    core_id uuid not null primary key,
    lastmodified timestamp not null,
    s_region varchar(512) not null,
    position_bounds_spoly spoly,
    planeid uuid not null,
    -- things to add
    ds_bg3_repo varchar(512) not null,
    em_filter_name varchar(10),
    access_estsize integer,
    target_name varchar(64),
    s_fov double precision not null,
    s_resolution double precision,
    -- number of pixels in horizontal direction
    s_xel1 integer,
    -- number of pixels in vertical direction
    s_xel2 integer,
    -- MJD for exposure start
    t_min numeric(12,6),
    -- MJD for exposure end
    t_max numeric(12,6),
    t_exptime integer,
    t_resolution integer,
    t_xel integer default 1,
    em_min double precision not null,
    em_max double precision not null,
    em_res_power double precision,
    em_xel integer,
    pol_states integer,
    pol_xel integer default 0,
    em_ucd varchar default 'em.wl',
    em_unit varchar default 'm'
)
;
