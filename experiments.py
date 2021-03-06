

phoz_LSSTgold = {
    "type" : "gal_clustering",
    "name" : "LSST_gold",
    "nzfi" : "curves_LSST/nz_shear_fiducial.txt",
    "bzfi" : "curves_LSST/bz_gold.txt",
    "szfi" : "curves_LSST/sz_gold.txt",
    "ezfi" : "curves_LSST/ez_gold.txt",
    "fsky" : 0.4
    }

spec_DESI = {
    "type" : "gal_clustering",
    "name" : "DESI",
    "nzfi" : "curves_spec/nz_DESI.txt",
    "bzfi" : "curves_spec/bz_DESI.txt",
    "szfi" : "curves_spec/sz_DESI.txt",
    "ezfi" : "curves_spec/ez_DESI.txt",
    "fsky" : 0.2
    }

spec_Euclid = {
    "type" : "gal_clustering",
    "name" : "Euclid",
    "nzfi" : "curves_spec/nz_Euclid.txt",
    "bzfi" : "curves_spec/bz_Euclid.txt",
    "szfi" : "curves_spec/sz_Euclid.txt",
    "ezfi" : "curves_spec/ez_Euclid.txt",
    "fsky" : 0.4
    }

spec_WFIRST = {
    "type" : "gal_clustering",
    "name" : "WFIRST",
    "nzfi" : "curves_spec/nz_WFIRST.txt",
    "bzfi" : "curves_spec/bz_WFIRST.txt",
    "szfi" : "curves_spec/sz_WFIRST.txt",
    "ezfi" : "curves_spec/ez_WFIRST.txt",
    "fsky" : 0.05
    }

im_HIRAX_32_6 = {
    "type" : "intensity_mapping",
    "name" : "HIRAX_32_6",
    "nzfi" : "curves_IM/nz_HI.txt",
    "bzfi" : "curves_IM/bz_HI.txt",
    "szfi" : "curves_IM/sz_HI.txt",
    "ezfi" : "curves_IM/ez_HI.txt",
    "tzfi" : "curves_IM/tz_HI.txt",
    "dish_size" : 6.,
    "t_inst" : 50.,
    "t_total" : 10000.,
    "n_dish" : 1024,
    "area_eff" : 1.0,
    "im_type" : "interferometer",
    "base_file" : "curves_IM/baseline_file_HIRAX_6m.txt",
    "base_min" : 6.,
    "base_max" : 180.,
    "fsky" : 0.4
}

im_HIRAX_32_7 = {
    "type" : "intensity_mapping",
    "name" : "HIRAX_32_7",
    "nzfi" : "curves_IM/nz_HI.txt",
    "bzfi" : "curves_IM/bz_HI.txt",
    "szfi" : "curves_IM/sz_HI.txt",
    "ezfi" : "curves_IM/ez_HI.txt",
    "tzfi" : "curves_IM/tz_HI.txt",
    "dish_size" : 6.,
    "t_inst" : 50.,
    "t_total" : 10000.,
    "n_dish" : 1024,
    "area_eff" : 1.0,
    "im_type" : "interferometer",
    "base_file" : "curves_IM/baseline_file_HIRAX_7m.txt",
    "base_min" : 6.,
    "base_max" : 180.,
    "fsky" : 0.4
}

im_SKA_SD = {
    "type" : "intensity_mapping",
    "name" : "SKA_SD",
    "nzfi" : "curves_IM/nz_HI.txt",
    "bzfi" : "curves_IM/bz_HI.txt",
    "szfi" : "curves_IM/sz_HI.txt",
    "ezfi" : "curves_IM/ez_HI.txt",
    "tzfi" : "curves_IM/tz_HI.txt",
    "dish_size" : 15.,
    "t_inst" : 25.,
    "t_total" : 10000.,
    "n_dish" : 197,
    "area_eff" : 1.0,
    "im_type" : "single_dish",
    "base_file" : "none",
    "base_min" : 0.,
    "base_max" : 15.,
    "fsky" : 0.4
}

im_SKA_IF = {
    "type" : "intensity_mapping",
    "name" : "SKA_IF",
    "nzfi" : "curves_IM/nz_HI.txt",
    "bzfi" : "curves_IM/bz_HI.txt",
    "szfi" : "curves_IM/sz_HI.txt",
    "ezfi" : "curves_IM/ez_HI.txt",
    "tzfi" : "curves_IM/tz_HI.txt",
    "dish_size" : 15.,
    "t_inst" : 25.,
    "t_total" : 10000.,
    "n_dish" : 197,
    "area_eff" : 1.0,
    "im_type" : "interferometer",
    "base_file" : "curves_IM/baseline_file_SKA.txt",
    "base_min" : 15.,
    "base_max" : 1000.,
    "fsky" : 0.4
}

im_SKA = {
    "type" : "intensity_mapping",
    "name" : "SKA",
    "nzfi" : "curves_IM/nz_HI.txt",
    "bzfi" : "curves_IM/bz_HI.txt",
    "szfi" : "curves_IM/sz_HI.txt",
    "ezfi" : "curves_IM/ez_HI.txt",
    "tzfi" : "curves_IM/tz_HI.txt",
    "dish_size" : 15.,
    "t_inst" : 25.,
    "t_total" : 10000.,
    "n_dish" : 197,
    "area_eff" : 1.0,
    "im_type" : "hybrid",
    "base_file" : "curves_IM/baseline_file_SKA.txt",
    "base_min" : 0.,
    "base_max" : 1000.,
    "fsky" : 0.4
}

im_MeerKAT_SD = {
    "type" : "intensity_mapping",
    "name" : "MeerKAT_SD",
    "nzfi" : "curves_IM/nz_HI.txt",
    "bzfi" : "curves_IM/bz_HI.txt",
    "szfi" : "curves_IM/sz_HI.txt",
    "ezfi" : "curves_IM/ez_HI.txt",
    "tzfi" : "curves_IM/tz_HI.txt",
    "dish_size" : 13.5,
    "t_inst" : 25.,
    "t_total" : 4000.,
    "n_dish" : 64,
    "area_eff" : 1.0,
    "im_type" : "single_dish",
    "base_file" : "none",
    "base_min" : 0.,
    "base_max" : 15.,
    "fsky" : 0.1
}

im_MeerKAT_IF = {
    "type" : "intensity_mapping",
    "name" : "MeerKAT_IF",
    "nzfi" : "curves_IM/nz_HI.txt",
    "bzfi" : "curves_IM/bz_HI.txt",
    "szfi" : "curves_IM/sz_HI.txt",
    "ezfi" : "curves_IM/ez_HI.txt",
    "tzfi" : "curves_IM/tz_HI.txt",
    "dish_size" : 13.5,
    "t_inst" : 25.,
    "t_total" : 4000.,
    "n_dish" : 64,
    "area_eff" : 1.0,
    "im_type" : "interferometer",
    "base_file" : "curves_IM/baseline_file_MeerKAT.txt",
    "base_min" : 15.,
    "base_max" : 1000.,
    "fsky" : 0.1
}

im_MeerKAT = {
    "type" : "intensity_mapping",
    "name" : "MeerKAT",
    "nzfi" : "curves_IM/nz_HI.txt",
    "bzfi" : "curves_IM/bz_HI.txt",
    "szfi" : "curves_IM/sz_HI.txt",
    "ezfi" : "curves_IM/ez_HI.txt",
    "tzfi" : "curves_IM/tz_HI.txt",
    "dish_size" : 13.5,
    "t_inst" : 25.,
    "t_total" : 4000.,
    "n_dish" : 64,
    "area_eff" : 1.0,
    "im_type" : "hybrid",
    "base_file" : "curves_IM/baseline_file_MeerKAT.txt",
    "base_min" : 0.,
    "base_max" : 1000.,
    "fsky" : 0.1
}

cmb_S3_opt = {
    "type" : "cmb_lensing",
    "name" : "S3_opt",
    "sigma_t" : 8.5,
    "beam_amin" : 1.4,
    "lmin" : 2,
    "fsky" : 0.4,
}

cmb_S4_opt = {
    "type" : "cmb_lensing",
    "name" : "S4_opt",
    "sigma_t" : 1.0,
    "beam_amin" : 1.4,
    "lmin" : 2,
    "fsky" : 0.4,

}

cmb_S3 = {
    "type" : "cmb_lensing",
    "name" : "S3",
    "sigma_t" : 8.5,
    "beam_amin" : 1.4,
    "lmin" : 30,
    "fsky" : 0.4,
}

cmb_S4 = {
    "type" : "cmb_lensing",
    "name" : "S4",
    "sigma_t" : 1.0,
    "beam_amin" : 1.4,
    "lmin" : 30,
    "fsky" : 0.4,
}

cmb_none = {
    "type" : "cmb_lensing",
    "name" : "none",
    "sigma_t" : 1.0,
    "beam_amin" : 1.4,
    "lmin" : 30,
    "fsky" : 0.4,
}
