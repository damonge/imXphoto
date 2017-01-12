

phoz_LSSTgold = {
    "name" : "LSST_gold",
    "nzfi" : "curves_LSST/nz_shear_fiducial.txt",
    "bzfi" : "curves_LSST/bz_gold.txt",
    "szfi" : "curves_LSST/sz_gold.txt",
    "ezfi" : "curves_LSST/ez_gold.txt",
    "fsky" : 0.4
    }

spec_DESI = {
    "name" : "DESI",
    "nzfi" : "curves_spec/nz_DESI.txt",
    "bzfi" : "curves_spec/bz_DESI.txt",
    "szfi" : "curves_spec/sz_DESI.txt",
    "ezfi" : "curves_spec/ez_DESI.txt",
    "fsky" : 0.2
    }

spec_Euclid = {
    "name" : "Euclid",
    "nzfi" : "curves_spec/nz_Euclid.txt",
    "bzfi" : "curves_spec/bz_Euclid.txt",
    "szfi" : "curves_spec/sz_Euclid.txt",
    "ezfi" : "curves_spec/ez_Euclid.txt",
    "fsky" : 0.4
    }

spec_WFIRST = {
    "name" : "WFIRST",
    "nzfi" : "curves_spec/nz_WFIRST.txt",
    "bzfi" : "curves_spec/bz_WFIRST.txt",
    "szfi" : "curves_spec/sz_WFIRST.txt",
    "ezfi" : "curves_spec/ez_WFIRST.txt",
    "fsky" : 0.05
    }

im_HIRAX_32_6 = {
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
    "fsky" : 0.4
}

im_HIRAX_32_7 = {
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
    "fsky" : 0.4
}

im_SKA_SD = {
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
    "fsky" : 0.4
}

im_SKA_IF = {
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
    "fsky" : 0.4
}

im_SKA = {
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
    "fsky" : 0.4
}

im_MeerKAT_SD = {
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
    "fsky" : 0.1
}

im_MeerKAT_IF = {
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
    "fsky" : 0.1
}

im_MeerKAT = {
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
    "fsky" : 0.1
}
