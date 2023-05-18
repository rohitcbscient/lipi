o = ospex()
o->set, spex_specfile='/media/rohit/MWA/20140914/rhessi/hsi_spectrum_20140914_014500.fits'
o->set, spex_drmfile='/media/rohit/MWA/20140914/rhessi/hsi_srm_20140914_014500.fits'
o -> set, spex_bk_time_interval=['14-Sep-2014 01:40:00.000', '14-Sep-2014 01:45:00.000']
o->set, spex_eband=[[3.00000, 6.00000], [6.00000, 12.0000], [12.0000, 25.0000], $
 [25.0000, 50.0000], [50.0000, 100.000], [100.000, 300.000]]
o->set, spex_fit_time_interval=[['14-Sep-2014 02:14:00.000','14-Sep-2014 02:14:30.000']]
o->set, spex_erange=[13,30]
o->set, fit_function='vth+thick2'
o->set, fit_comp_param=[1.0e-005,1.,1.,  .5, 3., 45., 4.5]
o->set, fit_comp_free = [0,1,1, 1,1,1,1]
o->set, spex_fit_manual=0
o->set, spex_autoplot_enable=1
o-> set, mcurvefit_itmax= 800L
o-> set, spex_uncert= 0.000100000
o->dofit, /all
o->fitsummary, file_text='fit_test.txt'
end