; OSPEX script created Tue Nov 29 17:56:14 2022 by OSPEX writescript method.              
;                                                                                         
;  Call this script with the keyword argument, obj=obj to return the                      
;  OSPEX object reference for use at the command line as well as in the GUI.              
;  For example:                                                                           
;     ospex_script_14_sep_2014_1, obj=obj                                                 
;                                                                                         
;  Note that this script simply sets parameters in the OSPEX object as they               
;  were when you wrote the script, and optionally restores fit results.                   
;  To make OSPEX do anything in this script, you need to add some action commands.        
;  For instance, the command                                                              
;     obj -> dofit, /all                                                                  
;  would tell OSPEX to do fits in all your fit time intervals.                            
;  See the OSPEX methods section in the OSPEX documentation at                            
;  http://hesperia.gsfc.nasa.gov/ssw/packages/spex/doc/ospex_explanation.htm              
;  for a complete list of methods and their arguments.                                    
;                                                                                         
pro ospex_script_14_sep_2014_1, obj=obj                                                   
if not is_class(obj,'SPEX',/quiet) then obj = ospex()                                     
obj-> set, $                                                                              
 spex_specfile= '/media/rohit/MWA/20140914/rhessi/hsi_spectrum_20140914_014000.fits'      
obj-> set, spex_drmfile= '/media/rohit/MWA/20140914/rhessi/hsi_srm_20140914_014000.fits'  
obj-> set, spex_source_angle= 90.0000                                                     
obj-> set, spex_source_xy= [940.435, -202.580]                                            
obj-> set, spex_fit_start_method= 'previous_start'                                        
obj-> set, spex_erange= [6.0000000D, 46.000000D]                                          
obj-> set, spex_fit_time_interval= ['14-Sep-2014 02:09:00.000', $                         
 '14-Sep-2014 02:11:00.000']                                                              
obj-> set, spex_bk_time_interval=['14-Sep-2014 01:40:00.000', '14-Sep-2014 01:48:00.000'] 
obj-> set, mcurvefit_itmax= 200L                                                          
obj-> set, spex_uncert= 0.00000                                                           
obj-> set, fit_function= 'vth+1pow'                                                       
obj-> set, fit_comp_params= [40.0121, 0.642378, 1.00000, 0.123726, 6.33191, 38.2334]      
obj-> set, fit_comp_minima= [1.00000e-20, 0.500000, 0.0100000, 1.00000e-10, 1.70000, $    
 5.00000]                                                                                 
obj-> set, fit_comp_maxima= [1.00000e+20, 8.00000, 10.0000, 1.00000e+10, 10.0000, 20.0000]
obj-> set, fit_comp_free_mask= [1B, 1B, 0B, 1B, 1B, 1B]                                   
obj-> set, fit_comp_spectrum= ['full', '']                                                
obj-> set, fit_comp_model= ['chianti', '']                                                
obj-> set, spex_autoplot_photons= 1                                                       
obj-> set, spex_autoplot_units= 'Flux'                                                    
obj-> set, spex_fitcomp_plot_bk= 1                                                        
obj-> set, spex_fitcomp_plot_err= 1                                                       
obj-> set, spex_fitcomp_plot_photons= 1                                                   
obj-> set, spex_eband= [[3.00000, 6.00000], [6.00000, 12.0000], [12.0000, 25.0000], $     
 [25.0000, 50.0000], [50.0000, 100.000], [100.000, 300.000]]                              
obj-> set, spex_tband= [['14-Sep-2014 01:40:00.000', '14-Sep-2014 01:50:00.000'], $       
 ['14-Sep-2014 01:50:00.000', '14-Sep-2014 02:00:00.000'], ['14-Sep-2014 02:00:00.000', $ 
 '14-Sep-2014 02:10:00.000'], ['14-Sep-2014 02:10:00.000', '14-Sep-2014 02:20:00.000']]   
obj -> restorefit, file='/media/rohit/MWA/20140914/rhessi/ospex_results_14_sep_2014.fits' 
end                                                                                       
