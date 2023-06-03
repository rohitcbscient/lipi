;+                                                                                           
; hsi_image script - created Mon Nov 28 11:35:25 2022 by hsi_params2script.pro               
; Includes all control parameters.                                                           
;                                                                                            
; Note: This script simply sets control parameters in the hsi_image object as they           
;  were when you wrote the script.  To make the object do anything in this script,           
;  you need to add some action commands.  For instance, the command                          
;    image= obj->getdata()                                                                   
;  would generate and return the image.                                                      
;                                                                                            
; For a complete list of control parameters look at the tables in                            
; http://hesperia.gsfc.nasa.gov/ssw/hessi/doc/hsi_params_all.htm                             
;                                                                                            
; There are several ways to use this script:                                                 
;                                                                                            
; A. Run this procedure via this command (compile first if it's not in your IDL path):       
;     hsi_image_script, obj=obj                                                              
;       Note: you can set additional parameters or override parameters in script via:        
;       hsi_image_script, obj=obj, param1=param1, param2=param2,...                          
; or                                                                                         
; B. Compile and execute it as a main program by:                                            
;    1. Comment out the "pro" line and save.                                                 
;    2. From the command line, type .run hsi_image_script.                                   
;       In the IDLDE, click Compile, then Execute.                                           
;    3. Use .GO to restart at the beginning of the script                                    
; or                                                                                         
; C. Execute it as a batch file via:                                                         
;    1. Comment out the "pro" and "end" lines and save.                                      
;    2. Run the commands via @hsi_image_script                                               
;                                                                                            
; Once you have run the script (via one of those 3 methods), you will have an                
; hsi_image object called obj that is set up as it was when you wrote the script.            
; You can proceed using obj from the command line or the hessi GUI.                          
; To use it in the GUI, type                                                                 
;   hessi,obj                                                                                
; To use it from the command line, here are a few examples of commands:                      
;  data = obj->getdata()   ; retrieve the last image made                                    
;  data = obj->getdata(use_single=0)  ; retrieve all images in cube                          
;  obj->plot               ; plot the last image                                             
;  obj->plotman            ; plot the last image in plotman                                  
;  obj->plotman, /choose   ; choose which image(s) in cube to plot in plotman                
;-                                                                                           

;pro hsi_image_script, obj=obj, _extra=_extra                                                 
 
search_network,/enable
for j=0b,15 do begin
j=j*1
for i=0b,0 do begin
ii=i*0+20
ie=ii+0
t1=strjoin('02:'+StrTrim(j,2)+':'+StrTrim(ii,2))
t2=strjoin('02:'+StrTrim(j+1,2)+':'+StrTrim(ie,2))
tname=strjoin(strsplit(t1,/extract,':'))+'_'+strjoin(strsplit(t2,/extract,':')) 
print,t1,t2,tname  
                                                                                                                                                                                                                                       
obj = hsi_image()

obj-> set, vis_input_fits= ''                                                                
obj-> set, cbe_normalize= 0                                                                  
obj-> set, im_time_bin= 0.0000000D                                                           
obj-> set, im_time_ref= 0.0000000D                                                           
obj-> set, full_info= 0B                                                                     
obj-> set, flare_id_nr= 0L                                                                   
obj-> set, image_algorithm= 'Clean'                                                          
obj-> set, nvis_min= 10                                                                      
obj-> set, noquintic_interp= 0                                                               
obj-> set, imaging_method= ''                                                                
obj-> set, profile_show_plot= 0                                                              
obj-> set, profile_ps_plot= 0                                                                
obj-> set, profile_jpeg_plot= 0                                                              
obj-> set, profile_plot_dir= ''                                                              
obj-> set, profile_plot_rate= 1                                                              
obj-> set, profile_plot_resid= 0                                                             
obj-> set, mc_ntrials= 0                                                                     
obj-> set, mc_show_plot= 0                                                                   
obj-> set, viscomp_show_plot= 0                                                              
obj-> set, viscomp_ps_plot= 0                                                                
obj-> set, viscomp_jpeg_plot= 0                                                              
obj-> set, viscomp_plot_dir= ''                                                              
obj-> set, skip_viscomp_chisq= 0                                                             
obj-> set, snr_chk= 2                                                                        
obj-> set, snr_thresh= 2.00000                                                               
obj-> set, snr_vis_file= ''                                                                  
obj-> set, old_cbe_normalize= 0                                                              
obj-> set, mpat_coord= 'ANNSEC'                                                              
obj-> set, factor_by= 1                                                                      
obj-> set, modpat_skip= 1                                                                    
obj-> set, r0_offset= 2560.00                                                                
obj-> set, image_dim= [128, 128]
obj-> set, pixel_scale= 1.00000                                                              
obj-> set, pixel_size= [2.00000, 2.00000]
obj-> set, im_time_interval= [strjoin('14-Sep-2014 '+t1), $                         
 strjoin('14-Sep-2014 '+t2)]        
obj-> set, cbe_filename= ''                                                                  
obj-> set, use_auto_time_bin= 1L                                                             
obj-> set, cbe_digital_quality= 0.950000                                                     
obj-> set, cbe_powers_of_two= 1L                                                             
obj-> set, cbe_time_bin_floor= [0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L]                          
obj-> set, user_flux_var= 0                                                                  
obj-> set, use_flux_var= 1L                                                                  
obj-> set, smoothing_time= 4.00000                                                           
obj-> set, srt_filename= 'SRT_2*.dat'                                                        
obj-> set, use_time_window= [0.00000, 0.00000, 0.00000]                                      
obj-> set, use_phz_stacker= 0L                                                               
obj-> set, cb_coef= 1.00000                                                                  
obj-> set, user_hook= 0                                                                      
obj-> set, imaging_power_law= 4.00000                                                        
obj-> set, use_local_average= 0B                                                             
obj-> set, local_average_frequency= [[16.0000, 16.0000, 16.0000, 16.0000, 16.0000, 16.0000, $
 16.0000, 16.0000, 4.00000], [16.0000, 16.0000, 16.0000, 16.0000, 16.0000, 16.0000, $        
 16.0000, 16.0000, 4.00000], [16.0000, 16.0000, 16.0000, 16.0000, 16.0000, 16.0000, $        
 16.0000, 16.0000, 4.00000]]                                                                 
obj-> set, auto_frequency= 1B                                                                
obj-> set, max_harmonic= 1                                                                   
obj-> set, r_threshold= [0.600000, 0.600000, 0.600000, 0.600000, 0.600000, 0.600000, $       
 0.600000, 0.600000, 0.600000]                                                               
obj-> set, cbe_demodulate_engine= 'old'                                                      
obj-> set, cbe_max_corr= 0.00000                                                             
obj-> set, cbe_multi_atten_threshold= -1.00000                                               
obj-> set, aspect_p_error_threshold= 0.800000                                                
obj-> set, cbe_fast_multi_energy= 0                                                          
obj-> set, phz_n_roll_bins_min= [6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L]                         
obj-> set, phz_n_roll_bins_max= [64L, 64L, 64L, 64L, 64L, 64L, 64L, 64L, 64L]                
obj-> set, phz_n_roll_bins_control= [0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L]                     
obj-> set, phz_report_roll_bins= 0                                                           
obj-> set, phz_n_phase_bins= 12L                                                             
obj-> set, phz_radius= 60.0000                                                               
obj-> set, reference_position_angle= 0.00000                                                 
obj-> set, use_reference_position_angle= 0                                                   
obj-> set, use_stagger_rolls= 0                                                              
obj-> set, use_phase_conjugation= 0                                                          
obj-> set, flare_xyoffset= [0.00000, 0.00000]                                                
obj-> set, use_flare_xyoffset= 1                                                             
obj-> set, front_segment= 1B                                                                 
obj-> set, rear_segment= 0B                                                                  
obj-> set, im_energy_binning= [6.0000, 12.0000]
obj-> set, time_bin_def= [1.00000, 2.00000, 4.00000, 8.00000, 8.00000, 16.0000, 32.0000, $   
 64.0000, 128.000]                                                                           
obj-> set, time_bin_min= 512L                                                                
obj-> set, det_index_mask= [0B, 0B, 1B, 0B, 0B, 1B, 0B, 1B, 1B]
obj-> set, poisson= 0B                                                                       
obj-> set, seed= 0.00000                                                                     
obj-> set, coincidence_flag= 0                                                               
obj-> set, rebin_method= ''                                                                  
obj-> set, rebin_poisson_flag= 0                                                             
obj-> set, rebin_seed= 0                                                                     
obj-> set, contig_energy_edges= 1                                                            
obj-> set, livetime_enable= 1                                                                
obj-> set, align512= 1                                                                       
obj-> set, sum_flag= 0                                                                       
obj-> set, sum_coincidence= 0                                                                
obj-> set, rear_no_anti= 1                                                                   
obj-> set, sp_dp_cutoff= 0.00000                                                             
obj-> set, decimation_correct= 1                                                             
obj-> set, rear_decimation_correct= 0                                                        
obj-> set, decim_apar= [1.00000, 2.00000, 1.00000, 4.00000, 400.000, 4.00000]                
obj-> set, use_cull= 1                                                                       
obj-> set, cull_frac= 0.500000                                                               
obj-> set, use_total_count= 0                                                                
obj-> set, clear_halfscale= 1                                                                
obj-> set, clear_front_fullscale= 0                                                          
obj-> set, gain_time_wanted= 0.0000000D                                                      
obj-> set, gain_generation= 1000                                                             
obj-> set, number_of_half_rotations= 0                                                       
obj-> set, time_unit= 1                                                                      
obj-> set, no_livetime= 0B                                                                   
obj-> set, ct_interpolate= 0B                                                                
obj-> set, deflt_atten_state= 1B                                                             
obj-> set, fr_deadtime_window= 0                                                             
obj-> set, min_time_4_off= 0.100000                                                          
obj-> set, evl_rollover_test= 12                                                             
obj-> set, dp_cutoff_max= 0.0500000                                                          
obj-> set, dp_cutoff_min= 0.000800000                                                        
obj-> set, dp_cutoff_coeff= 4.50000                                                          
obj-> set, dp_cutoff_xp= -0.900000                                                           
obj-> set, dp_extend_utlim= [8.3090880d+08, 8.5717440d+08, 8.4680640d+08, 8.5207680d+08, $   
 1.4516928d+09, 8.0464320d+08, 1.4516928d+09, 1.4516928d+09, 1.4516928d+09, 1.4516928d+09, $ 
 1.4516928d+09, 1.4516928d+09, 1.4516928d+09, 1.4516928d+09, 1.4516928d+09, 1.4516928d+09, $ 
 1.4516928d+09, 1.4516928d+09]                                                               
obj-> set, dp_extend_sec= [0.0100000, 0.0100000, 0.0100000, 0.0100000, 0.0200000, $          
 0.0100000, 0.0200000, 0.0200000, 0.0200000, 0.0200000, 0.0200000, 0.0200000, 0.0200000, $   
 0.0200000, 0.0200000, 0.0200000, 0.0200000, 0.0200000]                                      
obj-> set, dp_extend_def= [0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, $  
 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, $ 
 0.00000]                                                                                    
obj-> set, dp_append_def= [0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, $  
 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, $ 
 0.00000]                                                                                    
obj-> set, dp_prepend_def= [0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, $ 
 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, $ 
 0.00000]                                                                                    
obj-> set, dp_prepend_nvalid= [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]         
obj-> set, dp_append_nvalid= [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]          
obj-> set, extend_time_range= 2.00000                                                        
obj-> set, dp_enable= 1                                                                      
obj-> set, dp_lld= [51, 49, 52, 54, 52, 44, 65, 46, 40, 0, 0, 0, 0, 0, 0, 0, 0, 0]           
obj-> set, dp_uld= [6000L, 6000L, 6000L, 6000L, 6000L, 6000L, 6000L, 6000L, 6000L, 0L, 0L, $ 
 0L, 0L, 0L, 0L, 0L, 0L, 0L]                                                                 
obj-> set, ramp_epeak= 35.0000                                                               
obj-> set, ramp_ntest= 6                                                                     
obj-> set, ramp_emax= 300.000                                                                
obj-> set, ramp_emin= 20.0000                                                                
obj-> set, ramp_ebulk= 40.0000                                                               
obj-> set, atten0_test= 0.300000                                                             
obj-> set, no_csa_dropout= 1                                                                 
obj-> set, file_type= 'fits'                                                                 
obj-> set, check_bad_packet= 0B                                                              
obj-> set, adp_test= [0B, 0B, 0B, 0B, 0B, 0B, 0B, 0B, 0B, 0B]                                
obj-> set, packet_per_bunch= 5000L                                                           
obj-> set, aspect_mode= 0B                                                                   
obj-> set, aspect_cntl_level= 0                                                              
obj-> set, aspect_sim= 0B                                                                    
obj-> set, as_interpol= 'quad'                                                               
obj-> set, as_no_extrapol= 1                                                                 
obj-> set, ras_time_extension= [-1400.0000D, 1400.0000D]                                     
obj-> set, saszero= 0B                                                                       
obj-> set, as_spin_period= 4.0000000D                                                        
obj-> set, as_roll_offset= 0.0000000D                                                        
obj-> set, as_roll_solution= 'DBASE'                                                         
obj-> set, as_point_solution= 'SAS'                                                          
obj-> set, sc_sun_offset= [0.00000, 0.00000]                                                 
obj-> set, equat_ns= 0                                                                       
obj-> set, no_ptstop= 0                                                                      
obj-> set, pmtras_diagnostic= 0U                                                             
obj-> set, pmtras_method= 'ORIGINAL'                                                         
obj-> set, pmtras_time_option= 'CLIP2IMAGE'                                                  
obj-> set, pmtras_phase_lock= 0                                                              
obj-> set, pmtras_nstarlim= 0                                                                
obj-> set, pmtras_offset_phase= -10.0000                                                     
obj-> set, pmtras_phzlk_tol= 0.00200000                                                      
obj-> set, pmtras_stars_fraction= 0.100000                                                   
obj-> set, pmtras_add_planets= 0                                                             
obj-> set, pmtras_get_tdistance= 0                                                           
obj-> set, pmtras_period_estimator= 0                                                        
obj-> set, pmtras_deadtime= 0                                                                
obj-> set, flatfield= 1                                                                      
obj-> set, use_rate= 1                                                                       
obj-> set, weight= 1                                                                         
obj-> set, uniform_weighting= 0                                                              
obj-> set, spatial_frequency_weight= [1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, $
 1.00000, 1.00000, 1.00000]                                                                  
obj-> set, taper= 0.00000                                                                    
obj-> set, xy_pixel= 2185L                                                                   
obj-> set, psf_no_sum= 0B                                                                    
obj-> set, im_calc_error= 0B                                                                 
obj-> set, clean_niter= 100                                                                  
obj-> set, clean_more_iter= 0                                                                
obj-> set, clean_negative_max_test= 1                                                        
obj-> set, clean_chi_sq_min_test= 0                                                          
obj-> set, clean_frac= 0.0500000                                                             
obj-> set, clean_chi_sq_crit= -1.00000                                                       
obj-> set, clean_no_chi2= 1                                                                  
obj-> set, clean_show_maps= 1                                                                
obj-> set, clean_show_n_maps= 1                                                              
obj-> set, clean_show_map_xdim= 1024                                                         
obj-> set, clean_show_chi= 1                                                                 
obj-> set, clean_show_n_chi= 1                                                               
obj-> set, clean_box= 0                                                                      
obj-> set, clean_cw_list= 0                                                                  
obj-> set, clean_cw_nop= 0                                                                   
obj-> set, clean_cw_inverse= 0                                                               
obj-> set, clean_progress_bar= 1                                                             
obj-> set, clean_mark_box= 0                                                                 
obj-> set, clean_media_mode= 0                                                               
obj-> set, clean_beam_width_factor= 1.00000                                                  
obj-> set, clean_regress_combine= 'disable'                                                  
obj-> set, clean_input_dirty_map= 0.00000                                                    
obj-> set, clean_profile_no_resid= 1                                                         
obj-> set, nj_ferr= 0.00000                                                                  
obj-> set, nj_tol= 0.0300000                                                                 
obj-> set, nj_show_maps= 0                                                                   
obj-> set, nj_progress_bar= 1                                                                
obj-> set, vis_chi2lim= 1.00000e+09                                                          
obj-> set, vis_edit= 1B                                                                      
obj-> set, vis_conjugate= 1B                                                                 
obj-> set, vis_normalize= 1B                                                                 
obj-> set, vis_max_corr= 0.250000                                                            
obj-> set, vis_photon2el= 2                                                                  
obj-> set, vis_type= 'photon'                                                                
obj-> set, vis_out_filename= ''                                                              
obj-> set, vis_plotfit= 0B                                                                   
obj-> set, ge_silent= 1                                                                      
obj-> set, ge_auto_percent_lambda= 1                                                         
obj-> set, ge_percent_lambda= 0.0200000                                                      
obj-> set, ge_tol= 0.00100000                                                                
obj-> set, ge_max_iter= 1000L                                                                
obj-> set, vf_nophase= 0B                                                                    
obj-> set, vf_circle= 0B                                                                     
obj-> set, vf_noerr= 0B                                                                      
obj-> set, vf_absolute= 0B                                                                   
obj-> set, vf_maxiter= 2000L                                                                 
obj-> set, vf_multi= 0                                                                       
obj-> set, vf_loop= 0                                                                        
obj-> set, vf_noplotfit= 0                                                                   
obj-> set, uv_show_vismap= 1B                                                                
obj-> set, vis_cs_abort_numerical_issues= 0                                                  
obj-> set, vis_cs_sparseness= 0.500000                                                       
obj-> set, vis_wv_nscales= 3                                                                 
obj-> set, vis_wv_lam= 0.0500000                                                             
obj-> set, vis_wv_niter= 200                                                                 
obj-> set, vis_wv_silent= 0                                                                  
obj-> set, vis_wv_autolam= 1                 
obj->set, im_out_fits_filename= 'fitsimage_6-12keV_'+tname+'.fits'
if keyword_set(_extra) then obj->set, _extra=_extra   
obj->fitswrite 
endfor
endfor                                    
                                                                                             
end                                                                                          
