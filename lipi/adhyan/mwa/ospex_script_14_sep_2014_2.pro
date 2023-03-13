; OSPEX script created Wed Nov 30 09:05:58 2022 by OSPEX writescript method.                 
;                                                                                            
;  Call this script with the keyword argument, obj=obj to return the                         
;  OSPEX object reference for use at the command line as well as in the GUI.                 
;  For example:                                                                              
;     ospex_script_14_sep_2014_2, obj=obj                                                    
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
;pro ospex_script_14_sep_2014_2, obj=obj    

; 5 mins from 02:00 to 02:05, we average over 30 sec, then from 02:05 to 02:15 we average over 10 sec
; Total images = 10 + 6*10 = 70 images
;For 30 sec integration
for j=0b,5 do begin
j=j*1
for i=0b,1 do begin
ii=i*30
ie=ii+29

; Below is the 10 sec section of the for loop
;for j=6b,14 do begin
;j=j*1
;for i=0b,5 do begin
;ii=i*10
;ie=ii+9
t1_name=strjoin('02:'+string(j,format='(I3.2)')+':'+string(ii,format='(I3.2)'))
t2_name=strjoin('02:'+string(j,format='(I3.2)')+':'+string(ie,format='(I3.2)'))
t1=strjoin('02:'+StrTrim(j,2)+':'+StrTrim(ii,2))
t2=strjoin('02:'+StrTrim(j,2)+':'+StrTrim(ie,2))
tname=strjoin(strsplit(t1_name,/extract,': '))+'_'+strjoin(strsplit(t2_name,/extract,': '))
print,t1_name,'///',t2_name,'///',tname
if not is_class(obj,'SPEX',/quiet) then obj = ospex()   
obj-> set, spex_fit_manual=0, spex_autoplot_enable=0, spex_fitcomp_plot_resid=0, spex_fit_progbar=1                                     
obj-> set, $                                                                                 
 spex_specfile= '/media/rohit/MWA/20140914/rhessi/hsi_spectrum_20140914_014500.fits'         
obj-> set, spex_drmfile= '/media/rohit/MWA/20140914/rhessi/hsi_srm_20140914_014500.fits'    
obj-> set, spex_source_angle= 90.0000                                                        
obj-> set, spex_source_xy= [940.435, -202.580]                                               
obj-> set, spex_fit_manual= 2                                                                
obj-> set, spex_fit_start_method= 'previous_iter'                                            
obj-> set, spex_erange= [8.0000000D, 40.000000D]                                             
obj-> set, spex_fit_time_interval= [strjoin('14-Sep-2014 '+t1), strjoin('14-Sep-2014 '+t2)]
obj-> set, spex_bk_time_interval=['14-Sep-2014 01:40:00.000', '14-Sep-2014 01:45:00.000']    
obj-> set, mcurvefit_itmax= 800L                                                             
obj-> set, spex_uncert= 0.000100000                                                          
obj-> set, fit_function= 'bpow+vth'                                                          
obj-> set, fit_comp_params= [0.115391, 4.78865, 19.0751, 6.04168, 0.0820958, 0.842910, $     
 9.97226]                                                                                    
obj-> set, fit_comp_minima= [1.00000e-10, 1.70000, 10.0000, 1.70000, 1.00000e-20, 0.500000, $
 0.0100000]                                                                                  
obj-> set, fit_comp_maxima= [1.00000e+10, 10.0000, 30.0000, 10.0000, 1.00000e+20, 8.00000, $ 
 10.0000]                                                                                    
obj-> set, fit_comp_free_mask= [1B, 1B, 1B, 1B, 1B, 1B, 1B]                                  
obj-> set, fit_comp_spectrum= ['', 'full']                                                   
obj-> set, fit_comp_model= ['', 'chianti']                                                   
obj-> set, spex_autoplot_bksub= 0                                                            
obj-> set, spex_autoplot_photons= 1                                                          
obj-> set, spex_autoplot_units= 'Flux'                                                       
obj-> set, spex_fitcomp_plot_units= 'Flux'                                                   
obj-> set, spex_fitcomp_plot_bk= 1                                                           
obj-> set, spex_fitcomp_plot_err= 1                                                          
obj-> set, spex_fitcomp_plot_photons= 1                                                      
obj-> set, spex_fitcomp_plot_resid= 0                                                        
obj-> set, spex_eband= [[3.00000, 6.00000], [6.00000, 12.0000], [12.0000, 25.0000], $        
 [25.0000, 50.0000], [50.0000, 100.000], [100.000, 300.000]]                                 
obj-> set, spex_tband= [['14-Sep-2014 01:45:00.000', '14-Sep-2014 01:53:45.000'], $          
 ['14-Sep-2014 01:53:45.000', '14-Sep-2014 02:02:30.000'], ['14-Sep-2014 02:02:30.000', $    
 '14-Sep-2014 02:11:15.000'], ['14-Sep-2014 02:11:15.000', '14-Sep-2014 02:20:00.000']]      
;obj -> restorefit, file='/media/rohit/MWA/20140914/rhessi/ospex_results_14_sep_2014.fits' 
obj-> set, spex_fit_manual=0, spex_autoplot_enable=0, spex_fitcomp_plot_resid=0, spex_fit_progbar=1
obj -> dofit, /all
obj->fitsummary, file_text='fitsummary_'+tname+'.txt'
;obj->filewrite, /fits, /buildsrm, srmfile=srmfilename, specfile = spfilename, all_simplify=0, /create
;obj-> plot_spectrum, /show_fit, /use_fitted, spex_units='flux', /bksub, /photon, /show_err
ct_energy = obj-> getaxis(/ct_energy)
rate = obj->calc_summ(item='data_photon_flux',errors='data_photon_flux_errors')
bkg_rate = obj->calc_summ(item='background_photon_flux',errors='background_photon_flux_errors')
spectrum_fits = obj -> calc_func_components(spex_unit='flux',/photons)
err_rate = obj->calc_summ(errors='data_photon_flux_errors')
err_bkg_rate = obj->calc_summ(errors='background_photon_flux_errors')
resid = obj->calc_summ(item='resid')
total_fit = spectrum_fits.yvals[*,0]
pow = spectrum_fits.yvals[*,1]
vth_fit = spectrum_fits.yvals[*,2]
obj-> set, vth_fit= vth_fit
obj-> set, pow= pow
obj-> savefit,outfile='fitspec_'+tname+'.fits'
writefits,'fitspec_'+tname+'_vth_pow.fits',[[ct_energy],[rate],[bkg_rate],[err_rate],[err_bkg_rate],[resid],[vth_fit],[pow],[total_fit]]
endfor
endfor
end                                                                                          
