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
;For 1 min integration
for j=0b,5 do begin
j=j*1
for i=0b,0 do begin
ii=i
ie=ii+59

; Below is the 10 sec section of the for loop 14 & 5
;for j=5,15 do begin
;j=j*1
;for i=0b,2 do begin
;ii=i*20
;ie=ii+19
t1_name=strjoin('02:'+string(j,format='(I3.2)')+':'+string(ii,format='(I3.2)'))
t2_name=strjoin('02:'+string(j,format='(I3.2)')+':'+string(ie,format='(I3.2)'))
t1=strjoin('02:'+StrTrim(j,2)+':'+StrTrim(ii,2))
t2=strjoin('02:'+StrTrim(j,2)+':'+StrTrim(ie,2))
tname=strjoin(strsplit(t1_name,/extract,': '))+'_'+strjoin(strsplit(t2_name,/extract,': '))
print,t1_name,'///',t2_name,'///',tname
if not is_class(obj,'SPEX',/quiet) then obj = ospex()   
obj-> set, spex_fit_manual=0, spex_autoplot_enable=1, spex_fitcomp_plot_resid=0, spex_fit_progbar=1
obj-> set, spex_specfile= '/media/rohit/MWA/20140914/rhessi/hsi_spectrum_20140914_014500.fits'
obj-> set, spex_drmfile= '/media/rohit/MWA/20140914/rhessi/hsi_srm_20140914_014500.fits'
obj-> set, spex_erange=[9,50]
obj-> set, mcurvefit_itmax= 800L                                                             
obj-> set, spex_uncert= 0.000100000                                                          
obj-> set, fit_function= 'vth+thick2'
obj->set, fit_comp_free = [1,1,0, 1,1,0,1,1,1]
obj->set, fit_comp_params=[.005, 1.0, 1., 1.0, 1.0, 1000., 6.7, 15.0, 10000]
obj-> set, spex_eband= [[3.00000, 6.00000], [6.00000, 12.0000], [12.0000, 25.0000], $
 [25.0000, 50.0000], [50.0000, 100.000], [100.000, 300.000]]
obj-> set, spex_fit_time_interval= [strjoin('14-Sep-2014 '+t1), strjoin('14-Sep-2014 '+t2)]
obj-> set, spex_bk_time_interval=['14-Sep-2014 01:51:00.000', '14-Sep-2014 01:58:00.000']
obj -> dofit, /all
obj->fitsummary, file_text='fitsummary_'+tname+'.txt'
ct_energy = obj-> getaxis(/ct_energy)
rate = obj->calc_summ(item='data_photon_flux',errors=err_ph)
err_rate = err_ph
print,rate,err_ph
bkg_rate = obj->calc_summ(item='background_photon_flux',errors=back_err_ph)
print,bkg_rate,back_err_ph
spectrum_fits = obj -> calc_func_components(spex_unit='flux',/photons)
err_bkg_rate = back_err_ph
electron_flux = obj -> calc_summ(errors='thick_integ_flux')
resid = obj->calc_summ(item='resid')
total_fit = spectrum_fits.yvals[*,0]
pow = spectrum_fits.yvals[*,1]
vth_fit = spectrum_fits.yvals[*,2]
obj-> set, vth_fit= vth_fit
obj-> set, pow= pow
sum_params = obj -> get(/spex_summ_params)
obj-> savefit,outfile='fitspec_'+tname+'.fits'
s = obj->get(/spex_summ)
param_fit = s.spex_summ_params
print,rate,err_ph
writefits,'fitspec_'+tname+'_vth_pow.fits',[[ct_energy],[rate],[bkg_rate],[err_rate],[err_bkg_rate],[resid],[vth_fit],[pow],[total_fit]]
endfor
endfor
end                                                                                          
