;pro ospex2parameters, ospex_obj=ospex_obj
; ============================================================
; This programs take the ospex object and extracts the
; parameters of interest, for the
; 
; Event: 18-Feb-2021
; ============================================================



; Filename of the stored .sav file
;flpl='/home/battaglia/Documents/ETHZ/PhD/Flares/2021/02_February/18/A8/Data/SO_STIX/ospex/ospex-parameters_vth-thick2_Ec10keV_20210218_18:03:59-18:04:32_EarthUT-Time.sav'
;flpl='/home/battaglia/Documents/ETHZ/PhD/Flares/2021/02_February/18/A8/Data/SO_STIX/ospex/ospex-parameters_vth_20210218_18:04:32-18:04:46_EarthUT-Time.sav'
;flpl='/home/battaglia/Documents/ETHZ/PhD/Flares/2021/02_February/18/A8/Data/SO_STIX/ospex/ospex-parameters_vth-thick2_Ec10keV_L1data_20210218_18:03:59-18:04:32_EarthUT-Time.sav'
;flpl='/home/battaglia/Documents/ETHZ/PhD/Flares/2021/02_February/18/A8/Data/SO_STIX/ospex/ospex-parameters_vth_L1data_20210218_18:04:32-18:04:46_EarthUT-Time.sav'
;flpl='/home/afbattaglia/Documents/ETHZ/PhD/Flares/2021/02_February/18/A8/Data/SO_STIX/ospex/ospex-parameters_fixed-vth_high-from-DEM_20210218_18:03:19-18:03:39_EarthUT-Time.sav'
;flpl='/home/afbattaglia/Documents/ETHZ/PhD/Flares/2021/02_February/18/A8/Data/SO_STIX/ospex/ospex-parameters_fixed-vth_avg-from-DEM_20210218_18:03:19-18:03:39_EarthUT-Time.sav'
;flpl='/home/afbattaglia/Documents/ETHZ/PhD/Flares/2021/02_February/18/A8/Data/SO_STIX/ospex/ospex-parameters_fixed-vth_low-from-DEM_20210218_18:03:19-18:03:39_EarthUT-Time.sav'
flpl='/media/rohit/MWA/20140914/rhessi/counts_obs.sav'


; GET THE PARAMETERS
param_ospex = ospex_obj -> get(/spex_summ)
; Below, some useful parameters
; (see more: https://hesperia.gsfc.nasa.gov/ssw/packages/spex/doc/ospex_params_info.htm)
;
; param_ospex.spex_summ_chisq                Chi-square for each interval (ntime)
; param_ospex.spex_summ_fit_function         Fit function used
; param_ospex.spex_summ_params               Final fit parameters for each interval (nparams,ntime)
; param_ospex.spex_summ_resid                Residuals for each interval (nenergy,ntime)
; param_ospex.spex_summ_sigmas               Sigma in fit parameters for each interval (nparams,ntime)
; param_ospex.spex_summ_time_interval        Array of time intervals fit (2,ntime)

; Function used
funct_used = param_ospex.SPEX_SUMM_FIT_FUNCTION

; Energy axis
e_axis = average(param_ospex.spex_summ_energy,1)

; Width energy bin
e_bin_width = param_ospex.spex_summ_energy[1,*]-param_ospex.spex_summ_energy[0,*]

; Time interval
time_interval = param_ospex.SPEX_SUMM_TIME_INTERVAL

; Residuals
residuals = param_ospex.SPEX_SUMM_RESID

; Fit parameters
param_fit = param_ospex.spex_summ_params

; Sigma of the parameters
sigma_param = param_ospex.spex_summ_sigmas

; Chi square
chi2 = param_ospex.spex_summ_chisq

; Detector area (cm^2)
det_area = param_ospex.SPEX_SUMM_AREA

; Conversion  factors for each interval (nenergy,ntime), given in counts/photon
; To get photon spectrum, divide the count spectrum by this array
conv_counts_photons = param_ospex.spex_summ_conv

; Observed STIX spectrum (s-1 cm-2 keV-1)
; OSPEX gives back only the count rate (counts/sec)
stix_spectrum = param_ospex.SPEX_SUMM_CT_RATE/(det_area*e_bin_width)

; Error on the STIX spectrum
stix_spectrum_error = param_ospex.SPEX_SUMM_CT_ERROR/(det_area*e_bin_width)

; Photon spectrum (already in phot/cm^2/sec/keV)
photon_spectrum = param_ospex.spex_summ_ph_model



;
;; **************
;; Nonthermal energy
;func = 'thick2'
;integration_time = time_interval[1]-time_interval[0] ; seconds
;nth_param = param_fit[3:8]
;nth_power = calc_nontherm_electron_energy_flux(nth_param, func=func)
;print,''
;print,'Lower limit of the nonthermal energy [erg]'
;print,nth_power * integration_time
;print,''
;; **************


; GET THE SPECTRA OF THE FITS
spectrum_fits = ospex_obj -> calc_func_components(spex_unit='flux')
;sp2 = ospex_obj -> calc_func_components(spex_unit='flux')
;IDL> help, spectra, /str
;** Structure <656e058>, 4 tags, length=640, data length=640, refs=1:
;YVALS           DOUBLE    Array[17, 3]  ; [energy channels, fits spectra (0: total, 1: vth, 2: ...)]
;ID              STRING    Array[3]
;STRAT_ARR       STRING    Array[3]
;CT_ENERGY       FLOAT     Array[2, 17]

; Fit total (s-1 cm-2 keV-1)
;total_fit = spectrum_fits.yvals[*,0]

; Fit single thermal (s-1 cm-2 keV-1)
;vth_fit = spectrum_fits.yvals[*,1]
vth_fit = spectrum_fits.yvals[*,0]

; Fit thick target (s-1 cm-2 keV-1)
;thick2_fit = spectrum_fits.yvals[*,2]


; ****************************************************************************************
;save,filename=flpl,funct_used,e_axis,time_interval,stix_spectrum,stix_spectrum_error,photon_spectrum,$
;  residuals,param_fit,sigma_param,chi2,det_area,e_bin_width,conv_counts_photons,total_fit,vth_fit,thick2_fit
save,filename=flpl,funct_used,e_axis,time_interval,stix_spectrum,stix_spectrum_error,photon_spectrum,residuals,$
  param_fit,sigma_param,chi2,det_area,e_bin_width,conv_counts_photons,vth_fit
; ****************************************************************************************



end ; End of the script
