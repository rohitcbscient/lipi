
;download full sun level 1 data 

t0='2014/09/14 01:30:00'
t1='2014/09/14 06:00:00'

;searchfile = vso_search(t0,t1, wave='94 Angstrom',inst='aia')
;searchfile = vso_search(t0,t1, wave='magnetogram',inst='HMI')
;getfile=vso_get(searchfile,out_dir="/media/rohit/MWA/20140914/EUV/fits/")

;for n=0,5 do begin
;	  print, n
;  endfor

;make a full sun IDL map structure

;fitsfiles=findfile("/media/rohit/MWA/20140914/EUV/fits/*171*.fits")
fitsfiles=findfile("/sdata/20140914_hmi/hmi/hmi_20140914_*_magnetogram.fits")
;fitsfiles=findfile("/media/rohit/MWA/20140914/mwa_maps/*_240_*.fits")
reffile=fitsfiles(0)
;ref_time=02_32_18.84;300
;xrange=[-300,100]
;yrange=[-800,-400]
xrange=[600,800]
yrange=[-400,-200]
;xrange=[705,855]
;yrange=[-77,-277]
n=n_elements(fitsfiles)
i=10
print,fitsfiles(i)
read_sdo,fitsfiles(i),index,data
aia_prep,index,data,final_index,final_data,/normalize
;final_data=reverse(final_data) ; For HMI)
;final_data=reverse(final_data,2) ; For HMI
print,max(final_data)
index2map,final_index,final_data,fullsunmap
read_sdo,reffile,refindex,refdata
aia_prep,refindex,refdata,ref_final_index,ref_final_data,/normalize
;ref_final_data=reverse(ref_final_data) ; For HMI)
;ref_final_data=reverse(ref_final_data,2) ; For HMI
index2map,ref_final_index,ref_final_data,ref_fullsunmap
;sub_map,fullsunmap,submap,xrange=xrange,yrange=yrange
;sub_map,ref_fullsunmap,ref_submap,xrange=xrange,yrange=yrange
drot_map_=drot_map(fullsunmap,ref_map=ref_fullsunmap,/KEEP_LIMB)
sub_map,drot_map_,submap,xrange=xrange,yrange=yrange
gx_fov2box, '14-Sep-14 01:54:00', center_arcsec=[700,-300], size_pix=[334,334,167], dx_km=726, /cea, /uv, /euv, tmp_dir= '.', out_dir= '.'
;RESTORE, '/media/rohit/MWA/20140914/mwa_maps/mwa_240_2:48:34.0.fitsrot.sav'

end

