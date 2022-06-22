
;download full sun level 1 data 

t0='2013/04/23 20:30:00'
t1='2013/04/23 20:50:00'

;searchfile = vso_search(t0,t1, wave='94 Angstrom',inst='aia')
;searchfile = vso_search(t0,t1, wave='magnetogram',inst='HMI')
;getfile=vso_get(searchfile,out_dir="/media/rohit/MWA/20140914/EUV/fits/")

;for n=0,5 do begin
;	  print, n
;  endfor

;make a full sun IDL map structure

fitsfiles=findfile("/sdata/fits/*.193*.fits")
;fitsfiles=findfile("/sdata/20140914_hmi/hmi/hmi_20140914_*_magnetogram.fits")
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
;for i=10,n do begin 
for i=0,n do begin 
	print,fitsfiles(i)
	read_sdo,fitsfiles(i),index,data
	aia_prep,index,data,final_index,final_data,/normalize
	;final_data=reverse(final_data) ; For HMI)
	;final_data=reverse(final_data,2) ; For HMI
	print,max(final_data)
	index2map,final_index,final_data,fullsunmap
	;----------------------------------------------------------------
	read_sdo,reffile,refindex,refdata
	aia_prep,refindex,refdata,ref_final_index,ref_final_data,/normalize
	;ref_final_data=reverse(ref_final_data) ; For HMI)
	;ref_final_data=reverse(ref_final_data,2) ; For HMI
	index2map,ref_final_index,ref_final_data,ref_fullsunmap
	;sub_map,fullsunmap,submap,xrange=xrange,yrange=yrange
	;sub_map,ref_fullsunmap,ref_submap,xrange=xrange,yrange=yrange
	;---------------------------------------------------------------
	drot_map_=drot_map(ref_fullsunmap,ref_map=fullsunmap);,/KEEP_LIMB)
	dd1=drot_map_.data-final_data
	;index2map,ref_final_index,dd1,bdiff_map
	;---------------------------------------------------------------
	read_sdo,fitsfiles(i-10),indexd,datad
	aia_prep,indexd,datad,final_indexd,final_datad,/normalize
	index2map,final_indexd,final_datad,ref_fullsunmapd
	drot_map_=drot_map(ref_fullsunmapd,ref_map=fullsunmap);,/KEEP_LIMB)
	dd2=drot_map_.data-final_data
	;index2map,ref_final_index,dd2,ddiff_map
	;sub_map,drot_map_,submap,xrange=xrange,yrange=yrange
	;SAVE, bdiff_map.data, FILENAME = fitsfiles(i)+'_bdiff.sav'
	;SAVE, ddiff_map.data, FILENAME = fitsfiles(i)+'_ddiff.sav'
	SAVE, dd1, FILENAME = fitsfiles(i)+'_bdiff.sav'
	SAVE, dd2, FILENAME = fitsfiles(i)+'_ddiff.sav'
	;SAVE, submap, FILENAME = fitsfiles(i)+'submap.sav'
	endfor

;RESTORE, '/media/rohit/MWA/20140914/mwa_maps/mwa_240_2:48:34.0.fitsrot.sav'

end

