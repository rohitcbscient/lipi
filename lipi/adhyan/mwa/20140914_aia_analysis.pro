
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

fitsfiles=findfile("/media/rohit/MWA/20140914/EUV/fits/*171*.fits")
;fitsfiles=findfile("/media/rohit/MWA/20140914/mwa_maps/*_240_*.fits")
reffile=fitsfiles(300)
;ref_time=02_32_18.84;300
;xrange=[-300,100]
;yrange=[-800,-400]
xrange=[550,1150]
yrange=[-700,-100]
n=n_elements(fitsfiles)
for i=0,5 do begin 
	print,fitsfiles(i)
	read_sdo,fitsfiles(i),index,data
	aia_prep,index,data,final_index,final_data,/normalize
	print,max(final_data)
	index2map,final_index,final_data,fullsunmap
	read_sdo,reffile,refindex,refdata
	aia_prep,refindex,refdata,ref_final_index,ref_final_data,/normalize
	index2map,ref_final_index,ref_final_data,ref_fullsunmap
	;sub_map,fullsunmap,submap,xrange=xrange,yrange=yrange
	;sub_map,ref_fullsunmap,ref_submap,xrange=xrange,yrange=yrange
	drot_map_=drot_map(fullsunmap,ref_map=ref_fullsunmap,/KEEP_LIMB)
	sub_map,drot_map_,submap,xrange=xrange,yrange=yrange
	SAVE, drot_map_, FILENAME = fitsfiles(i)+'rot.sav'
	;SAVE, sub_map, FILENAME = fitsfiles(i)+'submap.sav'
	;SAVE, drot_map_, FILENAME = fitsfiles(i)+'rot.fits'
	endfor

;RESTORE, '/media/rohit/MWA/20140914/mwa_maps/mwa_240_2:48:34.0.fitsrot.sav'

end

