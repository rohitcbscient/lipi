
;download full sun level 1 data 

;t0='2016/04/09 18:40:00'
;t1='2016/04/09 18:50:00'

;searchfile = vso_search(t0,t1, wave='1700 Angstrom',inst='aia')
;getfile=vso_get(searchfile,out_dir="/media/rohit/VLA/20160409_EUV/full_sun/")

;for n=0,5 do begin
;	  print, n
;  endfor

;make a full sun IDL map structure

fitsfiles=findfile("/media/rohit/MWA/20201116_STIX_EUV/*94*.fits")
reffile=fitsfiles(140)
;xrange=[-300,100]
;yrange=[-800,-400]
xrange=[700,1100]
yrange=[-600,-200]
n=n_elements(fitsfiles)
for i=0,n do begin 
	print,fitsfiles(i)
	read_sdo,fitsfiles(i),index,data
	aia_prep,index,data,final_index,final_data,/normalize
	index2map,final_index,final_data,fullsunmap
	read_sdo,reffile,refindex,refdata
	aia_prep,refindex,refdata,ref_final_index,ref_final_data,/normalize
	index2map,ref_final_index,ref_final_data,ref_fullsunmap
	sub_map,fullsunmap,submap,xrange=xrange,yrange=yrange
	sub_map,ref_fullsunmap,ref_submap,xrange=xrange,yrange=yrange
	drot_map_=drot_map(submap,ref_map=ref_submap)
	SAVE, drot_map_, FILENAME = fitsfiles(i)+'rot.sav'
	SAVE, submap, FILENAME = fitsfiles(i)+'sub.sav'
	endfor

end

