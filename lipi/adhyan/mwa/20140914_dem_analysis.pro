
;download full sun level 1 data 

t0='2014/09/14 03:00:00'
t1='2014/09/14 05:30:00'

;searchfile = vso_search(t0,t1, wave='131 Angstrom',inst='aia')
;searchfile = vso_search(t0,t1, wave='magnetogram',inst='HMI')
;getfile=vso_get(searchfile,out_dir="/media/rohit/MWA/20140914/EUV/fits/")

;for n=0,5 do begin
;	  print, n
;  endfor

;make a full sun IDL map structure

fitsfiles=findfile("*.fits")
reffile=fitsfiles(0)
n=n_elements(fitsfiles)
aia_arr=FINDGEN(4096,4096,6)
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
	SAVE, fullsunmap, FILENAME = fitsfiles(i)+'_lev1.5.sav'
	aia_arr(*,*,i) = final_data
	endfor

save,aia_arr,filename='aia_arr_lev1.5.sav'


;RESTORE, '/media/rohit/MWA/20140914/mwa_maps/mwa_240_2:48:34.0.fitsrot.sav'

end

