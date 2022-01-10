

fitsfiles=findfile("F:\20140914_hmi\hmi\hmi_20140914_*_magnetogram.fitsrot.sav")
reffile=fitsfiles(0)
xrange=[570,770]
yrange=[-220,-420]
nz=335
n=n_elements(fitsfiles)
for i=0,n do begin
  print,fitsfiles(i)
  restore,fitsfiles(i)
  sub_map,drot_map_,submap,xrange=xrange,yrange=yrange
  gx_bz2lff,submap.data,nz,dr,Bout,alpha1=alpha1,seehafer=seehafer,sub_b00=sub_b00,sub_plan=sub_plan
  SAVE, Bout, FILENAME = fitsfiles(i)+'.cube.sav'
  endfor

;image=findgen(256,256)
;map=make_map(image)
;nc=100 
;loadct,3,ncolors=nc
;plot_map,submap,ncolors=nc
;alpha1=0.5
;

END
