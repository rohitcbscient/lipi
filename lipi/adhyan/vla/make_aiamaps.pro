
filename='ssw_cutout_20160409_184022_AIA_171_.fts' 
read_sdo, filename, index, data
;aia_respike, index, data, out_index, out_data
;aia_prep, out_index, out_data, final_index, final_data, /use_pnt_file
aia_prep, index, data, final_index, final_data, /use_pnt_file
writefits,'test.fits',out_data,out_index
end
