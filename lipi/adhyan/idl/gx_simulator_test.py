import pidly
idl = pidly.IDL('/sdata/ssw/gen/setup/ssw_idl') # Path to the sswidl

fitsfiles=idl('fitsfiles=findfile("/sdata/20140914_hmi/hmi/hmi_20140914_*_magnetogram.fits")')
idl('read_sdo,fitsfiles(0),index,data')

# To run entire file: idl('.r ~/my_git/lipi/adhyan/mwa/20140914_aia_analysis.pro')


