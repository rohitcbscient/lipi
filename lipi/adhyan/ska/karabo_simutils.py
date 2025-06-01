import numpy as np
from astropy.io import fits
from karabo.imaging.imager_oskar import OskarDirtyImager, OskarDirtyImagerConfig
from karabo.imaging import imager_wsclean

def get_end_freq(start_freq,bandwidth,samples):
    """Get end frequency of the radio band
    Caution:Give freq and bandwidth arguments in the same units
    Args:
        start_freq (float): Start Frequency
        bandwidth (float): End Frequency
        samples (int): number of samples
    """
    end_freq = start_freq+bandwidth
    freq_arr = np.linspace(start_freq,end_freq,samples)
    return end_freq,freq_arr

def get_telescope(telpath):
    """Get the SKA and it's precursor configurations

    Paramaters:
        telpath (str): Path to the telescope repo

    Returns:
        telescope_ska_combined: Paths for SKA-mid and SKA-low
        telescope_precursor_combined: Paths for SKA precursor 
    """
    telescope_meerkat = telpath+'meerkat.tm'
    telescope_mwaphase1 = telpath+'mwa.phase1.tm'
    telescope_ska1mid = telpath+'ska1mid.tm'
    telescope_ska1midAAs = telpath+'ska-mid-AAstar.tm'
    telescope_ska1midAA2 = telpath+'ska-mid-AA2.tm'
    telescope_ska1midAA1 = telpath+'ska-mid-AA1.tm'
    telescope_ska1midAA05 = telpath+'ska-mid-AA0.5.tm'
    telescope_ska1low = telpath+'ska1low.tm'
    telescope_ska1lowAAs = telpath+'ska-low-AAstar.tm'
    telescope_ska1lowAA2 = telpath+'ska-low-AA2.tm'
    telescope_ska1lowAA1 = telpath+'ska-low-AA1.tm'
    telescope_ska1lowAA05 = telpath+'ska-low-AA0.5.tm'
    telescope_ska_combined = {'skamid':[telescope_ska1mid,telescope_ska1midAAs,telescope_ska1midAA2,\
                            telescope_ska1midAA1,telescope_ska1midAA05],
                            'skalow':[telescope_ska1low,telescope_ska1lowAAs,telescope_ska1lowAA2,\
                            telescope_ska1lowAA1,telescope_ska1lowAA05]
                            }
    telescope_precursor_combined = {'mwa':[telescope_mwaphase1],'meerkat':[telescope_meerkat]}
    return telescope_ska_combined,telescope_precursor_combined


def get_telescope_config(telescope,info):
    """Function gives the relavent information about the telescope config
    Parameters
    ----------
        'telescope_name'(str): Path of the .tm file
    Returns
    -------
    `base_length` (array): Array of baselines
    `maxb` (float): Maximum baseline
    `medb` (float): Median baseline
    """
    layout=np.loadtxt(telescope+'/layout.txt')
    nant=len(layout)
    nb=int(nant*(nant-1)*0.5)
    base_length=[0]*nb
    k=0
    for i in range(nant):
        for j in range(i,nant):
            if(i!=j):
                base_length[k] = np.sqrt((layout[i][0]-layout[j][0])**2 + (layout[i][1]-layout[j][1])**2)
                k=k+1
    base_length=np.array(base_length); maxb = np.nanmax(base_length); medb = np.nanmedian(base_length)
    if(info):
        print('---------------------')
        print(telescope)
        print("Total Baselines: ",nb)
        print('Maximum Baseline (m):',np.round(maxb,2),'| Median Baseline (m):',np.round(medb,2))
    return base_length,maxb,medb

def ska_frequency_params(bandwidth_samples,info=None):
    """SKA Parameters 
    Note that all values are in MHz
    Parameters
    ----------
        `baseline_samples`: Number of samples for dividing the bandwidth
    Return
    ------
        ``
        ska1low,ska1mid
    """
    start_freq_ska1_low = 50
    start_freq_ska1_mid1 = 350
    start_freq_ska1_mid2 = 950
    start_freq_ska1_mid3 = 1650
    start_freq_ska1_mid4 = 2800
    start_freq_ska1_mid5a = 4600
    start_freq_ska1_mid5b = 8300
    #---------------------------
    bandwidth_ska1_low = 300
    bandwidth_ska1_mid1 = 700
    bandwidth_ska1_mid2 = 810
    bandwidth_ska1_mid3 = 1400
    bandwidth_ska1_mid4 = 2380
    bandwidth_ska1_mid5a = 3900
    bandwidth_ska1_mid5b = 5000
    #--------------------------
    channel_width_ska1_low = 5.4e-3 
    channel_width_ska1_mid1 = 5.4e-3
    channel_width_ska1_mid2 = 13.4e-3
    channel_width_ska1_mid3 = 13.4e-3
    channel_width_ska1_mid4 = 13.4e-3
    channel_width_ska1_mid5a = 80.6e-3
    channel_width_ska1_mid5b = 80.4e-3
    #---------------------------------
    end_freq_ska1_low, freq_ska1_low = get_end_freq(start_freq_ska1_low,bandwidth_ska1_low,bandwidth_samples)
    end_freq_ska1_mid1, freq_ska1_mid1 = get_end_freq(start_freq_ska1_mid1,bandwidth_ska1_mid1,bandwidth_samples)
    end_freq_ska1_mid2, freq_ska1_mid2 = get_end_freq(start_freq_ska1_mid2,bandwidth_ska1_mid2,bandwidth_samples)
    end_freq_ska1_mid3, freq_ska1_mid3 = get_end_freq(start_freq_ska1_mid3,bandwidth_ska1_mid3,bandwidth_samples)
    end_freq_ska1_mid4, freq_ska1_mid4 = get_end_freq(start_freq_ska1_mid4,bandwidth_ska1_mid4,bandwidth_samples)
    end_freq_ska1_mid5a, freq_ska1_mid5a = get_end_freq(start_freq_ska1_mid5a,bandwidth_ska1_mid5a,bandwidth_samples)
    end_freq_ska1_mid5b, freq_ska1_mid5b = get_end_freq(start_freq_ska1_mid5b,bandwidth_ska1_mid5b,bandwidth_samples)
    #----------------------------------------------------------------------------------------------------------------
    tellow_param = [[freq_ska1_low,start_freq_ska1_low,bandwidth_ska1_low,end_freq_ska1_low,channel_width_ska1_low]]
    freq_ska1_mid = np.vstack((freq_ska1_mid1,freq_ska1_mid2,freq_ska1_mid3,freq_ska1_mid4,freq_ska1_mid5a,freq_ska1_mid5b))
    start_freq_ska1_mid = np.vstack((start_freq_ska1_mid1,start_freq_ska1_mid2,start_freq_ska1_mid3,start_freq_ska1_mid4,start_freq_ska1_mid5a,start_freq_ska1_mid5b))
    end_freq_ska1_mid = np.vstack((end_freq_ska1_mid1,end_freq_ska1_mid2,end_freq_ska1_mid3,end_freq_ska1_mid4,end_freq_ska1_mid5a,end_freq_ska1_mid5b))
    bandwidth_ska1_mid = np.vstack((bandwidth_ska1_mid1,bandwidth_ska1_mid2,bandwidth_ska1_mid3,bandwidth_ska1_mid4,bandwidth_ska1_mid5a,bandwidth_ska1_mid5b))
    channel_width_ska1_mid = np.vstack((channel_width_ska1_mid1,channel_width_ska1_mid2,channel_width_ska1_mid3,channel_width_ska1_mid4,channel_width_ska1_mid5a,channel_width_ska1_mid5b))
    telmid_param = [freq_ska1_mid,start_freq_ska1_mid,bandwidth_ska1_mid,end_freq_ska1_mid,channel_width_ska1_mid]
    tel_params = {'low':tellow_param,'mid':telmid_param}
    if(info):
        print('---- Array Information (MHz) ---------')
        print('(1) frequency array, (2) start frequency, '\
            '(3) bandwidth, (4) end frequency, (5) uniform channel width')
    return tel_params


def get_solar_skydata(skymodel_path,sm_save_str,start_frequency_hz_,ra_sun_center,dec_sun_center,\
                skymodel_cellsize,sm_fscale,sm_threshold, save_sources,source_cut,\
                add_sm, randfile, randsize, point_flux, flux_rand_low, flux_rand_high,angle_rand):
    """Get skymodel data with optional features of adding point sources

    Parameters:
        skymodel_path (str): path of the skymodel
        sm_save_str (str): path to save skymodel
        start_frequency_hz_ (float): start frequency in Hz
        ra_sun_center (float): ra of the center in degree
        dec_sun_center (float): dec of the center in degree
        sm_fscale (float): spectral index to scale flux density with frequency
        sm_threshold (float): cut in flux density of skymodel w.r.t maximum flux density
        save_sources (bool): save skymodel added sources (not entire skymodel)
        source_cut (int): cut on the number of source
        randfile (str): Path to file cointaining random sources 
        randsize (int): size of randomly (uniform) added sources
        point_flux (float): flux density of point source
        flux_rand_low (float): lower end of the random source distribution
        flux_rand_high (float): higher end of the random source distribution
        angle_rand (float): angular extent of the random source distribution w.r.t solar center 

    Returns:
        `sky_array` (array): Array of the skymodel (disk + added sources)
        `solar_map` (image): 2-D array of the disk skymodel at reference frequency (240 MHz)
        `solar_map_jy` (image): Frequency scaled flux density map in Jansky
    """
    solar_map,smh=fits.getdata(skymodel_path,header=True)
    solar_map_jy=solar_map/np.nanmax(solar_map)\
                    *20*1.e4*(start_frequency_hz_/2.4e8)**sm_fscale
    ra_grid,dec_grid=np.meshgrid((np.arange(smh['NAXIS1'])-smh['NAXIS1']/2)*skymodel_cellsize/3600.,
                                (np.arange(smh['NAXIS2'])-smh['NAXIS2']/2)*skymodel_cellsize/3600.)
    ra_grid=ra_grid+ra_sun_center;dec_grid=dec_grid+dec_sun_center
    idx=np.where(solar_map>sm_threshold*np.nanmax(solar_map))
    sky_model_ra=ra_grid[idx];sky_model_dec=dec_grid[idx];flux=solar_map_jy[idx]
    sky_data = np.array([sky_model_ra, sky_model_dec, flux,np.zeros(len(flux)), \
            np.zeros(len(flux)),np.zeros(len(flux)), np.zeros(len(flux)),np.zeros(len(flux)), \
        np.zeros(len(flux)),np.zeros(len(flux)), np.zeros(len(flux)),np.zeros(len(flux))]).T
    sky_data=sky_data[0:source_cut,:]
    add_source = np.zeros((1,12))
    if(add_sm=='point'):
        size_=1
        add_source = np.zeros((size_,12))
        add_source[0][0] = ra_sun_center # RA
        add_source[0][1] = dec_sun_center # DEC
        add_source[0][2] = point_flux # Flux
        sky_data1 = np.vstack((sky_data,add_source))
    elif(add_sm=='random'):
        if(randfile==''):
            size_=randsize
            add_source = np.zeros((size_,12))
            ra_point_array = np.random.uniform(low=ra_sun_center-angle_rand, high=ra_sun_center+angle_rand, size=size_)
            dec_point_array = np.random.uniform(low=dec_sun_center-angle_rand, high=dec_sun_center+angle_rand, size=size_)
            flux_array = np.random.uniform(low=flux_rand_low, high=flux_rand_high, size=size_)
            add_source[:,0] = ra_point_array # RA
            add_source[:,1] = dec_point_array # DEC
            add_source[:,2] = flux_array # Flux
            sky_data1 = np.vstack((sky_data,add_source))
        else:
            sky_data = np.loadtxt(randfile)
            sky_data1 = sky_data
    else:
        sky_data1 = sky_data
    if(randfile==''):
        save_sources = False
    if(save_sources):
        np.savetxt(sm_save_str+'.sm', add_source.T)
    return sky_data1,solar_map,solar_map_jy

def make_images(vis,imager_str,imgsize,img_str,cellsize_rad,niter,weight,maxuv,minuv,nchan,vis_str):
    """Make either `oskar` dirty or `wsclean` images

    Parameters:
        vis (object): visibility instance
        imager_str (_type_): _description_
        imgsize (_type_): _description_
        img_str (_type_): _description_
        cellsize_rad (_type_): _description_
        niter (_type_): _description_
        weight (_type_): _description_
        maxuv (_type_): _description_
        minuv (_type_): _description_
        nchan (_type_): _description_
        vis_str (_type_): _description_

    Returns:
        image (image_instance): _description_
    """
    dirty_image='no_image';restored='no_image'
    if(imager_str=='oskar'):
        dirty_imager = OskarDirtyImager(
        OskarDirtyImagerConfig(
            imaging_npixel=imgsize,
            imaging_cellsize=cellsize_rad,
        ))
        dirty_image = dirty_imager.create_dirty_image(vis,output_fits_path=img_str+'oskar.fits')
    if(imager_str=='wsclean'):
        img_cmd = 'wsclean \
                -size '+str(imgsize)+' '+str(imgsize)+' \
                -name '+img_str+' \
                -scale '+str(cellsize_rad)+'rad -niter '+str(niter)+' -mgain 0.8 \
                -weight '+weight+'\
                -maxuv-l '+str(maxuv)+' -minuv-l '+str(minuv)+'\
                -channels-out '+str(nchan)+' '+vis_str+'.ms'
        print(img_cmd)
        try:
            restored = imager_wsclean.create_image_custom_command(command=img_cmd)
        except:
            pass # doing nothing on exception
    return restored,dirty_image