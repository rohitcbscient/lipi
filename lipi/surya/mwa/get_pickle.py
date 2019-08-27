import numpy as np
import glob
import matplotlib.pyplot as plt
import commands
import pickle
import os
import urllib
import urllib2
import json
import sqlite3 as sq
import time
import os
time1= time.time()

##############################################################################################################################################
def get_freq(MSNAME,tbo):			# Return the start frequency and channel width of the first spectral channel in the data set
	tbo.open(MSNAME+'/SPECTRAL_WINDOW')
	RefFreq = tbo.getcol('REF_FREQUENCY')[0]
	RefFreq = RefFreq/(10**6)		# Conversion from Hz to MHz
	spec_resol = tbo.getcol('CHAN_WIDTH')[0][0]/(10**3)
	tbo.close()
	return RefFreq,spec_resol

def calc_central_freq(RefFreq,chan_list,spec_resol):	#Calculate central frequency for each spectral channel in the channel list
	central_freq = np.zeros(len(chan_list))
	for i in range(len(chan_list)):
		central_freq[i] = RefFreq + (chan_list[i]+1)*spec_resol/(1000)
	return central_freq

def get_amp(MSNAME,tile1,tile2,POL_LIST,mso,casalogo):		# Extraction of amplitude of crosscorrelations
	mso.open(MSNAME)
        casalogo.post("#### Tile1 %03d; Tile2 %03d" % (tile1, tile2));
        print '#### Tile1 %03d; Tile2 %03d; POL_LIST %s' % (tile1, tile2, POL_LIST)
        mso.selectinit(datadescid=0) # Reset any earlier selections
        mso.select({'antenna1':[tile1],'antenna2':[tile2]})
        mso.selectpolarization(POL_LIST)
	#amp=ms.getdata(['amplitude'])
	amp=mso.getdata(["amplitude","axis_info"],ifraxis=True)
	mso.close()
	return amp['amplitude']

def get_phases(MSNAME,tile1,tile2,POL_LIST,mso,casalogo):	# Extraction of phase of crosscorrelations data along with u,v and w for each time slice.
        mso.open(MSNAME)
        casalogo.post("#### Tile1 %03d; Tile2 %03d" % (tile1, tile2));
        print '#### Tile1 %03d; Tile2 %03d; POL_LIST %s' % (tile1, tile2, POL_LIST)
        mso.selectinit(datadescid=0) # Reset any earlier selections
        mso.select({'antenna1':[tile1],'antenna2':[tile2]})
        mso.selectpolarization(POL_LIST)
        phases=mso.getdata(['phase'])['phase']
        u=mso.getdata(['u'])['u']
        v=mso.getdata(['v'])['v']
        w=mso.getdata(['w'])['w']
	mso.close()
        return phases,u,v,w

def calc_ncross(auto_t1,auto_t2,cross):		# Calculation of normalised crosscorelations
	ncross = [None]*4
	j = 0
	k = 0		#Counter variables for auto_t1 and auto_t2
	for i in range(4):		
		ncross[i] = cross[i]/np.sqrt(auto_t1[j]*auto_t2[k])
		if (i%2==0):
			k+=1
		else:
			j+=1
			k = 0
	ncross = np.array(ncross)
	return ncross

# Function to call a JSON web service and return a dictionary:
def get_metafit1(obs_id):
        BASEURL = 'http://mwa-metadata01.pawsey.org.au/metadata/obs/'
        obs_id = obs_id
        params = {'obs_id': obs_id}
        result = urllib2.urlopen(url = BASEURL, data = urllib.urlencode(params))
        result_set = json.loads(result.read())
        azi_pointing = result_set['metadata']['azimuth_pointing']
        ele_pointing = result_set['metadata']['elevation_pointing']
        gridpoint_number = result_set['metadata']['gridpoint_number']
        print azi_pointing, ele_pointing, gridpoint_number
        return azi_pointing, ele_pointing, gridpoint_number

def getmeta(BASEURL,service='obs', params=None):
  """Given a JSON web service ('obs', find, or 'con') and a set of parameters as
     a Python dictionary, return a Python dictionary containing the result.
  """
  if params:
    data = urllib.urlencode(params)  # Turn the dictionary into a string with encoded 'name=value' pairs
  else:
    data = ''
  #             Validate the service name
  if service.strip().lower() in ['obs', 'find', 'con']:
    service = service.strip().lower()
  else:
    print "invalid service name: %s" % service
    return
  #             Get the data
  try:
    result = json.load(urllib2.urlopen(BASEURL + service + '?' + data))
  except urllib2.HTTPError as error:
    print "HTTP error from server: code=%d, response:\n %s" % (error.code, error.read())
    return
  except urllib2.URLError as error:
    print "URL or network error: %s" % error.reason
    return
  #            Return the result dictionary
  return result

def update(c,tablename,MSNAME,obs_id,baseline,starttime,stoptime,RefFreq,spec_resol,int_time,u_mean,v_mean,w_mean,azi_pointing,ele_pointing,gridpoint_no):
	names = c.execute('SELECT MSNAME FROM %s'% (tablename))
	names = names.fetchall()		#names returns list of filenames. Each entry of list is of form (unicode(MSNAME),)
	baselines = c.execute('SELECT Baseline FROM %s'% (tablename))
	baselines = baselines.fetchall()
	stopfreq = RefFreq + 30.72		#30.72MHz is the bandwidth of data in one MS file.
	row = [MSNAME,obs_id,baseline,str(starttime),str(stoptime),str(RefFreq),str(stopfreq),str(spec_resol),str(int_time),str(u_mean),str(v_mean),str(w_mean),str(azi_pointing),str(ele_pointing),str(gridpoint_no)]
	c.execute('INSERT INTO %s VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)'% (tablename),row)

#print ms.__doc__
def datetime(MSNAME,mso,qao):		# Returning observation integration time (in seconds),start date,start time, end date and end time (in UT)
	mso.open(MSNAME)
	a= mso.summary()
	start=(qao.time(qao.quantity(a['BeginTime'],'d'),form="ymd"))[0]
	end=(qao.time(qao.quantity(a['EndTime'],'d'),form="ymd"))[0]
	int_time = a['scan_1']['0']['IntegrationTime']

	starttime = start.split('/')[0] + '/' + start.split('/')[1] + '/' + start.split('/')[2] + ' ' + start.split('/')[3]
	endtime = end.split('/')[0] + '/'+ end.split('/')[1] + '/' + end.split('/')[2] + ' ' + end.split('/')[3]
	mso.close()	
	return int_time,starttime,endtime

def get_TileID_phase1(TileName):
	"""
	From MWA PHASE-I. 
	A dictionary to return the TileID given the TileName as specified 
	in the MS ('Tile???MWA') Use this function rather than using the IDs 
	directly simply to reduce the possibility of human error.
	
	Divya 02Dec2015
	"""
	TileDict = {'Tile011MWA' : 0, 'Tile012MWA' : 1, 'Tile013MWA' : 2, 'Tile014MWA' : 3, 'Tile015MWA' : 4, 'Tile016MWA' : 5, 'Tile017MWA' : 6, 'Tile018MWA' : 7, 'Tile021MWA' : 8, 'Tile022MWA' : 9, 'Tile023MWA' : 10, 'Tile024MWA' : 11, 'Tile025MWA' : 12, 'Tile026MWA' : 13, 'Tile027MWA' : 14, 'Tile028MWA' : 15, 'Tile031MWA' : 16, 'Tile033MWA' : 18, 'Tile034MWA' : 19, 'Tile035MWA' : 20, 'Tile036MWA' : 21, 'Tile037MWA' : 22, 'Tile038MWA' : 23, 'Tile041MWA' : 24, 'Tile042MWA' : 25, 'Tile043MWA' : 26, 'Tile044MWA' : 27, 'Tile045MWA' : 28, 'Tile046MWA' : 29, 'Tile047MWA' : 30, 'Tile048MWA' : 31, 'Tile051MWA' : 32, 'Tile052MWA' : 33, 'Tile053MWA' : 34, 'Tile054MWA' : 35, 'Tile055MWA' : 36, 'Tile056MWA' : 37, 'Tile057MWA' : 38, 'Tile058MWA' : 39, 'Tile061MWA' : 40, 'Tile062MWA' : 41, 'Tile063MWA' : 42, 'Tile064MWA' : 43, 'Tile065MWA' : 44, 'Tile066MWA' : 45, 'Tile067MWA' : 46, 'Tile068MWA' : 47, 'Tile071MWA' : 48, 'Tile072MWA' : 49, 'Tile073MWA' : 50, 'Tile074MWA' : 51, 'Tile075MWA' : 52, 'Tile076MWA' : 53, 'Tile077MWA' : 54, 'Tile078MWA' : 55, 'Tile081MWA' : 56, 'Tile082MWA' : 57, 'Tile083MWA' : 58, 'Tile084MWA' : 59, 'Tile085MWA' : 60, 'Tile086MWA' : 61, 'Tile087MWA' : 62, 'Tile088MWA' : 63, 'Tile091MWA' : 64, 'Tile092MWA' : 65, 'Tile093MWA' : 66, 'Tile094MWA' : 67, 'Tile095MWA' : 68, 'Tile096MWA' : 69, 'Tile097MWA' : 70, 'Tile098MWA' : 71, 'Tile101MWA' : 72, 'Tile102MWA' : 73, 'Tile103MWA' : 74, 'Tile104MWA' : 75, 'Tile105MWA' : 76, 'Tile106MWA' : 77, 'Tile107MWA' : 78, 'Tile108MWA' : 79, 'Tile112MWA' : 81, 'Tile113MWA' : 82, 'Tile114MWA' : 83, 'Tile115MWA' : 84, 'Tile116MWA' : 85, 'Tile117MWA' : 86, 'Tile118MWA' : 87, 'Tile121MWA' : 88, 'Tile122MWA' : 89, 'Tile123MWA' : 90, 'Tile124MWA' : 91, 'Tile125MWA' : 92, 'Tile126MWA' : 93, 'Tile127MWA' : 94, 'Tile128MWA' : 95, 'Tile131MWA' : 96, 'Tile132MWA' : 97, 'Tile133MWA' : 98, 'Tile134MWA' : 99, 'Tile135MWA' : 100, 'Tile136MWA' : 101, 'Tile137MWA' : 102, 'Tile138MWA' : 103, 'Tile141MWA' : 104, 'Tile142MWA' : 105, 'Tile143MWA' : 106, 'Tile144MWA' : 107, 'Tile145MWA' : 108, 'Tile146MWA' : 109, 'Tile147MWA' : 110, 'Tile148MWA' : 111, 'Tile151MWA' : 112, 'Tile152MWA' : 113, 'Tile153MWA' : 114, 'Tile154MWA' : 115, 'Tile155MWA' : 116, 'Tile156MWA' : 117, 'Tile157MWA' : 118, 'Tile158MWA' : 119, 'Tile161MWA' : 120, 'Tile162MWA' : 121, 'Tile163MWA' : 122, 'Tile164MWA' : 123, 'Tile165MWA' : 124, 'Tile166MWA' : 125, 'Tile167MWA' : 126, 'Tile168MWA' : 127}
	return TileDict[TileName]

def get_TileID_phase2(TileName):
	"""
	From MWA PHASE-II. 
	A dictionary to return the TileID given the TileName as specified 
	in the MS ('Tile???MWA') Use this function rather than using the IDs 
	directly simply to reduce the possibility of human error.
	
	Rohit 
	"""
	TileDict = {'Tile051MWA' : 0, 'Tile052MWA' : 1, 'Tile053MWA' : 2, 'Tile054MWA' : 3, 'Tile055MWA' : 4, 'Tile056MWA' : 5, 'Tile057MWA' : 6, 'Tile058MWA' : 7, 'Tile071MWA' : 8, 'Tile072MWA' : 9, 'Tile073MWA' : 10, 'Tile074MWA' : 11, 'Tile075MWA' : 12, 'Tile076MWA' : 13, 'Tile077MWA' : 14, 'Tile078MWA' : 15, 'Tile101MWA' : 16, 'Tile102MWA' : 17 ,'Tile103MWA' : 18, 'Tile104MWA' : 19, 'Tile105MWA' : 20, 'Tile106MWA' : 21, 'Tile107MWA' : 22, 'Tile108MWA' : 23, 'Tile111MWA' : 24, 'Tile112MWA' : 25, 'Tile113MWA' : 26, 'Tile114MWA' : 27, 'Tile115MWA' : 28, 'Tile116MWA' : 29, 'Tile117MWA' : 30, 'Tile118MWA' : 31, 'Tile121MWA' : 32, 'Tile122MWA' : 33, 'Tile123MWA' : 34, 'Tile124MWA' : 35, 'Tile125MWA' : 36, 'Tile126MWA' : 37, 'Tile127MWA' : 38, 'Tile128MWA' : 39, 'Tile131MWA' : 40, 'Tile132MWA' : 41, 'Tile133MWA' : 42, 'Tile134MWA' : 43, 'Tile135MWA' : 44, 'Tile136MWA' : 45, 'Tile137MWA' : 46, 'Tile138MWA' : 47, 'Tile141MWA' : 48, 'Tile142MWA' : 49, 'Tile143MWA' : 50, 'Tile144MWA' : 51, 'Tile145MWA' : 52, 'Tile146MWA' : 53, 'Tile147MWA' : 54, 'Tile148MWA' : 55, 'Tile151MWA' : 56, 'Tile152MWA' : 57, 'Tile153MWA' : 58, 'Tile154MWA' : 59, 'Tile155MWA' : 60, 'Tile156MWA' : 61, 'Tile157MWA' : 62, 'Tile158MWA' : 63, 'Tile161MWA' : 64, 'Tile162MWA' : 65, 'Tile163MWA' : 66, 'Tile164MWA' : 67, 'Tile165MWA' : 68, 'Tile166MWA' : 69, 'Tile167MWA' : 70, 'Tile168MWA' : 71, 'LBA1MWA' : 72, 'LBA2MWA' : 73, 'LBA3MWA' : 74, 'LBA4MWA' : 75, 'LBA5MWA' : 76, 'LBA6MWA' : 77, 'LBA7MWA' : 78, 'LBA8MWA' : 79, 'LBB1MWA' : 80, 'LBB2MWA' : 81, 'LBB3MWA':82, 'LBB4MWA': 83, 'LBB5MWA' : 84, 'LBB6MWA' : 85, 'LBB7MWA' : 86, 'LBB8MWA' : 87, 'LBC1MWA' : 88, 'LBC2MWA' : 89, 'LBC3MWA' : 90, 'LBC4MWA' : 91, 'LBC5MWA' : 92, 'LBC6MWA' : 93, 'LBC7MWA' : 94, 'LBC8MWA' : 95, 'LBD1MWA' : 96, 'LBD2MWA' : 97, 'LBD3MWA' : 98, 'LBD4MWA' : 99, 'LBD5MWA' : 100, 'LBD6MWA' : 101, 'LBD7MWA' : 102, 'LBD8MWA' : 103, 'LBE1MWA' : 104, 'LBE2MWA' : 105, 'LBE3MWA' : 106, 'LBE4MWA' : 107, 'LBE5MWA' : 108, 'LBE6MWA' : 109, 'LBE7MWA' : 110, 'LBE8MWA' : 111, 'LBF1MWA' : 112, 'LBF2MWA' : 113, 'LBF3MWA' : 114, 'LBF4MWA' : 115, 'LBF5MWA' : 116, 'LBF6MWA' : 117, 'LBF7MWA' : 118, 'LBF8MWA' : 119, 'LBG1MWA' : 120, 'LBG2MWA' : 121, 'LBG3MWA' : 122, 'LBG4MWA' : 123, 'LBG5MWA' : 124, 'LBG6MWA' : 125, 'LBG7MWA' : 126, 'LBG8MWA' : 127}
	return TileDict[TileName]
def baseline_from_tile(tilenames_list,phase):		# Generate list of baselines given the tile names
	baseline_list = []
	tile1 = []
	tile2 = []
	for i in range(len(tilenames_list)): 
		for j in range(i+1,len(tilenames_list)):
			if (phase==1):
				tile1.append(get_TileID_phase1(tilenames_list[i]))
				tile2.append(get_TileID_phase1(tilenames_list[j]))
			if (phase==2):
				tile1.append(get_TileID_phase2(tilenames_list[i]))
                                tile2.append(get_TileID_phase2(tilenames_list[j]))
			baseline_list.append(tilenames_list[i] + '-' + tilenames_list[j])
	tile1 = np.array(tile1)
	tile2 = np.array(tile2)
	baseline_list = np.array(baseline_list)	
	return tile1,tile2,baseline_list


#conn = sq.connect(database)	 	#Creating a connection to the database 
#c = conn.cursor()			#Placing cursor in database for execution of commands
#c.execute('DROP TABLE %s'% (tablename))
#c.execute('CREATE TABLE %s (MSNAME varchar2, Obs_id varchar2, Baseline varchar2, Starttime varchar2, Stoptime  varchar2, Startfreq  varchar2, Stopfreq varchar2, Freq_res  varchar2, Int_time  varchar2, u_mean varchar2, v_mean varchar2, w_mean varchar2,Azimuth_pointing varchar2, Elevation_pointing varchar2,Gridpoint_number varchar2)'% (tablename) )	#Needed to create database for first time	(Assuming that we are dealing with only 4-5 mins of data)

def main(INDATA_DIR,OUTDATA_DIR,DB_DIR,WORKING_DIR,MS_LIST,POL_LIST,CPOL_LIST,MWA_PHASE,azi_pointing,ele_pointing,tilenames_list,ms,qa,tb,casalog):
    mso=ms
    qao=qa
    tbo=tb
    casalogo=casalog
    tile1,tile2,baseline_list = baseline_from_tile(tilenames_list,MWA_PHASE)
    os.chdir(INDATA_DIR)
    #MS_LIST = sorted(glob.glob('1208325736-b187-188.MS'))
    #MS_LIST=['1208326216-b062-063.MS','1208326216-b187-188.MS']
    print MS_LIST
    os.chdir(DB_DIR)
    auto_t1 = [[0]*len(baseline_list)]*len(MS_LIST)
    auto_t2 = [[0]*len(baseline_list)]*len(MS_LIST)
    cross_t12 = [[0]*len(baseline_list)]*len(MS_LIST)
    u = [[0]*len(baseline_list)]*len(MS_LIST)
    v = [[0]*len(baseline_list)]*len(MS_LIST)
    w = [[0]*len(baseline_list)]*len(MS_LIST)
    phase_ncross = [[0]*len(baseline_list)]*len(MS_LIST)
    os.chdir(INDATA_DIR)
    for i in range(len(MS_LIST)):
            os.chdir(INDATA_DIR)	
            MSNAME = MS_LIST[i]
            int_time,start_time,end_time = datetime(MSNAME,mso,qao)
            RefFreq,spec_resol = get_freq(MSNAME,tbo)
            #obs_id=MSNAME.split('.')[0]
            obs_id=MSNAME.split('-')[0]
            print "Observation id : \t",obs_id
            #obsinfo = getmeta(BASEURL,service='obs', params={'obs_id':obs_id})
            #azi_pointing,ele_pointing,gridpoint_number = get_metafit1(obs_id)
            ph_ra,ph_dec=azi_pointing,ele_pointing
            for j in range(len(baseline_list)):		
                    print "Autocorrelations"
                    auto_t1[i][j] = get_amp(MS_LIST[i],tile1[j],tile1[j],POL_LIST,mso,casalogo)			# Autocorrelation from tile 1 - XX and YY pol.
                    print auto_t1[i][j].shape
                    auto_t2[i][j] = get_amp(MS_LIST[i],tile2[j],tile2[j],POL_LIST,mso,casalogo)			# Autocorrelation from tile 2	
                    print "Crosscorrelations.... size: ", auto_t2[i][j].shape
                    cross_t12[i][j] = get_amp(MS_LIST[i],tile1[j],tile2[j],CPOL_LIST,mso,casalogo)	# Crosscorrelations - XX,XY,YX and YY pols. in sequence
                    phase_ncross[i][j],u[i][j],v[i][j],w[i][j] = get_phases(MS_LIST[i],tile1[j],tile2[j],CPOL_LIST,mso,casalogo)
                    print "MSNAME :",MSNAME,"\nTile 1 :",tile1[j],"\nTile 2 :",tile2[j]
                    chan_list = np.arange(len(auto_t1[i][j][0]))		# List of all channels - can be refined to select only "good" channels
                    central_freq = calc_central_freq(RefFreq,chan_list,spec_resol)
                    ncross = calc_ncross(auto_t1[i][j],auto_t2[i][j],cross_t12[i][j])		
                    print ncross
                    DS_FILE = MSNAME.split('.MS')[0] + '_T'+ '%03d'% (tile1[j])+'-'+'%03d'% (tile2[j])+'.DS.dat.p'
                    #update(c,tablename,MSNAME,obs_id,baseline_list[j],start_time,end_time,RefFreq,spec_resol,int_time,np.mean(u[i][j]),np.mean(v[i][j]),np.mean(w[i][j]),azi_pointing,ele_pointing,gridpoint_number)
                    os.chdir(OUTDATA_DIR)		
                    pickle.dump([central_freq,chan_list,auto_t1[i][j],auto_t2[i][j],cross_t12[i][j],ncross,phase_ncross[i][j],u[i][j],v[i][j],w[i][j],azi_pointing,ele_pointing,start_time,end_time],open(DS_FILE,"wb"))
                    os.chdir(DB_DIR)
    os.chdir(WORKING_DIR)
##############################################################################################################################################
