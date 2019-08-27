import numpy as np
from mpi4py import MPI
import os,sys
import glob
#from surya.mwa import mwa_fluxcal as mfl


def __MAIN__():
	__MPI_MAIN__()


def STEP_1_MASTER_PROC_READ_INPUT(DR_DIR,WORKING_DIR):
	os.chdir(DS_DIR)
	#DS_LIST = list(set(sorted(glob.glob('1*.DS.dat.p'))) - set(sorted(glob.glob('1*%b084-085_T000-008.DS.dat.p'))))
	#DS_LIST=sorted(glob.glob('1079445*-%b229-232_T063-089.DS.dat.p'))
	DS_LIST = sorted(glob.glob('1*.p'))
	print 'Total files :',len(DS_LIST)
	os.chdir(WORKING_DIR)
	return DS_LIST	


def MPI_MAIN_FLUX(mfl,DS_LIST,DS_DIR,BEAM_DIR,WORKING_DIR,HASLAM_DIR,mwa_phase,array_lon,array_lat,array_elevation,rec_path,grd_path,spidx_path,ifsun,star_ra,star_dec):
	comm = MPI.COMM_WORLD
	rank = comm.Get_rank()
	stat = MPI.Status()
	num_of_processor = comm.Get_size()
	if rank == 0:
		#DS_LIST=STEP_1_MASTER_PROC_READ_INPUT()
		DS_LIST_DIV=np.array_split(np.array(DS_LIST),num_of_processor-1)
		# Send data to slave processors
		for indx in xrange(1,num_of_processor):
			comm.send(DS_LIST_DIV[indx-1], dest=indx, tag=indx)
		#comm.Barrier()
		# Recieve Data from slave Processors
		#comm.Barrier()
		
	else:
		# Recieve Data from Master Processor
		DS_LIST	= comm.recv(source=0, tag=rank)
		#comm.Barrier()
		print 'STARTING RANK: ',rank
		for DS in DS_LIST:
			mfl.tsun_computation(DS,DS_DIR,BEAM_DIR,WORKING_DIR,HASLAM_DIR,mwa_phase,array_lon,array_lat,array_elevation,rec_path,grd_path,spidx_path,ifsun,star_ra,star_dec)
			#tsun_computation_sources(DS)
		print 'FINISHING RANK: ',rank
		#comm.Barrier()
		# Send data to Master Processor
		
		
		
