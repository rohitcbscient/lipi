import matplotlib.pyplot as plt



def flagcmd_tf(vis_,infile):
    '''
    Telescope flags 
    '''
    flagcmd(vis=vis_,inpmode='list',inpfile=infile)
    

def flagscan(vis_,scan0):
    '''
    Initial flags of the system configuration and first scans of Calibration and Sun 
    '''
    flagdata(vis=vis_,scan=scan0) 
    
def flagquack(vis_,quackinterval_):
    '''
    Flag quackinterval (time) from the start of the scan
    '''
    flagdata(vis=vis_,mode='quack',quackmode='beg',quackinterval=quackinterval_)

    
        

