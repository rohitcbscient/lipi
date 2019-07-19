# This program construct the parameter array for the input to the gsfit

def gsfit_var(A,d,T0,eps,kappa,nn,emin,emax,ebrk,del1,del2,nth,nb,B,theta,fs,delf,ed,nf,paf,lca,bda,delmu,a4,fc,fw,mat,q):
    '''
    INPUTS
    A: Area in cm^-2
    d: Depth in cm
    T0: Temperature in K
    eps: 
    kappa: 
    nn: Number of integration nodes in energy
    emin: E_min in MeV
    emax: E_max in MeV
    ebrk: E_break in MeV 
    del1: 
    del2: 
    nth: thermal electron density cm^-3
    B: Mangetic field 
    theta: Viewing angle in degrees
    fs: starting frequencies to calculate spectrum in Hz
    delf: logarithmic step in frequency
    ed: distribution over energy
    nf: Number of frequencies
    paf: distribution over pitch-angle 
    lca: loss cone boundary in degrees
    bda: beam direction in GAU and SGA
    delmu: delta mu
    a4: a_4 in SGA
    fc: 
    fw:
    mat: mathcing on
    q: Q-optimization on

    NOte: There is no 24th element
    OUTPUTS
    '''
    pl=29*[0]
    pl[0]=A
    pl[1]=d
    pl[2]=T0
    pl[3]=eps
    pl[4]=kappa
    pl[5]=nn
    pl[6]=emin
    pl[7]=emax
    pl[8]=ebrk
    pl[9]=del1
    pl[10]del2
    pl[11]=nth
    pl[12]=nb
    pl[13]=B
    pl[14]=theta
    pl[15]=fs
    pl[16]=delf
    pl[17]=ed
    pl[18]=nf
    pl[19]=paf
    pl[20]=lca
    pl[21]=bda
    pl[22]=delmu
    pl[23]=a4
    pl[25]=fc
    pl[26]=fw
    pl[27]=mat
    pl[28]=q
    return pl

def get_emission(pl,n):
    '''
    Now, we create two-dimensional array which contains information 
    about all volume elements along the line-of-sight.
    
    INPUT
    pl: list of the parameters required for fitting
    n: Number of nodes along line-of-sight

    NOTE
    Vieweing angle varies from 80 to 110 degrees

    OUTPUT
    res: left and right polarized emission intensity in SFU
    '''
    parms=n*[None]
    for i in range(n):
        parms[i]=pl[:]
        parms[i][1]=parms[i][1]/n # But the length of an elementary interval should equal the source depth
        # divided by the number of intervals (all intervals are assumed to have the same length
        parms[i][14]=80.0+30.0*i/(n-1) # Viewing angles varies from 80 to 110 degrees
    res=MWTransfer.GET_MW(params)
    return res


def get_Tb(res,f):
    '''
    Get Brightness Temperatures from emission array
    INPUT
    res: Emission Array in SFU
    
    NOTE
    bmin=np.sqrt(ut.rad2sec(pl[0]/(7.25e7)**2))
    bmin=bmax

    OUTPUT
    Brightness temperature
    '''
    S=np.array(res[1])+np.array(res[2])
    bmin=np.sqrt(ut.rad2sec(pl[0]/(7.25e7)**2))
    bmax=np.sqrt(ut.rad2sec(pl[0]/(7.25e7)**2))
    Tb=ut.flux2Tb(S,bmin,bmax,f)
    return Tb



