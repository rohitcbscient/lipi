import numpy as np
import sys

# --- Constants ---
mdtor = np.pi / 180.0   # !dpi/180d0
ruser=1.05
thuser=90.0
phuser=0.0
coorduser='observer'
gridtype = 'PLANEOFSKY'

# Example: assume gridtype etc. are defined externally
# gridtype = 'PLANEOFSKY'
# xxmin, xxmax, yymin, yymax, ngrid, ngy, cmer, phio must be defined before

if gridtype.upper() == 'PLANEOFSKY':

    xxr = xxmax - xxmin
    yyr = yymax - yymin

    if ngy == ngrid:
        rrr = max([xxr, yyr])
        dx = rrr / float(ngrid)
        dy = dx
    else:
        dx = xxr / float(ngrid)
        dy = yyr / float(ngy)

    # nx=fix(round(xxr/dx))
    # ny=fix(round(yyr/dy))
    nx = int(np.fix(np.round(xxr / dx)))
    ny = int(np.fix(np.round(yyr / dy)))

    if nx != ngrid:
        print("NX changed", nx, ngrid)
    if ny != ngy:
        print("NY changed", ny, ngy)

    xxmaxin = xxmax
    yymaxin = yymax

    xxmax = xxmin + float(nx) * dx
    yymax = yymin + float(ny) * dy

    # make sure we are not trying to evaluate right at origin
    if abs(xxmin) < 1e-8: xxmin = 1e-8
    if abs(yymin) < 1e-8: yymin = 1e-8
    if abs(xxmax) < 1e-8: xxmax = 1e-8
    if abs(yymax) < 1e-8: yymax = 1e-8

    if abs(xxmaxin - xxmax) > 1e-7:
        print("XMAX changed", xxmaxin, xxmax)
    if abs(yymaxin - yymax) > 1e-7:
        print("YMAX changed", yymaxin, yymax)

    # dindgen(nx) => np.arange(nx, dtype=float)
    x1d = xxmin + np.arange(nx, dtype=float) * dx + dx / 2.0
    y1d = yymin + np.arange(ny, dtype=float) * dy + dy / 2.0

    xcenter = (xxmax + xxmin) / 2.0
    ycenter = (yymax + yymin) / 2.0

    # x = x1d # replicate(1.d0,ny)
    # y = replicate(1.d0,nx) # y1d
    x = np.outer(x1d, np.ones(ny))
    y = np.outer(np.ones(nx), y1d)

    r = np.sqrt(x**2 + y**2)
    th = np.arccos(y / r)

    west = np.where(x > 0.0)
    if west[0].size > 0:
        th[west] = 2.0 * np.pi - th[west]

    ph = x * 0.0 + (cmer - phio) * mdtor

    rheight = 0.0
    


# ; special case AZEQUI set
if azequi == 1:
    # ;
    # ; step in even elongation
    # ;
    elong = r * np.pi / 180.0       # !dtor
    r = distobs * np.tan(elong)
    # ;
    # ; the purpose of this is just to pick plane of sky intersections for the lines of sight
    # ; that are equally spaced in elongation rather than r_pos. The points in 3D space will be 
    # ; calculated the same as usual, dealing with non-parallel lines of sight in for_intlos
    # ;
    # ; alpha = th
    # ;  --- remember though this is a polar angle, 0 to 360 counterclockwise
    # ;
    # ; e_y = -elong*sin(alpha)
    # ; e_z = elong*cos(alpha)
    # ;
    # ; to be sorted out in for_plot.pro
    # ;
# ;

# ; put losoffset into same dimensions as r, th, ph
# ;  and get rid of 'NULL' for LOS integral/Physical Diagnostic (replace with 0.)
# ;
if isinstance(losoffset, (list, np.ndarray)) and np.isscalar(losoffset[0]) == False:
    losoffuse = np.array(losoffset, dtype=float) + r * 0.0
else:
    losoffuse = r * 0.0

# ;

# ; USERINPUT or USER grid type
if gridtype.upper() in ['USERINPUT', 'USER']:

    # ;
    # ; ruser, thuser, phuser (in DEGREES)  
    # ;  are spherical coords
    # ;  r should be radius of point
    # ;  th should be colatitude of point
    # ;  ph should be longitude of point 
    # ; colatitude/longitude are in two frames defined by COORDUSER
    # ; and described in detail in FOR_GRIDDEFAULTS
    # ;
    # ; First we will get arrays in the right dimension and change to radians
    # ;
    test_r = np.shape(ruser)
    test_th = np.shape(thuser)
    test_ph = np.shape(phuser)

    if len(test_r) != len(test_th) or len(test_r) != len(test_ph):
        print('ruser, thuser, phuser must have same dimensions')
        print('stopping in for_get_grid -- reset and try again.')
        raise SystemExit

    if len(test_r) == 0:
        r3D = np.zeros(1, dtype=float)
        th3D = np.zeros(1, dtype=float)
        ph3D = np.zeros(1, dtype=float)

        r3D[0] = ruser
        th3D[0] = thuser * mdtor
        ph3D[0] = phuser * mdtor

        if thuser > 180:
            print('thuser should be 0–180 degrees (colatitude)')
    else:
        r3D = np.array(ruser, dtype=float)
        th3D = np.array(thuser, dtype=float) * mdtor
        ph3D = np.array(phuser, dtype=float) * mdtor

        test = np.where(thuser > 180)
        if test[0].size > 0:
            print('thuser should be 0–180 degrees (colatitude)')


