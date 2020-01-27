
import mayavi.mlab as mlab
import matplotlib.pyplot as plt
import numpy as np
from mayavi.modules.scalar_cut_plane import ScalarCutPlane

plt.style.use('~/vla_scripts/plt_style.py')

def xyz_line_plot(T,freq,xc,yc,xl,xr,yl,yr,t):
    z,x,y = np.mgrid[freq[0]:freq[-1]:96j,xl:xr:24j,yl:yr:25j]
    mlab.figure(1,size = (1824,1224), bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5))
    iso = mlab.contour3d(z,x,y,T,vmin=T.min(),vmax=T.max(),opacity=0.1, colormap='spring')
    iso.contour.number_of_contours = 15
    mlab.view(220,240,distance=400)
    #mlab.orientation_axes()
    mlab.outline(extent=[freq[0],freq[-1],xl,xr,yl,yr])
    mlab.axes(color=(0.5,0.5,0.5), line_width=1, xlabel='FREQUENCY (MHz)', ylabel='X', zlabel='Y',nb_labels=5)
    cl=mlab.colorbar(title='Temperature (MK)', orientation='vertical')
    cl.label_text_property.font_family = 'courier'
    cl.label_text_property.font_size = 4
    cl.data_range=(0.5,9.5)
    mlab.title(t,size=0.5,line_width=1.0)
    # Create the points
    zz=np.linspace(freq[0],freq[-1],96)
    src = mlab.plot3d(freq,xc,yc,zz,opacity=1.0,line_width=200,colormap='copper',tube_radius=2.)
    # Connect them
    #src.mlab_source.dataset.lines = connections
    #mlab.savefig(filename='3d_plots_1s/20120225_vla_'+"%03d" % l+'_2.png', figure=mlab.gcf(), magnification=1)
    scp = ScalarCutPlane() # set scp as ScalarCutPlane() module
    iso.add_module(scp) # add module to the scene
    scp.implicit_plane.normal = (1, 0, 0) # set normal to Ox axis
    # set origin to (i=10, j=25, k=25) i.e. integers for a structured grid
    scp.implicit_plane.origin = (10, 25, 25)
    # set origin to (x=1.0, y=2.5, z=2.5) i.e. reals for unstructured grids
    # scp.implicit_plane.origin = (1.0, 2.5, 2.5)
    scp.implicit_plane.widget.enabled = False
    scp.actor.property.diffuse = 0.0 # set some color properties
    scp.actor.property.ambient = 1.0 # 
    scp.actor.property.opacity = 1.0 #
    scp.module_manager.scalar_lut_manager.data_range = [0, 1]
    mlab.show()
    mlab.close()

def xyz_plot(T,freq,xl,xr,yl,yr,t):
    z,x,y = np.mgrid[freq[0]:freq[-1]:96j,xl:xr:24j,yl:yr:25j]
    mlab.figure(1,size = (1824,1224), bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5))
    iso = mlab.contour3d(z,x,y,T,vmin=T.min(),vmax=T.max(),opacity=0.3, colormap='gist_rainbow')
    iso.contour.number_of_contours = 15
    mlab.view(220,240,distance=600)
    #mlab.orientation_axes()
    mlab.outline(extent=[freq[0],freq[-1],xl,xr,yl,yr])
    mlab.axes(color=(0.5,0.5,0.5), line_width=1, xlabel='FREQUENCY (MHz)', ylabel='X', zlabel='Y',nb_labels=5)
    cl=mlab.colorbar(title='Temperature (MK)', orientation='vertical')
    cl.label_text_property.font_family = 'courier'
    cl.label_text_property.font_size = 4
    cl.data_range=(0.5,6.5)
    mlab.title(t,size=0.5,line_width=1.0)
    #mlab.savefig(filename='3d_plots_1s/20120225_vla_'+"%03d" % l+'_2.png', figure=mlab.gcf(), magnification=1)
    mlab.show()
    mlab.close()
