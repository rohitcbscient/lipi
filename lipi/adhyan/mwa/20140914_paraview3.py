# trace generated using paraview version 5.7.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
import glob
import sys
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()
listvtk=sorted(glob.glob('/sdata/20140914_hmi/hmi/hmi*.vtk'))[0:150][::10]
for i in range(11,len(listvtk)):
    print listvtk[i]
    # create a new 'Legacy VTK Reader'
    hmi_20140914_015530_magnetogramfitsrotsavcubesav_0001vtk = LegacyVTKReader(FileNames=[listvtk[i]])
    # create a new 'Stream Tracer'
    streamTracer1 = StreamTracer(Input=hmi_20140914_015530_magnetogramfitsrotsavcubesav_0001vtk,SeedType='Point Source')
    streamTracer1.Vectors = ['POINTS', 'Magnetic Field']
    streamTracer1.MaximumStreamlineLength = 1664.0
    # init the 'High Resolution Line Source' selected for 'SeedType'
    streamTracer1.SeedType.Center = [177.0, 223.0, 20]
    streamTracer1.SeedType.Radius = 50.0
    # Create a new 'SpreadSheet View'
    spreadSheetView1 = CreateView('SpreadSheetView')
    spreadSheetView1.ColumnToSort = ''
    spreadSheetView1.BlockSize = 150
    # uncomment following to set a specific view size
    # spreadSheetView1.ViewSize = [400, 400]
    # show data in view
    streamTracer1Display_1 = Show(streamTracer1, spreadSheetView1)
    # export view
    ExportView(listvtk[i]+'.csv', view=spreadSheetView1)

