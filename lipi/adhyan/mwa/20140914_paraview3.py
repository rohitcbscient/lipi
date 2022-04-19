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
listvtk=sorted(glob.glob('/sdata/20140914_hmi/2014-09-14/hmi.*.vtk'))#[0:150][::10]
for i in range(len(listvtk)):
    print listvtk[i]
    # create a new 'Legacy VTK Reader'
    hmi_20140914_015530_magnetogramfitsrotsavcubesav_0001vtk = LegacyVTKReader(FileNames=[listvtk[i]])
    # create a new 'Stream Tracer'
    streamTracer1 = StreamTracer(Input=hmi_20140914_015530_magnetogramfitsrotsavcubesav_0001vtk,SeedType='Point Source')
    streamTracer1.Vectors = ['POINTS', 'Magnetic Field']
    streamTracer1.MaximumStreamlineLength = 1500
    # init the 'High Resolution Line Source' selected for 'SeedType'
    #streamTracer1.SeedType.Center = [85.0, 183.0, 1.5]
    #streamTracer1.SeedType.Center = [80.0, 140.0, 1.5]
    #streamTracer1.SeedType.Center = [183.0, 123.0, 180]
    #streamTracer1.SeedType.Center = [149.0, 105.0, 40]
    streamTracer1.SeedType.Center = [114.0, 168.0, 5] # X-ray
    #streamTracer1.SeedType.Center = [208.0, 117.0,110] # Large open loop
    streamTracer1.SeedType.Radius = 30.0
    # Create a new 'SpreadSheet View'
    spreadSheetView1 = CreateView('SpreadSheetView')
    spreadSheetView1.ColumnToSort = ''
    spreadSheetView1.BlockSize = 250
    # uncomment following to set a specific view size
    # spreadSheetView1.ViewSize = [400, 400]
    # show data in view
    streamTracer1Display_1 = Show(streamTracer1, spreadSheetView1)
    # export view
    ExportView(listvtk[i]+'.x-ray.csv', view=spreadSheetView1)

