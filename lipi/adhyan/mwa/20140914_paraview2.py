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

# find source

listvtk=sorted(glob.glob('/sdata/20140914_hmi/hmi/hmi*.vtk'))
#for i in range(len(listvtk)):
for i in range(2):
    print(listvtk[i])
    hmi_20140914_vtk = FindSource(listvtk[i])
    # create a new 'Stream Tracer'
    streamTracer1 = StreamTracer(Input=hmi_20140914_vtk,SeedType='High Resolution Line Source')
    streamTracer1.Vectors = ['POINTS', 'Magnetic Field']
    streamTracer1.MaximumStreamlineLength = 416.0
    # init the 'High Resolution Line Source' selected for 'SeedType'
    streamTracer1.SeedType.Point2 = [416.0, 416.0, 249.0]
    # toggle 3D widget visibility (only when running from the GUI)
    Show3DWidgets(proxy=streamTracer1.SeedType)
    # Properties modified on streamTracer1
    streamTracer1.MaximumStreamlineLength = 1664.0
    streamTracer1.SeedType = 'Point Source'
    # Rescale transfer function
    magnitudeLUT.RescaleTransferFunction(0.9931426738473721, 1038.9114228660312)
    # Rescale transfer function
    magnitudePWF.RescaleTransferFunction(0.9931426738473721, 1038.9114228660312)
    # get layout
    layout1 = GetLayout()
    # split cell
    layout1.SplitHorizontal(0, 0.5)
    # set active view
    SetActiveView(None)
    # Create a new 'SpreadSheet View'
    spreadSheetView1 = CreateView('SpreadSheetView')
    spreadSheetView1.ColumnToSort = ''
    spreadSheetView1.BlockSize = 1024
    # uncomment following to set a specific view size
    # spreadSheetView1.ViewSize = [400, 400]
    # show data in view
    streamTracer1Display_1 = Show(streamTracer1, spreadSheetView1)
    # export view
    ExportView(listvtk[i]+'.csv', view=spreadSheetView1)

sys.exit()
# set active view
SetActiveView(renderView1)

# Properties modified on magnitudeLUT
magnitudeLUT.RGBPoints = [0.9931426738473721, 0.231373, 0.298039, 0.752941, 516.8260498046875, 0.865003, 0.865003, 0.865003, 1038.9114228660312, 0.705882, 0.0156863, 0.14902]

# Properties modified on magnitudeLUT
magnitudeLUT.RGBPoints = [0.9931426738473721, 0.231373, 0.298039, 0.752941, 507.4472351074219, 0.865003, 0.865003, 0.865003, 1038.9114228660312, 0.705882, 0.0156863, 0.14902]

# Properties modified on magnitudeLUT
magnitudeLUT.RGBPoints = [0.9931426738473721, 0.231373, 0.298039, 0.752941, 501.1947021484375, 0.865003, 0.865003, 0.865003, 1038.9114228660312, 0.705882, 0.0156863, 0.14902]

# Properties modified on magnitudeLUT
magnitudeLUT.RGBPoints = [0.9931426738473721, 0.231373, 0.298039, 0.752941, 491.81591796875, 0.865003, 0.865003, 0.865003, 1038.9114228660312, 0.705882, 0.0156863, 0.14902]

# Properties modified on magnitudeLUT
magnitudeLUT.RGBPoints = [0.9931426738473721, 0.231373, 0.298039, 0.752941, 463.6795959472656, 0.865003, 0.865003, 0.865003, 1038.9114228660312, 0.705882, 0.0156863, 0.14902]

# Properties modified on magnitudeLUT
magnitudeLUT.RGBPoints = [0.9931426738473721, 0.231373, 0.298039, 0.752941, 460.5533142089844, 0.865003, 0.865003, 0.865003, 1038.9114228660312, 0.705882, 0.0156863, 0.14902]

# Properties modified on magnitudeLUT
magnitudeLUT.RGBPoints = [0.9931426738473721, 0.231373, 0.298039, 0.752941, 438.66949462890625, 0.865003, 0.865003, 0.865003, 1038.9114228660312, 0.705882, 0.0156863, 0.14902]

# Properties modified on magnitudeLUT
magnitudeLUT.RGBPoints = [0.9931426738473721, 0.231373, 0.298039, 0.752941, 435.5432434082031, 0.865003, 0.865003, 0.865003, 1038.9114228660312, 0.705882, 0.0156863, 0.14902]

# Properties modified on magnitudeLUT
magnitudeLUT.RGBPoints = [0.9931426738473721, 0.231373, 0.298039, 0.752941, 429.2907409667969, 0.865003, 0.865003, 0.865003, 1038.9114228660312, 0.705882, 0.0156863, 0.14902]

# Properties modified on magnitudeLUT
magnitudeLUT.RGBPoints = [0.9931426738473721, 0.231373, 0.298039, 0.752941, 423.0382080078125, 0.865003, 0.865003, 0.865003, 1038.9114228660312, 0.705882, 0.0156863, 0.14902]

# Properties modified on magnitudeLUT
magnitudeLUT.RGBPoints = [0.9931426738473721, 0.231373, 0.298039, 0.752941, 413.659423828125, 0.865003, 0.865003, 0.865003, 1038.9114228660312, 0.705882, 0.0156863, 0.14902]

# Properties modified on magnitudeLUT
magnitudeLUT.RGBPoints = [0.9931426738473721, 0.231373, 0.298039, 0.752941, 398.02813720703125, 0.865003, 0.865003, 0.865003, 1038.9114228660312, 0.705882, 0.0156863, 0.14902]

# Properties modified on magnitudeLUT
magnitudeLUT.RGBPoints = [0.9931426738473721, 0.231373, 0.298039, 0.752941, 391.7756042480469, 0.865003, 0.865003, 0.865003, 1038.9114228660312, 0.705882, 0.0156863, 0.14902]

# Properties modified on magnitudeLUT
magnitudeLUT.RGBPoints = [0.9931426738473721, 0.231373, 0.298039, 0.752941, 388.64935302734375, 0.865003, 0.865003, 0.865003, 1038.9114228660312, 0.705882, 0.0156863, 0.14902]

# Properties modified on magnitudeLUT
magnitudeLUT.RGBPoints = [0.9931426738473721, 0.231373, 0.298039, 0.752941, 382.3968200683594, 0.865003, 0.865003, 0.865003, 1038.9114228660312, 0.705882, 0.0156863, 0.14902]

# Properties modified on magnitudeLUT
magnitudeLUT.RGBPoints = [0.9931426738473721, 0.231373, 0.298039, 0.752941, 379.27056884765625, 0.865003, 0.865003, 0.865003, 1038.9114228660312, 0.705882, 0.0156863, 0.14902]

# Properties modified on magnitudeLUT
magnitudeLUT.RGBPoints = [0.9931426738473721, 0.231373, 0.298039, 0.752941, 373.0180358886719, 0.865003, 0.865003, 0.865003, 1038.9114228660312, 0.705882, 0.0156863, 0.14902]

# Properties modified on magnitudeLUT
magnitudeLUT.RGBPoints = [0.9931426738473721, 0.231373, 0.298039, 0.752941, 369.89178466796875, 0.865003, 0.865003, 0.865003, 1038.9114228660312, 0.705882, 0.0156863, 0.14902]

# Properties modified on magnitudeLUT
magnitudeLUT.RGBPoints = [0.9931426738473721, 0.231373, 0.298039, 0.752941, 357.3867492675781, 0.865003, 0.865003, 0.865003, 1038.9114228660312, 0.705882, 0.0156863, 0.14902]

# Properties modified on magnitudeLUT
magnitudeLUT.RGBPoints = [0.9931426738473721, 0.231373, 0.298039, 0.752941, 338.6291809082031, 0.865003, 0.865003, 0.865003, 1038.9114228660312, 0.705882, 0.0156863, 0.14902]

# Properties modified on magnitudeLUT
magnitudeLUT.RGBPoints = [0.9931426738473721, 0.231373, 0.298039, 0.752941, 326.1241455078125, 0.865003, 0.865003, 0.865003, 1038.9114228660312, 0.705882, 0.0156863, 0.14902]

# Properties modified on magnitudeLUT
magnitudeLUT.RGBPoints = [0.9931426738473721, 0.231373, 0.298039, 0.752941, 322.9978942871094, 0.865003, 0.865003, 0.865003, 1038.9114228660312, 0.705882, 0.0156863, 0.14902]

# Properties modified on magnitudeLUT
magnitudeLUT.RGBPoints = [0.9931426738473721, 0.231373, 0.298039, 0.752941, 291.73529052734375, 0.865003, 0.865003, 0.865003, 1038.9114228660312, 0.705882, 0.0156863, 0.14902]

# Properties modified on magnitudeLUT
magnitudeLUT.RGBPoints = [0.9931426738473721, 0.231373, 0.298039, 0.752941, 279.2302551269531, 0.865003, 0.865003, 0.865003, 1038.9114228660312, 0.705882, 0.0156863, 0.14902]

# Properties modified on magnitudeLUT
magnitudeLUT.RGBPoints = [0.9931426738473721, 0.231373, 0.298039, 0.752941, 166.68492126464844, 0.865003, 0.865003, 0.865003, 1038.9114228660312, 0.705882, 0.0156863, 0.14902]

# Properties modified on magnitudeLUT
magnitudeLUT.RGBPoints = [0.9931426738473721, 0.231373, 0.298039, 0.752941, 147.92735290527344, 0.865003, 0.865003, 0.865003, 1038.9114228660312, 0.705882, 0.0156863, 0.14902]

# Properties modified on magnitudeLUT
magnitudeLUT.RGBPoints = [0.9931426738473721, 0.231373, 0.298039, 0.752941, 141.67483520507812, 0.865003, 0.865003, 0.865003, 1038.9114228660312, 0.705882, 0.0156863, 0.14902]

# Properties modified on magnitudeLUT
magnitudeLUT.RGBPoints = [0.9931426738473721, 0.231373, 0.298039, 0.752941, 138.548583984375, 0.865003, 0.865003, 0.865003, 1038.9114228660312, 0.705882, 0.0156863, 0.14902]

# Properties modified on magnitudeLUT
magnitudeLUT.RGBPoints = [0.9931426738473721, 0.231373, 0.298039, 0.752941, 135.4223175048828, 0.865003, 0.865003, 0.865003, 1038.9114228660312, 0.705882, 0.0156863, 0.14902]

# Properties modified on magnitudeLUT
magnitudeLUT.RGBPoints = [0.9931426738473721, 0.231373, 0.298039, 0.752941, 132.2960662841797, 0.865003, 0.865003, 0.865003, 1038.9114228660312, 0.705882, 0.0156863, 0.14902]

# Properties modified on magnitudeLUT
magnitudeLUT.RGBPoints = [0.9931426738473721, 0.231373, 0.298039, 0.752941, 129.1697998046875, 0.865003, 0.865003, 0.865003, 1038.9114228660312, 0.705882, 0.0156863, 0.14902]

# Properties modified on magnitudeLUT
magnitudeLUT.RGBPoints = [0.9931426738473721, 0.231373, 0.298039, 0.752941, 104.15971374511719, 0.865003, 0.865003, 0.865003, 1038.9114228660312, 0.705882, 0.0156863, 0.14902]

# Properties modified on magnitudeLUT
magnitudeLUT.RGBPoints = [0.9931426738473721, 0.231373, 0.298039, 0.752941, 88.5284194946289, 0.865003, 0.865003, 0.865003, 1038.9114228660312, 0.705882, 0.0156863, 0.14902]

# Properties modified on magnitudeLUT
magnitudeLUT.RGBPoints = [0.9931426738473721, 0.231373, 0.298039, 0.752941, 79.1496353149414, 0.865003, 0.865003, 0.865003, 1038.9114228660312, 0.705882, 0.0156863, 0.14902]

# Properties modified on magnitudeLUT
magnitudeLUT.RGBPoints = [0.9931426738473721, 0.231373, 0.298039, 0.752941, 76.02337646484375, 0.865003, 0.865003, 0.865003, 1038.9114228660312, 0.705882, 0.0156863, 0.14902]

# Properties modified on magnitudeLUT
magnitudeLUT.RGBPoints = [0.9931426738473721, 0.231373, 0.298039, 0.752941, 54.139556884765625, 0.865003, 0.865003, 0.865003, 1038.9114228660312, 0.705882, 0.0156863, 0.14902]

# Properties modified on magnitudeLUT
magnitudeLUT.RGBPoints = [0.9931426738473721, 0.231373, 0.298039, 0.752941, 44.760780334472656, 0.865003, 0.865003, 0.865003, 1038.9114228660312, 0.705882, 0.0156863, 0.14902]

# Properties modified on magnitudeLUT
magnitudeLUT.RGBPoints = [0.9931426738473721, 0.231373, 0.298039, 0.752941, 47.88703918457031, 0.865003, 0.865003, 0.865003, 1038.9114228660312, 0.705882, 0.0156863, 0.14902]

# Properties modified on magnitudeLUT
magnitudeLUT.RGBPoints = [0.9931426738473721, 0.231373, 0.298039, 0.752941, 54.139556884765625, 0.865003, 0.865003, 0.865003, 1038.9114228660312, 0.705882, 0.0156863, 0.14902]

# Properties modified on magnitudeLUT
magnitudeLUT.RGBPoints = [0.9931426738473721, 0.231373, 0.298039, 0.752941, 60.3920783996582, 0.865003, 0.865003, 0.865003, 1038.9114228660312, 0.705882, 0.0156863, 0.14902]

# Properties modified on magnitudeLUT
magnitudeLUT.RGBPoints = [0.9931426738473721, 0.231373, 0.298039, 0.752941, 82.27589416503906, 0.865003, 0.865003, 0.865003, 1038.9114228660312, 0.705882, 0.0156863, 0.14902]

# Properties modified on magnitudeLUT
magnitudeLUT.RGBPoints = [0.9931426738473721, 0.231373, 0.298039, 0.752941, 85.40215301513672, 0.865003, 0.865003, 0.865003, 1038.9114228660312, 0.705882, 0.0156863, 0.14902]

# Properties modified on magnitudeLUT
magnitudeLUT.RGBPoints = [0.9931426738473721, 0.231373, 0.298039, 0.752941, 88.5284194946289, 0.865003, 0.865003, 0.865003, 1038.9114228660312, 0.705882, 0.0156863, 0.14902]

# Properties modified on magnitudeLUT
magnitudeLUT.RGBPoints = [0.9931426738473721, 0.231373, 0.298039, 0.752941, 94.78093719482422, 0.865003, 0.865003, 0.865003, 1038.9114228660312, 0.705882, 0.0156863, 0.14902]

# Properties modified on magnitudeLUT
magnitudeLUT.RGBPoints = [0.9931426738473721, 0.231373, 0.298039, 0.752941, 97.90719604492188, 0.865003, 0.865003, 0.865003, 1038.9114228660312, 0.705882, 0.0156863, 0.14902]

# Properties modified on magnitudeLUT
magnitudeLUT.RGBPoints = [0.9931426738473721, 0.231373, 0.298039, 0.752941, 94.78093719482422, 0.865003, 0.865003, 0.865003, 1038.9114228660312, 0.705882, 0.0156863, 0.14902]

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=streamTracer1.SeedType)

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [-211.483082122164, -698.8938192824869, 1234.2504164475292]
renderView1.CameraFocalPoint = [207.99999999999977, 208.0000000000002, 124.49999999999959]
renderView1.CameraViewUp = [-0.1241764108644562, 0.7908180075438999, 0.5993221987621669]
renderView1.CameraParallelScale = 319.4186124821157

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
