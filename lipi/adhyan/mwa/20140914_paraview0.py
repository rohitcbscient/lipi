# trace generated using paraview version 5.7.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'Legacy VTK Reader'
hmi_20140914_015445_magnetogramfitsrotsavcubesav_0000vtk = LegacyVTKReader(FileNames=['/sdata/20140914_hmi/hmi/hmi_20140914_015445_magnetogram.fitsrot.sav.cube.sav_0000.vtk'])

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1817, 865]

# show data in view
hmi_20140914_015445_magnetogramfitsrotsavcubesav_0000vtkDisplay = Show(hmi_20140914_015445_magnetogramfitsrotsavcubesav_0000vtk, renderView1)

# trace defaults for the display properties.
hmi_20140914_015445_magnetogramfitsrotsavcubesav_0000vtkDisplay.Representation = 'Outline'
hmi_20140914_015445_magnetogramfitsrotsavcubesav_0000vtkDisplay.ColorArrayName = ['POINTS', '']
hmi_20140914_015445_magnetogramfitsrotsavcubesav_0000vtkDisplay.OSPRayScaleArray = 'Magnitude'
hmi_20140914_015445_magnetogramfitsrotsavcubesav_0000vtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
hmi_20140914_015445_magnetogramfitsrotsavcubesav_0000vtkDisplay.SelectOrientationVectors = 'Magnetic Field'
hmi_20140914_015445_magnetogramfitsrotsavcubesav_0000vtkDisplay.ScaleFactor = 33.4
hmi_20140914_015445_magnetogramfitsrotsavcubesav_0000vtkDisplay.SelectScaleArray = 'Magnitude'
hmi_20140914_015445_magnetogramfitsrotsavcubesav_0000vtkDisplay.GlyphType = 'Arrow'
hmi_20140914_015445_magnetogramfitsrotsavcubesav_0000vtkDisplay.GlyphTableIndexArray = 'Magnitude'
hmi_20140914_015445_magnetogramfitsrotsavcubesav_0000vtkDisplay.GaussianRadius = 1.67
hmi_20140914_015445_magnetogramfitsrotsavcubesav_0000vtkDisplay.SetScaleArray = ['POINTS', 'Magnitude']
hmi_20140914_015445_magnetogramfitsrotsavcubesav_0000vtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
hmi_20140914_015445_magnetogramfitsrotsavcubesav_0000vtkDisplay.OpacityArray = ['POINTS', 'Magnitude']
hmi_20140914_015445_magnetogramfitsrotsavcubesav_0000vtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'
hmi_20140914_015445_magnetogramfitsrotsavcubesav_0000vtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
hmi_20140914_015445_magnetogramfitsrotsavcubesav_0000vtkDisplay.PolarAxes = 'PolarAxesRepresentation'
hmi_20140914_015445_magnetogramfitsrotsavcubesav_0000vtkDisplay.ScalarOpacityUnitDistance = 1.7320542659082296

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
hmi_20140914_015445_magnetogramfitsrotsavcubesav_0000vtkDisplay.ScaleTransferFunction.Points = [0.06843648464, 0.0, 0.5, 0.0, 1226.9029256, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
hmi_20140914_015445_magnetogramfitsrotsavcubesav_0000vtkDisplay.OpacityTransferFunction.Points = [0.06843648464, 0.0, 0.5, 0.0, 1226.9029256, 1.0, 0.5, 0.0]

# reset view to fit data
renderView1.ResetCamera()

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Slice'
slice1 = Slice(Input=hmi_20140914_015445_magnetogramfitsrotsavcubesav_0000vtk)
slice1.SliceType = 'Plane'
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [167.0, 166.5, 167.0]

# Properties modified on slice1.SliceType
slice1.SliceType.Origin = [167.0, 166.5, 0.5]
slice1.SliceType.Normal = [0.0, 0.0, 1.0]

# show data in view
slice1Display = Show(slice1, renderView1)

# get color transfer function/color map for 'Magnitude'
magnitudeLUT = GetColorTransferFunction('Magnitude')

# trace defaults for the display properties.
slice1Display.Representation = 'Surface'
slice1Display.ColorArrayName = ['POINTS', 'Magnitude']
slice1Display.LookupTable = magnitudeLUT
slice1Display.OSPRayScaleArray = 'Magnitude'
slice1Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice1Display.SelectOrientationVectors = 'Magnetic Field'
slice1Display.ScaleFactor = 33.4
slice1Display.SelectScaleArray = 'Magnitude'
slice1Display.GlyphType = 'Arrow'
slice1Display.GlyphTableIndexArray = 'Magnitude'
slice1Display.GaussianRadius = 1.67
slice1Display.SetScaleArray = ['POINTS', 'Magnitude']
slice1Display.ScaleTransferFunction = 'PiecewiseFunction'
slice1Display.OpacityArray = ['POINTS', 'Magnitude']
slice1Display.OpacityTransferFunction = 'PiecewiseFunction'
slice1Display.DataAxesGrid = 'GridAxesRepresentation'
slice1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
slice1Display.ScaleTransferFunction.Points = [2.2847777665, 0.0, 0.5, 0.0, 1129.04292455, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
slice1Display.OpacityTransferFunction.Points = [2.2847777665, 0.0, 0.5, 0.0, 1129.04292455, 1.0, 0.5, 0.0]

# show color bar/color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# get opacity transfer function/opacity map for 'Magnitude'
magnitudePWF = GetOpacityTransferFunction('Magnitude')

# create a new 'Stream Tracer'
streamTracer1 = StreamTracer(Input=slice1,
    SeedType='High Resolution Line Source')
streamTracer1.Vectors = ['POINTS', 'Magnetic Field']
streamTracer1.MaximumStreamlineLength = 334.0

# init the 'High Resolution Line Source' selected for 'SeedType'
streamTracer1.SeedType.Point2 = [334.0, 333.0, 0.0]

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=streamTracer1.SeedType)

# Properties modified on streamTracer1
streamTracer1.MaximumStreamlineLength = 1336.0
streamTracer1.SeedType = 'Point Source'

# show data in view
streamTracer1Display = Show(streamTracer1, renderView1)

# trace defaults for the display properties.
streamTracer1Display.Representation = 'Surface'
streamTracer1Display.ColorArrayName = ['POINTS', 'Magnitude']
streamTracer1Display.LookupTable = magnitudeLUT
streamTracer1Display.OSPRayScaleArray = 'Magnitude'
streamTracer1Display.OSPRayScaleFunction = 'PiecewiseFunction'
streamTracer1Display.SelectOrientationVectors = 'Normals'
streamTracer1Display.ScaleFactor = 0.11270481944084168
streamTracer1Display.SelectScaleArray = 'Magnitude'
streamTracer1Display.GlyphType = 'Arrow'
streamTracer1Display.GlyphTableIndexArray = 'Magnitude'
streamTracer1Display.GaussianRadius = 0.005635240972042084
streamTracer1Display.SetScaleArray = ['POINTS', 'Magnitude']
streamTracer1Display.ScaleTransferFunction = 'PiecewiseFunction'
streamTracer1Display.OpacityArray = ['POINTS', 'Magnitude']
streamTracer1Display.OpacityTransferFunction = 'PiecewiseFunction'
streamTracer1Display.DataAxesGrid = 'GridAxesRepresentation'
streamTracer1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
streamTracer1Display.ScaleTransferFunction.Points = [659.4791934384598, 0.0, 0.5, 0.0, 694.3051613635297, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
streamTracer1Display.OpacityTransferFunction.Points = [659.4791934384598, 0.0, 0.5, 0.0, 694.3051613635297, 1.0, 0.5, 0.0]

# hide data in view
Hide(slice1, renderView1)

# show color bar/color legend
streamTracer1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(slice1)

# show data in view
slice1Display = Show(slice1, renderView1)

# show color bar/color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

# set active source
SetActiveSource(streamTracer1)

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [372.08057390430804, -685.194220548648, 1195.383344461329]
renderView1.CameraFocalPoint = [167.00000000000006, 166.50000000000006, 166.99999999999991]
renderView1.CameraViewUp = [-0.24624489096980598, 0.7217598363757127, 0.6468587112083851]
renderView1.CameraParallelScale = 288.96409811601166

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).