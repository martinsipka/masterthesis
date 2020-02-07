#### import the simple module from the paraview
from paraview.simple import *
import os
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

cwd = os.getcwd()
outputDir = "/home/martin/praha/DiplomovaPraca/Tex/Praca/img/Snapshots/" + cwd.split('/')[-2] + "/"
inputDir = cwd + "/"

# create a new 'XDMF Reader'
btens = XDMFReader(FileNames=[inputDir + 'Btens0.xdmf', inputDir + 'Btens10.xdmf', inputDir + 'Btens20.xdmf', inputDir + 'Btens30.xdmf', inputDir + 'Btens40.xdmf', inputDir + 'Btens50.xdmf', inputDir + 'Btens60.xdmf', inputDir + 'Btens70.xdmf', inputDir + 'Btens80.xdmf', inputDir + 'Btens90.xdmf', inputDir + 'Btens100.xdmf', inputDir + 'Btens110.xdmf', inputDir + 'Btens120.xdmf', inputDir + 'Btens130.xdmf', inputDir + 'Btens140.xdmf', inputDir + 'Btens150.xdmf',inputDir + 'Btens160.xdmf', inputDir + 'Btens170.xdmf', inputDir + 'Btens180.xdmf', inputDir + 'Btens190.xdmf', inputDir + 'Btens200.xdmf'])

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# create a new 'XDMF Reader'
nvect = XDMFReader(FileNames=[inputDir + 'nvect0.xdmf', inputDir + 'nvect10.xdmf', inputDir + 'nvect20.xdmf', inputDir + 'nvect30.xdmf', inputDir + 'nvect40.xdmf', inputDir + 'nvect50.xdmf', inputDir + 'nvect60.xdmf', inputDir + 'nvect70.xdmf', inputDir + 'nvect80.xdmf', inputDir + 'nvect90.xdmf', inputDir + 'nvect100.xdmf', inputDir + 'nvect110.xdmf', inputDir + 'nvect120.xdmf', inputDir + 'nvect130.xdmf', inputDir + 'nvect140.xdmf', inputDir + 'nvect150.xdmf', inputDir + 'nvect160.xdmf', inputDir + 'nvect170.xdmf', inputDir + 'nvect180.xdmf', inputDir + 'nvect190.xdmf', inputDir + 'nvect200.xdmf'])

# create a new 'XDMF Reader'
velo = XDMFReader(FileNames=[inputDir + 'velo0.xdmf', inputDir + 'velo10.xdmf', inputDir + 'velo20.xdmf', inputDir + 'velo30.xdmf', inputDir + 'velo40.xdmf', inputDir + 'velo50.xdmf', inputDir + 'velo60.xdmf', inputDir + 'velo70.xdmf', inputDir + 'velo80.xdmf', inputDir + 'velo90.xdmf', inputDir + 'velo100.xdmf', inputDir + 'velo110.xdmf',  inputDir + 'velo120.xdmf',  inputDir + 'velo130.xdmf',  inputDir + 'velo140.xdmf',  inputDir + 'velo150.xdmf',  inputDir + 'velo160.xdmf',  inputDir + 'velo170.xdmf',  inputDir + 'velo180.xdmf',  inputDir + 'velo190.xdmf',  inputDir + 'velo200.xdmf'])

# Properties modified on velo
velo.GridStatus = ['mesh']

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1481, 786]

# show data in view
veloDisplay = Show(velo, renderView1)
# trace defaults for the display properties.
veloDisplay.Representation = 'Surface'

# reset view to fit data
renderView1.ResetCamera()

#changing interaction mode based on data extents
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [1.0, 0.5, 10000.0]
renderView1.CameraFocalPoint = [1.0, 0.5, 0.0]

# Properties modified on nvect
nvect.GridStatus = ['mesh']

# show data in view
nvectDisplay = Show(nvect, renderView1)
# trace defaults for the display properties.
nvectDisplay.Representation = 'Surface'

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(velo)

# create a new 'Glyph'
glyph1 = Glyph(Input=velo,
    GlyphType='Arrow')

# Properties modified on glyph1
glyph1.ScaleMode = 'vector'

# show data in view
glyph1Display = Show(glyph1, renderView1)
# trace defaults for the display properties.
glyph1Display.Representation = 'Surface'

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(glyph1Display, ('POINTS', 'GlyphVector', 'Magnitude'))

# rescale color and/or opacity maps used to include current data range
glyph1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'GlyphVector'
glyphVectorLUT = GetColorTransferFunction('GlyphVector')

# set active source
SetActiveSource(velo)

# hide data in view
Hide(glyph1, renderView1)

# show data in view
veloDisplay = Show(velo, renderView1)

# destroy glyph1
Delete(glyph1)
del glyph1

# set active source
SetActiveSource(nvect)

# set active source
SetActiveSource(velo)

# set scalar coloring
ColorBy(veloDisplay, ('POINTS', 'v', 'Magnitude'))

# rescale color and/or opacity maps used to include current data range
veloDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
veloDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'v'
vLUT = GetColorTransferFunction('v')

# hide color bar/color legend
veloDisplay.SetScalarBarVisibility(renderView1, False)

# set active source
SetActiveSource(nvect)

# create a new 'Glyph'
glyph1 = Glyph(Input=nvect,
    GlyphType='Arrow')

# Properties modified on glyph1
glyph1.ScaleMode = 'vector'

# show data in view
glyph1Display = Show(glyph1, renderView1)
# trace defaults for the display properties.
glyph1Display.Representation = 'Surface'

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(glyph1Display, ('POINTS', 'GlyphVector', 'Magnitude'))

# rescale color and/or opacity maps used to include current data range
glyph1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, True)

# Rescale transfer function
glyphVectorLUT.RescaleTransferFunction(0.764306282017, 1.2)

# get opacity transfer function/opacity map for 'GlyphVector'
glyphVectorPWF = GetOpacityTransferFunction('GlyphVector')

# Rescale transfer function
glyphVectorPWF.RescaleTransferFunction(0.764306282017, 1.2)

# set active source
SetActiveSource(velo)

# hide data in view
Hide(nvect, renderView1)

# Rescale transfer function
vLUT.RescaleTransferFunction(0.0, 2.0)

# get opacity transfer function/opacity map for 'v'
vPWF = GetOpacityTransferFunction('v')

# Rescale transfer function
vPWF.RescaleTransferFunction(0.0, 2.0)

# set active source
SetActiveSource(glyph1)

# hide color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, False)

# Properties modified on renderView1
renderView1.OrientationAxesVisibility = 0

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [1.0, 0.5, 10000.0]
renderView1.CameraFocalPoint = [1.0, 0.5, 0.0]
renderView1.CameraParallelScale = 0.7636322578716578

animationScene1.GoToNext()

animationScene1.GoToPrevious()

# save screenshot
SaveScreenshot(outputDir + 'time0s.png', renderView1, ImageResolution=[1481, 786],
    FontScaling='Scale fonts proportionally',
    OverrideColorPalette='WhiteBackground',
    StereoMode='No change',
    TransparentBackground=1,
    ImageQuality=100)

animationScene1.GoToNext()

animationScene1.GoToNext()



# save screenshot
SaveScreenshot(outputDir + 'time1s.png', renderView1, ImageResolution=[1481, 786],
    FontScaling='Scale fonts proportionally',
    OverrideColorPalette='WhiteBackground',
    StereoMode='No change',
    TransparentBackground=1,
    ImageQuality=100)

animationScene1.GoToNext()

animationScene1.GoToNext()

# save screenshot
SaveScreenshot(outputDir + 'time2s.png', renderView1, ImageResolution=[1481, 786],
    FontScaling='Scale fonts proportionally',
    OverrideColorPalette='WhiteBackground',
    StereoMode='No change',
    TransparentBackground=1,
    ImageQuality=100)

animationScene1.GoToNext()

animationScene1.GoToNext()

# save screenshot
SaveScreenshot(outputDir + 'time3s.png', renderView1, ImageResolution=[1481, 786],
    FontScaling='Scale fonts proportionally',
    OverrideColorPalette='WhiteBackground',
    StereoMode='No change',
    TransparentBackground=1,
    ImageQuality=100)

animationScene1.GoToNext()

animationScene1.GoToNext()

# save screenshot
SaveScreenshot(outputDir + 'time4s.png', renderView1, ImageResolution=[1481, 786],
    FontScaling='Scale fonts proportionally',
    OverrideColorPalette='WhiteBackground',
    StereoMode='No change',
    TransparentBackground=1,
    ImageQuality=100)

animationScene1.GoToNext()

animationScene1.GoToNext()

# save screenshot
SaveScreenshot(outputDir + 'time5s.png', renderView1, ImageResolution=[1481, 786],
    FontScaling='Scale fonts proportionally',
    OverrideColorPalette='WhiteBackground',
    StereoMode='No change',
    TransparentBackground=1,
    ImageQuality=100)

animationScene1.GoToNext()

animationScene1.GoToNext()

# save screenshot
SaveScreenshot(outputDir + 'time6s.png', renderView1, ImageResolution=[1481, 786],
    FontScaling='Scale fonts proportionally',
    OverrideColorPalette='WhiteBackground',
    StereoMode='No change',
    TransparentBackground=1,
    ImageQuality=100)

animationScene1.GoToNext()

animationScene1.GoToNext()

# save screenshot
SaveScreenshot(outputDir + 'time7s.png', renderView1, ImageResolution=[1481, 786],
    FontScaling='Scale fonts proportionally',
    OverrideColorPalette='WhiteBackground',
    StereoMode='No change',
    TransparentBackground=1,
    ImageQuality=100)

animationScene1.GoToNext()

animationScene1.GoToNext()

# save screenshot
SaveScreenshot(outputDir + 'time8s.png', renderView1, ImageResolution=[1481, 786],
    FontScaling='Scale fonts proportionally',
    OverrideColorPalette='WhiteBackground',
    StereoMode='No change',
    TransparentBackground=1,
    ImageQuality=100)

animationScene1.GoToNext()

animationScene1.GoToNext()

# save screenshot
SaveScreenshot(outputDir + 'time9s.png', renderView1, ImageResolution=[1481, 786],
    FontScaling='Scale fonts proportionally',
    OverrideColorPalette='WhiteBackground',
    StereoMode='No change',
    TransparentBackground=1,
    ImageQuality=100)

animationScene1.GoToNext()

animationScene1.GoToNext()


# save screenshot
SaveScreenshot(outputDir + 'time10s.png', renderView1, ImageResolution=[1481, 786],
    FontScaling='Scale fonts proportionally',
    OverrideColorPalette='WhiteBackground',
    StereoMode='No change',
    TransparentBackground=1,
    ImageQuality=100)


# set active source
SetActiveSource(velo)

# show data in view
veloDisplay = Show(velo, renderView1)

# show color bar/color legend
veloDisplay.SetScalarBarVisibility(renderView1, True)

# reset view to fit data
renderView1.ResetCamera()

# set active source
SetActiveSource(glyph1)

# show data in view
glyph1Display = Show(glyph1, renderView1)

# show color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, True)

# get color legend/bar for glyphVectorLUT in view renderView1
glyphVectorLUTColorBar = GetScalarBar(glyphVectorLUT, renderView1)

# change scalar bar placement
glyphVectorLUTColorBar.WindowLocation = 'AnyLocation'
glyphVectorLUTColorBar.Position = [0.14947670492910192, 0.07760814249363868]
glyphVectorLUTColorBar.ScalarBarLength = 0.3299999999999999

# get color legend/bar for vLUT in view renderView1
vLUTColorBar = GetScalarBar(vLUT, renderView1)

# change scalar bar placement
vLUTColorBar.Position = [0.8366576637407157, 0.04417302798982188]
vLUTColorBar.ScalarBarLength = 0.33

# Properties modified on glyphVectorLUTColorBar
glyphVectorLUTColorBar.AutoOrient = 0
glyphVectorLUTColorBar.Orientation = 'Horizontal'
glyphVectorLUTColorBar.ScalarBarLength = 0.33

# set active source
SetActiveSource(velo)

# Properties modified on vLUTColorBar
vLUTColorBar.AutoOrient = 0
vLUTColorBar.Orientation = 'Horizontal'

# Properties modified on nLUTColorBar
vLUTColorBar.Title = '|v| [m/s]'
vLUTColorBar.ComponentTitle = ''
vLUTColorBar.AutomaticLabelFormat = 0
vLUTColorBar.LabelFormat = '%-#6.1f'
vLUTColorBar.RangeLabelFormat = '%-#6.1f'

# change scalar bar placement
vLUTColorBar.Position = [0.6216137744767055, 0.33027989821882897]
vLUTColorBar.ScalarBarLength = 0.33000000000000007

# change scalar bar placement
glyphVectorLUTColorBar.Position = [0.0742032410533423, 0.3153689567430026]

# change scalar bar placement
glyphVectorLUTColorBar.Position = [0.08365631330182305, 0.32045801526717566]

# change scalar bar placement
glyphVectorLUTColorBar.Position = [0.08298109385550301, 0.386615776081425]

# Properties modified on nLUTColorBar
glyphVectorLUTColorBar.Title = '|n|'
glyphVectorLUTColorBar.ComponentTitle = ''
glyphVectorLUTColorBar.AutomaticLabelFormat = 0
glyphVectorLUTColorBar.LabelFormat = '%-#6.1f'
glyphVectorLUTColorBar.RangeLabelFormat = '%-#6.1f'

# change scalar bar placement
vLUTColorBar.Position = [0.614861580013505, 0.3824427480916025]

# get layout
layout1 = GetLayout()

# split cell
layout1.SplitVertical(0, 0.5)

# set active view
SetActiveView(None)

# set active view
SetActiveView(renderView1)

# resize frame
layout1.SetSplitFraction(0, 0.400739827374)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [1.095491848901468, -1.4263053347930088, 4.319751617610021]
renderView1.CameraFocalPoint = [1.095491848901468, -1.4263053347930088, 0.0]
renderView1.CameraParallelScale = 0.706124836011424

# save screenshot
SaveScreenshot(outputDir + 'colorBars.png', renderView1, ImageResolution=[1481, 297],
    OverrideColorPalette='WhiteBackground',
    TransparentBackground=1)

#### saving camera placements for all active views


#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
