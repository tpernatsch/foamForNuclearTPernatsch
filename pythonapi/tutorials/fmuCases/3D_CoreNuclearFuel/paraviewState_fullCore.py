# state file generated using paraview version 5.10.0-RC1

# uncomment the following three lines to ensure this script works in future versions
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 10

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# folderName = "2024-02-29_SUCCESS_steadyState_144pins_fullCore"
folderName = "steadyState"

def formatFromHoneycomb(data: str):
    xPos, yPos = [], []
    for line in data.split('\n'):
        d = line.split(";")
        if (len(d) >= 3):
            if (d[2] != "e"):
                xPos.append(float(d[0]))
                yPos.append(float(d[1]))
    return(xPos, yPos)

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1411, 784]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.CenterOfRotation = [0.0, 0.0, 0.0]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [2.976275972667933, -4.43997884037497, 3.6681641557451075]
renderView1.CameraFocalPoint = [0.0, 0.0, 0.0]
renderView1.CameraViewUp = [-0.3091357913223522, 0.4834560761173773, 0.8189659852452309]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 1.6891352608047858

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='GeN-Foam + 144 OFFBEAT')
layout1.AssignView(0, renderView1)
layout1.SetSize(1411, 784)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'OpenFOAMReader'
corefoam = OpenFOAMReader(registrationName='Core.foam', FileName=f'/home/thomas-guilbaud/Simulation/GeN-Foam/3D_CorePinsCouplingGFOB/{folderName}/Core/Core.foam')
corefoam.MeshRegions = ['/fluidRegion/internalMesh']
corefoam.CellArrays = ['CRDisp', 'E', 'Kd', 'Re', 'T', 'T.cladAvForNeutronics', 'T.fuelAvForNeutronics', 'TClad', 'TCool', 'TFuel', 'TStruct', 'TStructFromTH', 'TStructMech', 'Tsurface.nuclearFuelFMU', 'U', 'alpha', 'alpha.passiveStructure', 'alpha.powerModel', 'alphaRhoMagU', 'defaultExternalSourceFlux', 'defaultFlux', 'defaultPrec', 'dgdt', 'disp', 'externalSourceFlux0', 'externalSourceFlux1', 'flux0', 'flux1', 'fluxStar0', 'fluxStar1', 'fuelDisp', 'fuelDispVector', 'gapWidth', 'heatFlux.structure', 'htc', 'isPorous', 'meshDisp', 'mu', 'nu', 'oneGroupFlux', 'p', 'p_rgh', 'powerDensity', 'powerDensityNeutronics', 'powerDensityNeutronicsToLiquid', 'prec0', 'prec1', 'prec2', 'prec3', 'prec4', 'prec5', 'precStar0', 'precStar1', 'precStar2', 'precStar3', 'precStar4', 'precStar5', 'rho', 'rhoCool', 'secondaryPowerDenisty']

# create a new 'Plot Over Line'
corePlotOverLineZ = PlotOverLine(registrationName='CorePlotOverLineZ', Input=corefoam)
corePlotOverLineZ.Point1 = [-0.723, 0.083, -1] # 0.0]
corePlotOverLineZ.Point2 = [-0.723, 0.083, 1] # 1.8473526239395142]

# create a new 'Slice'
coreSlice = Slice(registrationName='CoreSlice', Input=corefoam)
coreSlice.SliceType = 'Plane'
coreSlice.HyperTreeGridSlicer = 'Plane'
coreSlice.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
coreSlice.SliceType.Origin = [0.0, 0.0, 0.0]
coreSlice.SliceType.Normal = [0.0, 0.0, 1.0]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
coreSlice.HyperTreeGridSlicer.Origin = [0.0, 0.0, 0.0]

# create a new 'Warp By Scalar'
coreWarpByScalar1 = WarpByScalar(registrationName='CoreWarpByScalar1', Input=coreSlice)
coreWarpByScalar1.Scalars = ['POINTS', 'Tsurface.nuclearFuelFMU']
coreWarpByScalar1.ScaleFactor = 0.0012009159329527696
coreWarpByScalar1.UseNormal = 1

# create a new 'Transform'
coreTransform = Transform(registrationName='CoreTransform', Input=coreWarpByScalar1)
coreTransform.Transform = 'Transform'

# init the 'Transform' selected for 'Transform'
coreTransform.Transform.Translate = [0.0, 0.0, -0.5]

# create a new 'Plot Over Line'
corePlotOverLineX = PlotOverLine(registrationName='CorePlotOverLineX', Input=corefoam)
corePlotOverLineX.Point1 = [-1.0, 0.0, 0]
corePlotOverLineX.Point2 = [1.0, 0.0, 0]

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from corefoam
corefoamDisplay = Show(corefoam, renderView1, 'UnstructuredGridRepresentation')

# get color transfer function/color map for 'p'
pLUT = GetColorTransferFunction('p')
pLUT.RGBPoints = [81877.46875, 0.231373, 0.298039, 0.752941, 91615.3359375, 0.865003, 0.865003, 0.865003, 101353.203125, 0.705882, 0.0156863, 0.14902]
pLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'p'
pPWF = GetOpacityTransferFunction('p')
pPWF.Points = [81877.46875, 0.0, 0.5, 0.0, 101353.203125, 1.0, 0.5, 0.0]
pPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
corefoamDisplay.Representation = 'Surface'
corefoamDisplay.ColorArrayName = ['POINTS', '']
corefoamDisplay.LookupTable = pLUT
corefoamDisplay.Opacity = 0.2
corefoamDisplay.SelectTCoordArray = 'None'
corefoamDisplay.SelectNormalArray = 'None'
corefoamDisplay.SelectTangentArray = 'None'
corefoamDisplay.OSPRayScaleArray = 'p'
corefoamDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
corefoamDisplay.SelectOrientationVectors = 'U'
corefoamDisplay.ScaleFactor = 0.2
corefoamDisplay.SelectScaleArray = 'p'
corefoamDisplay.GlyphType = 'Arrow'
corefoamDisplay.GlyphTableIndexArray = 'p'
corefoamDisplay.GaussianRadius = 0.01
corefoamDisplay.SetScaleArray = ['POINTS', 'p']
corefoamDisplay.ScaleTransferFunction = 'PiecewiseFunction'
corefoamDisplay.OpacityArray = ['POINTS', 'p']
corefoamDisplay.OpacityTransferFunction = 'PiecewiseFunction'
corefoamDisplay.DataAxesGrid = 'GridAxesRepresentation'
corefoamDisplay.PolarAxes = 'PolarAxesRepresentation'
corefoamDisplay.ScalarOpacityFunction = pPWF
corefoamDisplay.ScalarOpacityUnitDistance = 0.10231208221624212
corefoamDisplay.OpacityArrayName = ['POINTS', 'p']

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
corefoamDisplay.ScaleTransferFunction.Points = [81877.46875, 0.0, 0.5, 0.0, 101353.203125, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
corefoamDisplay.OpacityTransferFunction.Points = [81877.46875, 0.0, 0.5, 0.0, 101353.203125, 1.0, 0.5, 0.0]

# show data from coreTransform
coreTransformDisplay = Show(coreTransform, renderView1, 'GeometryRepresentation')

# get color transfer function/color map for 'TsurfacenuclearFuelFMU'
tsurfacenuclearFuelFMULUT = GetColorTransferFunction('TsurfacenuclearFuelFMU')
tsurfacenuclearFuelFMULUT.RGBPoints = [300, 0.231373, 0.298039, 0.752941, 548.8034820556641, 0.865003, 0.865003, 0.865003, 800.1378784179688, 0.705882, 0.0156863, 0.14902]
tsurfacenuclearFuelFMULUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
coreTransformDisplay.Representation = 'Surface'
coreTransformDisplay.ColorArrayName = ['CELLS', 'Tsurface.nuclearFuelFMU']
coreTransformDisplay.LookupTable = tsurfacenuclearFuelFMULUT
coreTransformDisplay.SelectTCoordArray = 'None'
coreTransformDisplay.SelectNormalArray = 'None'
coreTransformDisplay.SelectTangentArray = 'None'
coreTransformDisplay.OSPRayScaleArray = 'p'
coreTransformDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
coreTransformDisplay.SelectOrientationVectors = 'U'
coreTransformDisplay.ScaleFactor = 0.2
coreTransformDisplay.SelectScaleArray = 'p'
coreTransformDisplay.GlyphType = 'Arrow'
coreTransformDisplay.GlyphTableIndexArray = 'p'
coreTransformDisplay.GaussianRadius = 0.01
coreTransformDisplay.SetScaleArray = ['POINTS', 'p']
coreTransformDisplay.ScaleTransferFunction = 'PiecewiseFunction'
coreTransformDisplay.OpacityArray = ['POINTS', 'p']
coreTransformDisplay.OpacityTransferFunction = 'PiecewiseFunction'
coreTransformDisplay.DataAxesGrid = 'GridAxesRepresentation'
coreTransformDisplay.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
coreTransformDisplay.ScaleTransferFunction.Points = [91615.34375, 0.0, 0.5, 0.0, 91631.34375, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
coreTransformDisplay.OpacityTransferFunction.Points = [91615.34375, 0.0, 0.5, 0.0, 91631.34375, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for tsurfacenuclearFuelFMULUT in view renderView1
tsurfacenuclearFuelFMULUTColorBar = GetScalarBar(tsurfacenuclearFuelFMULUT, renderView1)
tsurfacenuclearFuelFMULUTColorBar.WindowLocation = 'Upper Right Corner'
tsurfacenuclearFuelFMULUTColorBar.Title = 'Tsurface.nuclearFuelFMU'
tsurfacenuclearFuelFMULUTColorBar.ComponentTitle = ''

# set color bar visibility
tsurfacenuclearFuelFMULUTColorBar.Visibility = 1

# show color legend
coreTransformDisplay.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'TsurfacenuclearFuelFMU'
tsurfacenuclearFuelFMUPWF = GetOpacityTransferFunction('TsurfacenuclearFuelFMU')
tsurfacenuclearFuelFMUPWF.Points = [300, 0.0, 0.5, 0.0, 800.1378784179688, 1.0, 0.5, 0.0]
tsurfacenuclearFuelFMUPWF.ScalarRangeInitialized = 1

# ----------------------------------------------------------------
# restore active source
SetActiveSource(corePlotOverLineX)
# ----------------------------------------------------------------

# Inner fuel
xPosInner, yPosInner = formatFromHoneycomb(
    """-1.012383697024009;-1.7535;e
-1.012383697024009;-1.5865;e
-1.012383697024009;-1.4195000000000002;e
-1.012383697024009;-1.2525;e
-1.012383697024009;-1.0855000000000001;e
-1.012383697024009;-0.9185000000000001;e
-1.012383697024009;-0.7515000000000001;e
-1.012383697024009;-0.5845;e
-1.012383697024009;-0.4175;e
-1.012383697024009;-0.2505;e
-1.012383697024009;-0.08350000000000002;e
-1.012383697024009;0.08350000000000002;e
-1.012383697024009;0.25050000000000006;e
-1.012383697024009;0.4175;e
-1.012383697024009;0.5845;e
-0.8677574545920077;-1.67;e
-0.8677574545920077;-1.5030000000000001;e
-0.8677574545920077;-1.336;e
-0.8677574545920077;-1.169;e
-0.8677574545920077;-1.002;e
-0.8677574545920077;-0.835;e
-0.8677574545920077;-0.668;e
-0.8677574545920077;-0.501;e
-0.8677574545920077;-0.33399999999999996;e
-0.8677574545920077;-0.16699999999999998;e
-0.8677574545920077;0;e
-0.8677574545920077;0.16700000000000004;e
-0.8677574545920077;0.3340000000000001;e
-0.8677574545920077;0.501;e
-0.8677574545920077;0.668;e
-0.7231312121600064;-1.5865;e
-0.7231312121600064;-1.4195;e
-0.7231312121600064;-1.2525000000000002;e
-0.7231312121600064;-1.0855000000000001;e
-0.7231312121600064;-0.9185000000000001;e
-0.7231312121600064;-0.7515000000000001;e
-0.7231312121600064;-0.5845;e
-0.7231312121600064;-0.41750000000000004;e
-0.7231312121600064;-0.25050000000000006;e
-0.7231312121600064;-0.08350000000000002;e
-0.7231312121600064;0.08349999999999996;e
-0.7231312121600064;0.2505;e
-0.7231312121600064;0.41750000000000004;e
-0.7231312121600064;0.5845;e
-0.7231312121600064;0.7515000000000001;e
-0.5785049697280051;-1.5030000000000001;e
-0.5785049697280051;-1.336;e
-0.5785049697280051;-1.169;e
-0.5785049697280051;-1.002;e
-0.5785049697280051;-0.835;e
-0.5785049697280051;-0.668;e
-0.5785049697280051;-0.501;e
-0.5785049697280051;-0.334;F
-0.5785049697280051;-0.167;F
-0.5785049697280051;0;F
-0.5785049697280051;0.16699999999999998;F
-0.5785049697280051;0.334;F
-0.5785049697280051;0.5010000000000001;e
-0.5785049697280051;0.6679999999999999;e
-0.5785049697280051;0.835;e
-0.43387872729600385;-1.4195;e
-0.43387872729600385;-1.2525;e
-0.43387872729600385;-1.0855000000000001;e
-0.43387872729600385;-0.9185000000000001;e
-0.43387872729600385;-0.7515000000000001;e
-0.43387872729600385;-0.5845;e
-0.43387872729600385;-0.4175;F
-0.43387872729600385;-0.2505;F
-0.43387872729600385;-0.08349999999999999;F
-0.43387872729600385;0.08350000000000002;F
-0.43387872729600385;0.2505;F
-0.43387872729600385;0.41750000000000004;F
-0.43387872729600385;0.5845;e
-0.43387872729600385;0.7515000000000001;e
-0.43387872729600385;0.9185000000000001;e
-0.28925248486400257;-1.336;e
-0.28925248486400257;-1.169;e
-0.28925248486400257;-1.002;e
-0.28925248486400257;-0.8350000000000001;e
-0.28925248486400257;-0.668;e
-0.28925248486400257;-0.501;F
-0.28925248486400257;-0.334;F
-0.28925248486400257;-0.167;F
-0.28925248486400257;0;F
-0.28925248486400257;0.167;F
-0.28925248486400257;0.33399999999999996;F
-0.28925248486400257;0.501;F
-0.28925248486400257;0.668;e
-0.28925248486400257;0.835;e
-0.28925248486400257;1.002;e
-0.14462624243200128;-1.2525;e
-0.14462624243200128;-1.0855;e
-0.14462624243200128;-0.9185000000000001;e
-0.14462624243200128;-0.7515000000000001;e
-0.14462624243200128;-0.5845;F
-0.14462624243200128;-0.41750000000000004;F
-0.14462624243200128;-0.2505;F
-0.14462624243200128;-0.0835;e
-0.14462624243200128;0.0835;e
-0.14462624243200128;0.2505;F
-0.14462624243200128;0.4175;F
-0.14462624243200128;0.5845;F
-0.14462624243200128;0.7515000000000001;e
-0.14462624243200128;0.9185;e
-0.14462624243200128;1.0855000000000001;e
0;-1.169;e
0;-1.002;e
0;-0.8350000000000001;e
0;-0.668;F
0;-0.501;F
0;-0.334;F
0;-0.167;e
0;0;e
0;0.167;e
0;0.334;F
0;0.501;F
0;0.668;F
0;0.8350000000000001;e
0;1.002;e
0;1.169;e
0.14462624243200128;-1.0855000000000001;e
0.14462624243200128;-0.9185;e
0.14462624243200128;-0.7515000000000001;e
0.14462624243200128;-0.5845;F
0.14462624243200128;-0.4175;F
0.14462624243200128;-0.2505;F
0.14462624243200128;-0.0835;e
0.14462624243200128;0.0835;e
0.14462624243200128;0.2505;F
0.14462624243200128;0.41750000000000004;F
0.14462624243200128;0.5845;F
0.14462624243200128;0.7515000000000001;e
0.14462624243200128;0.9185000000000001;e
0.14462624243200128;1.0855;e
0.14462624243200128;1.2525;e
0.28925248486400257;-1.002;e
0.28925248486400257;-0.835;e
0.28925248486400257;-0.668;e
0.28925248486400257;-0.501;F
0.28925248486400257;-0.33399999999999996;F
0.28925248486400257;-0.167;F
0.28925248486400257;0;F
0.28925248486400257;0.167;F
0.28925248486400257;0.334;F
0.28925248486400257;0.501;F
0.28925248486400257;0.668;e
0.28925248486400257;0.8350000000000001;e
0.28925248486400257;1.002;e
0.28925248486400257;1.169;e
0.28925248486400257;1.336;e
0.43387872729600385;-0.9185000000000001;e
0.43387872729600385;-0.7515000000000001;e
0.43387872729600385;-0.5845;e
0.43387872729600385;-0.41750000000000004;F
0.43387872729600385;-0.2505;F
0.43387872729600385;-0.08350000000000002;F
0.43387872729600385;0.08349999999999999;F
0.43387872729600385;0.2505;F
0.43387872729600385;0.4175;F
0.43387872729600385;0.5845;e
0.43387872729600385;0.7515000000000001;e
0.43387872729600385;0.9185000000000001;e
0.43387872729600385;1.0855000000000001;e
0.43387872729600385;1.2525;e
0.43387872729600385;1.4195;e
0.5785049697280051;-0.835;e
0.5785049697280051;-0.6679999999999999;e
0.5785049697280051;-0.5010000000000001;e
0.5785049697280051;-0.334;F
0.5785049697280051;-0.16699999999999998;F
0.5785049697280051;0;F
0.5785049697280051;0.167;F
0.5785049697280051;0.334;F
0.5785049697280051;0.501;e
0.5785049697280051;0.668;e
0.5785049697280051;0.835;e
0.5785049697280051;1.002;e
0.5785049697280051;1.169;e
0.5785049697280051;1.336;e
0.5785049697280051;1.5030000000000001;e
0.7231312121600064;-0.7515000000000001;e
0.7231312121600064;-0.5845;e
0.7231312121600064;-0.41750000000000004;e
0.7231312121600064;-0.2505;e
0.7231312121600064;-0.08349999999999996;e
0.7231312121600064;0.08350000000000002;e
0.7231312121600064;0.25050000000000006;e
0.7231312121600064;0.41750000000000004;e
0.7231312121600064;0.5845;e
0.7231312121600064;0.7515000000000001;e
0.7231312121600064;0.9185000000000001;e
0.7231312121600064;1.0855000000000001;e
0.7231312121600064;1.2525000000000002;e
0.7231312121600064;1.4195;e
0.7231312121600064;1.5865;e
0.8677574545920077;-0.668;e
0.8677574545920077;-0.501;e
0.8677574545920077;-0.3340000000000001;e
0.8677574545920077;-0.16700000000000004;e
0.8677574545920077;0;e
0.8677574545920077;0.16699999999999998;e
0.8677574545920077;0.33399999999999996;e
0.8677574545920077;0.501;e
0.8677574545920077;0.668;e
0.8677574545920077;0.835;e
0.8677574545920077;1.002;e
0.8677574545920077;1.169;e
0.8677574545920077;1.336;e
0.8677574545920077;1.5030000000000001;e
0.8677574545920077;1.67;e
1.012383697024009;-0.5845;e
1.012383697024009;-0.4175;e
1.012383697024009;-0.25050000000000006;e
1.012383697024009;-0.08350000000000002;e
1.012383697024009;0.08350000000000002;e
1.012383697024009;0.2505;e
1.012383697024009;0.4175;e
1.012383697024009;0.5845;e
1.012383697024009;0.7515000000000001;e
1.012383697024009;0.9185000000000001;e
1.012383697024009;1.0855000000000001;e
1.012383697024009;1.2525;e
1.012383697024009;1.4195000000000002;e
1.012383697024009;1.5865;e
1.012383697024009;1.7535;e"""
)

# Outer fuel
xPosOuter, yPosOuter = formatFromHoneycomb(
    """-1.012383697024009;-1.7535;e
-1.012383697024009;-1.5865;e
-1.012383697024009;-1.4195000000000002;e
-1.012383697024009;-1.2525;e
-1.012383697024009;-1.0855000000000001;e
-1.012383697024009;-0.9185000000000001;e
-1.012383697024009;-0.7515000000000001;e
-1.012383697024009;-0.5845;e
-1.012383697024009;-0.4175;e
-1.012383697024009;-0.2505;F
-1.012383697024009;-0.08350000000000002;F
-1.012383697024009;0.08350000000000002;F
-1.012383697024009;0.25050000000000006;F
-1.012383697024009;0.4175;e
-1.012383697024009;0.5845;e
-0.8677574545920077;-1.67;e
-0.8677574545920077;-1.5030000000000001;e
-0.8677574545920077;-1.336;e
-0.8677574545920077;-1.169;e
-0.8677574545920077;-1.002;e
-0.8677574545920077;-0.835;e
-0.8677574545920077;-0.668;e
-0.8677574545920077;-0.501;F
-0.8677574545920077;-0.33399999999999996;F
-0.8677574545920077;-0.16699999999999998;F
-0.8677574545920077;0;F
-0.8677574545920077;0.16700000000000004;F
-0.8677574545920077;0.3340000000000001;F
-0.8677574545920077;0.501;F
-0.8677574545920077;0.668;e
-0.7231312121600064;-1.5865;e
-0.7231312121600064;-1.4195;e
-0.7231312121600064;-1.2525000000000002;e
-0.7231312121600064;-1.0855000000000001;e
-0.7231312121600064;-0.9185000000000001;e
-0.7231312121600064;-0.7515000000000001;F
-0.7231312121600064;-0.5845;F
-0.7231312121600064;-0.41750000000000004;F
-0.7231312121600064;-0.25050000000000006;F
-0.7231312121600064;-0.08350000000000002;F
-0.7231312121600064;0.08349999999999996;F
-0.7231312121600064;0.2505;F
-0.7231312121600064;0.41750000000000004;F
-0.7231312121600064;0.5845;F
-0.7231312121600064;0.7515000000000001;F
-0.5785049697280051;-1.5030000000000001;e
-0.5785049697280051;-1.336;e
-0.5785049697280051;-1.169;e
-0.5785049697280051;-1.002;e
-0.5785049697280051;-0.835;F
-0.5785049697280051;-0.668;F
-0.5785049697280051;-0.501;F
-0.5785049697280051;-0.334;e
-0.5785049697280051;-0.167;e
-0.5785049697280051;0;e
-0.5785049697280051;0.16699999999999998;e
-0.5785049697280051;0.334;e
-0.5785049697280051;0.5010000000000001;F
-0.5785049697280051;0.6679999999999999;F
-0.5785049697280051;0.835;F
-0.43387872729600385;-1.4195;e
-0.43387872729600385;-1.2525;e
-0.43387872729600385;-1.0855000000000001;e
-0.43387872729600385;-0.9185000000000001;F
-0.43387872729600385;-0.7515000000000001;F
-0.43387872729600385;-0.5845;F
-0.43387872729600385;-0.4175;e
-0.43387872729600385;-0.2505;e
-0.43387872729600385;-0.08349999999999999;e
-0.43387872729600385;0.08350000000000002;e
-0.43387872729600385;0.2505;e
-0.43387872729600385;0.41750000000000004;e
-0.43387872729600385;0.5845;F
-0.43387872729600385;0.7515000000000001;F
-0.43387872729600385;0.9185000000000001;F
-0.28925248486400257;-1.336;e
-0.28925248486400257;-1.169;e
-0.28925248486400257;-1.002;F
-0.28925248486400257;-0.8350000000000001;F
-0.28925248486400257;-0.668;F
-0.28925248486400257;-0.501;e
-0.28925248486400257;-0.334;e
-0.28925248486400257;-0.167;e
-0.28925248486400257;0;e
-0.28925248486400257;0.167;e
-0.28925248486400257;0.33399999999999996;e
-0.28925248486400257;0.501;e
-0.28925248486400257;0.668;F
-0.28925248486400257;0.835;F
-0.28925248486400257;1.002;F
-0.14462624243200128;-1.2525;e
-0.14462624243200128;-1.0855;e
-0.14462624243200128;-0.9185000000000001;F
-0.14462624243200128;-0.7515000000000001;F
-0.14462624243200128;-0.5845;e
-0.14462624243200128;-0.41750000000000004;e
-0.14462624243200128;-0.2505;e
-0.14462624243200128;-0.0835;e
-0.14462624243200128;0.0835;e
-0.14462624243200128;0.2505;e
-0.14462624243200128;0.4175;e
-0.14462624243200128;0.5845;e
-0.14462624243200128;0.7515000000000001;F
-0.14462624243200128;0.9185;F
-0.14462624243200128;1.0855000000000001;e
0;-1.169;e
0;-1.002;F
0;-0.8350000000000001;F
0;-0.668;e
0;-0.501;e
0;-0.334;e
0;-0.167;e
0;0;e
0;0.167;e
0;0.334;e
0;0.501;e
0;0.668;e
0;0.8350000000000001;F
0;1.002;F
0;1.169;e
0.14462624243200128;-1.0855000000000001;e
0.14462624243200128;-0.9185;F
0.14462624243200128;-0.7515000000000001;F
0.14462624243200128;-0.5845;e
0.14462624243200128;-0.4175;e
0.14462624243200128;-0.2505;e
0.14462624243200128;-0.0835;e
0.14462624243200128;0.0835;e
0.14462624243200128;0.2505;e
0.14462624243200128;0.41750000000000004;e
0.14462624243200128;0.5845;e
0.14462624243200128;0.7515000000000001;F
0.14462624243200128;0.9185000000000001;F
0.14462624243200128;1.0855;e
0.14462624243200128;1.2525;e
0.28925248486400257;-1.002;F
0.28925248486400257;-0.835;F
0.28925248486400257;-0.668;F
0.28925248486400257;-0.501;e
0.28925248486400257;-0.33399999999999996;e
0.28925248486400257;-0.167;e
0.28925248486400257;0;e
0.28925248486400257;0.167;e
0.28925248486400257;0.334;e
0.28925248486400257;0.501;e
0.28925248486400257;0.668;F
0.28925248486400257;0.8350000000000001;F
0.28925248486400257;1.002;F
0.28925248486400257;1.169;e
0.28925248486400257;1.336;e
0.43387872729600385;-0.9185000000000001;F
0.43387872729600385;-0.7515000000000001;F
0.43387872729600385;-0.5845;F
0.43387872729600385;-0.41750000000000004;e
0.43387872729600385;-0.2505;e
0.43387872729600385;-0.08350000000000002;e
0.43387872729600385;0.08349999999999999;e
0.43387872729600385;0.2505;e
0.43387872729600385;0.4175;e
0.43387872729600385;0.5845;F
0.43387872729600385;0.7515000000000001;F
0.43387872729600385;0.9185000000000001;F
0.43387872729600385;1.0855000000000001;e
0.43387872729600385;1.2525;e
0.43387872729600385;1.4195;e
0.5785049697280051;-0.835;F
0.5785049697280051;-0.6679999999999999;F
0.5785049697280051;-0.5010000000000001;F
0.5785049697280051;-0.334;e
0.5785049697280051;-0.16699999999999998;e
0.5785049697280051;0;e
0.5785049697280051;0.167;e
0.5785049697280051;0.334;e
0.5785049697280051;0.501;F
0.5785049697280051;0.668;F
0.5785049697280051;0.835;F
0.5785049697280051;1.002;e
0.5785049697280051;1.169;e
0.5785049697280051;1.336;e
0.5785049697280051;1.5030000000000001;e
0.7231312121600064;-0.7515000000000001;F
0.7231312121600064;-0.5845;F
0.7231312121600064;-0.41750000000000004;F
0.7231312121600064;-0.2505;F
0.7231312121600064;-0.08349999999999996;F
0.7231312121600064;0.08350000000000002;F
0.7231312121600064;0.25050000000000006;F
0.7231312121600064;0.41750000000000004;F
0.7231312121600064;0.5845;F
0.7231312121600064;0.7515000000000001;F
0.7231312121600064;0.9185000000000001;e
0.7231312121600064;1.0855000000000001;e
0.7231312121600064;1.2525000000000002;e
0.7231312121600064;1.4195;e
0.7231312121600064;1.5865;e
0.8677574545920077;-0.668;e
0.8677574545920077;-0.501;F
0.8677574545920077;-0.3340000000000001;F
0.8677574545920077;-0.16700000000000004;F
0.8677574545920077;0;F
0.8677574545920077;0.16699999999999998;F
0.8677574545920077;0.33399999999999996;F
0.8677574545920077;0.501;F
0.8677574545920077;0.668;e
0.8677574545920077;0.835;e
0.8677574545920077;1.002;e
0.8677574545920077;1.169;e
0.8677574545920077;1.336;e
0.8677574545920077;1.5030000000000001;e
0.8677574545920077;1.67;e
1.012383697024009;-0.5845;e
1.012383697024009;-0.4175;e
1.012383697024009;-0.25050000000000006;F
1.012383697024009;-0.08350000000000002;F
1.012383697024009;0.08350000000000002;F
1.012383697024009;0.2505;F
1.012383697024009;0.4175;e
1.012383697024009;0.5845;e
1.012383697024009;0.7515000000000001;e
1.012383697024009;0.9185000000000001;e
1.012383697024009;1.0855000000000001;e
1.012383697024009;1.2525;e
1.012383697024009;1.4195000000000002;e
1.012383697024009;1.5865;e
1.012383697024009;1.7535;e"""
)


xPos, yPos = xPosInner+xPosOuter, yPosInner+yPosOuter

for i in range(len(xPos)):
    # create a new 'OpenFOAMReader'
    fuelPinfoam = OpenFOAMReader(registrationName='FuelPin.foam', FileName=f'/home/thomas-guilbaud/Simulation/GeN-Foam/3D_CorePinsCouplingGFOB/{folderName}/FuelPin{i}/FuelPin.foam')
    fuelPinfoam.MeshRegions = ['internalMesh']
    # fuelPinfoam.CellArrays = ['Bu', 'Bu_kgU', 'DCyl', 'DEpsilonP', 'DEpsilonPEq', 'DLambda', 'DSigmaY', 'DepsilonCreep', 'DepsilonCreepEq', 'DepsilonCreepIrrEq', 'DepsilonCreepPrimEq', 'DepsilonCreepThEq', 'N_Am241', 'N_Am243', 'N_Cm242', 'N_Cm243', 'N_Cm244', 'N_Np237', 'N_Pu239', 'N_Pu240', 'N_Pu241', 'N_Pu242', 'N_U235', 'N_U236', 'N_U238', 'Q', 'T', 'activeYield', 'betaFraction', 'divSigmaExp', 'epsEl', 'epsTh', 'epsilon', 'epsilonCreep', 'epsilonCreepEq', 'epsilonCreepIrr', 'epsilonCreepIrrEq', 'epsilonCreepPrim', 'epsilonCreepPrimEq', 'epsilonCreepTh', 'epsilonCreepThEq', 'epsilonCyl', 'epsilonDensification', 'epsilonP', 'epsilonPEq', 'epsilonRecoveredRelocation', 'epsilonRelocation', 'epsilonSwelling', 'fastFluence', 'fastFlux', 'formFactor', 'gapWidth', 'gradD', 'gradT', 'grainRadius', 'interfaceP', 'intergranularGasSwelling', 'intragranularGasSwelling', 'lhgr', 'nCracks', 'neutronFlux0', 'oxygenMetalRatio', 'plasticN', 'poreVelocity', 'porosity', 'sigma', 'sigmaCyl', 'sigmaDev', 'sigmaEq', 'sigmaExp', 'sigmaHyd', 'sigmaY', 'sliceID']
    fuelPinfoam.CellArrays = ['DEpsilonP', 'DEpsilonPEq', 'DLambda', 'DSigmaY', 'DepsilonCreep', 'DepsilonCreepEq', 'DepsilonCreepIrrEq', 'DepsilonCreepPrimEq', 'DepsilonCreepThEq', 'Q', 'T', 'divSigmaExp', 'epsEl', 'epsTh', 'epsilon', 'epsilonCreep', 'epsilonCreepEq', 'epsilonCreepIrr', 'epsilonCreepIrrEq', 'epsilonCreepPrim', 'epsilonCreepPrimEq', 'epsilonCreepTh', 'epsilonCreepThEq', 'epsilonCyl', 'epsilonDensification', 'epsilonP', 'epsilonPEq', 'epsilonSwelling', 'fastFluence', 'fastFlux', 'gapWidth', 'gradD', 'gradT', 'grainRadius', 'interfaceP', 'intergranularGasSwelling', 'intragranularGasSwelling', 'lhgr', 'nCracks', 'neutronFlux0', 'plasticN', 'sigma', 'sigmaCyl', 'sigmaDev', 'sigmaEq', 'sigmaExp', 'sigmaHyd', 'sigmaY', 'sliceID']

    # create a new 'Transform'
    fuelPinFoamTransformed = Transform(registrationName=f'Transform{i}', Input=fuelPinfoam)
    fuelPinFoamTransformed.Transform = 'Transform'

    # init the 'Transform' selected for 'Transform'
    fuelPinFoamTransformed.Transform.Scale = [10.0, 10.0, 1.0]
    fuelPinFoamTransformed.Transform.Translate = [xPos[i], yPos[i], 0.0]

    # show data from corefoam
    fuelPinfoamDisplay = Show(fuelPinFoamTransformed, renderView1, 'UnstructuredGridRepresentation')

    # get color transfer function/color map for 'T'
    tLUT = GetColorTransferFunction('T')
    tLUT.RGBPoints = [300, 0.231373, 0.298039, 0.752941, 1369.2445526123047, 0.865003, 0.865003, 0.865003, 2441.02001953125, 0.705882, 0.0156863, 0.14902]
    tLUT.ScalarRangeInitialized = 1.0

    # get opacity transfer function/opacity map for 'T'
    tPWF = GetOpacityTransferFunction('T')
    tPWF.Points = [300, 0.0, 0.5, 0.0, 2441.02001953125, 1.0, 0.5, 0.0]
    tPWF.ScalarRangeInitialized = 1

    # trace defaults for the display properties.
    fuelPinfoamDisplay.Representation = 'Surface'
    fuelPinfoamDisplay.ColorArrayName = ['CELLS', 'T']
    fuelPinfoamDisplay.LookupTable = tLUT
    fuelPinfoamDisplay.SelectTCoordArray = 'None'
    fuelPinfoamDisplay.SelectNormalArray = 'None'
    fuelPinfoamDisplay.SelectTangentArray = 'None'
    fuelPinfoamDisplay.OSPRayScaleArray = 'Bu'
    fuelPinfoamDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    fuelPinfoamDisplay.SelectOrientationVectors = 'DCyl'
    fuelPinfoamDisplay.ScaleFactor = 0.1847000002861023
    fuelPinfoamDisplay.SelectScaleArray = 'Bu'
    fuelPinfoamDisplay.GlyphType = 'Arrow'
    fuelPinfoamDisplay.GlyphTableIndexArray = 'Bu'
    fuelPinfoamDisplay.GaussianRadius = 0.009235000014305115
    fuelPinfoamDisplay.SetScaleArray = ['POINTS', 'Bu']
    fuelPinfoamDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    fuelPinfoamDisplay.OpacityArray = ['POINTS', 'Bu']
    fuelPinfoamDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    fuelPinfoamDisplay.DataAxesGrid = 'GridAxesRepresentation'
    fuelPinfoamDisplay.PolarAxes = 'PolarAxesRepresentation'
    fuelPinfoamDisplay.ScalarOpacityFunction = tPWF
    fuelPinfoamDisplay.ScalarOpacityUnitDistance = 0.15594667949912208
    fuelPinfoamDisplay.OpacityArrayName = ['POINTS', 'Bu']

    # show color legend
    fuelPinfoamDisplay.SetScalarBarVisibility(renderView1, True)

# GetActiveView().UseOffscreenRenderingForScreenshots = 0
# GetActiveView().UseInteractiveRenderingForScreenshots = 0

# Show()
# Render()
# SaveScreenshot('greenSphereScreenshot.png', renderView1, ImageResolution=[1670, 1090])


if __name__ == '__main__':
    # generate extracts
    SaveExtracts(ExtractsOutputDirectory='extracts')
