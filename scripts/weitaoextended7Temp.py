# trace generated using paraview version 5.7.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'CSV Reader'
pararesults0 = CSVReader(FileName=['/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults001.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults002.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults003.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults004.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults005.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults006.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults007.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults008.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults009.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults010.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults011.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults012.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults013.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults014.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults015.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults016.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults017.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults018.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults019.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults020.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults021.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults022.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults023.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults024.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults025.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults026.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults027.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults028.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults029.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults030.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults031.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults032.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults033.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults034.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults035.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults036.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults037.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults038.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults039.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults040.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults041.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults042.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults043.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults044.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults045.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults046.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults047.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults048.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults049.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults050.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults051.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults052.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults053.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults054.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults055.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults056.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults057.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults058.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults059.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults060.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults061.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults062.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults063.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults064.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults065.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults066.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults067.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults068.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults069.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults070.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults071.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults072.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults073.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults074.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults075.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults076.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults077.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults078.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults079.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults080.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults081.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults082.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults083.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults084.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults085.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults086.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults087.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults088.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults089.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults090.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults091.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults092.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults093.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults094.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults095.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults096.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults097.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults098.csv', '/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Raw_DataFiles/wtresults021621DcTrWt/pararesults099.csv'])
pararesults0.DetectNumericColumns = 1
pararesults0.UseStringDelimiter = 1
pararesults0.HaveHeaders = 1
pararesults0.FieldDelimiterCharacters = ','
pararesults0.AddTabFieldDelimiter = 0
pararesults0.MergeConsecutiveDelimiters = 0

# get animation scene
animationScene1 = GetAnimationScene()

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# Properties modified on pararesults0
pararesults0.FieldDelimiterCharacters = ' '

# Create a new 'SpreadSheet View'
spreadSheetView1 = CreateView('SpreadSheetView')
spreadSheetView1.UseCache = 0
spreadSheetView1.ViewSize = [400, 400]
spreadSheetView1.CellFontSize = 9
spreadSheetView1.HeaderFontSize = 9
spreadSheetView1.SelectionOnly = 0
spreadSheetView1.GenerateCellConnectivity = 0
spreadSheetView1.ColumnToSort = ''
spreadSheetView1.InvertOrder = 0
spreadSheetView1.BlockSize = 1024
spreadSheetView1.HiddenColumnLabels = []
spreadSheetView1.FieldAssociation = 'Point Data'

# show data in view
pararesults0Display = Show(pararesults0, spreadSheetView1)

# trace defaults for the display properties.
pararesults0Display.CompositeDataSetIndex = [0]

# get layout
layout1 = GetLayoutByName("Layout #1")

# add view to a layout so it's visible in UI
AssignViewToLayout(view=spreadSheetView1, layout=layout1, hint=0)

# create a new 'Table To Points'
tableToPoints1 = TableToPoints(Input=pararesults0)
tableToPoints1.XColumn = 'CLV3Sig1'
tableToPoints1.YColumn = 'CLV3Sig1'
tableToPoints1.ZColumn = 'CLV3Sig1'
tableToPoints1.a2DPoints = 0
tableToPoints1.KeepAllDataArrays = 0

# Properties modified on tableToPoints1
tableToPoints1.XColumn = 'x'
tableToPoints1.YColumn = 'y'
tableToPoints1.ZColumn = 'z'

# show data in view
tableToPoints1Display = Show(tableToPoints1, spreadSheetView1)

# trace defaults for the display properties.
tableToPoints1Display.CompositeDataSetIndex = [0]

# hide data in view
Hide(pararesults0, spreadSheetView1)

# update the view to ensure updated data information
spreadSheetView1.Update()

# find view
renderView1 = FindViewOrCreate('RenderView1', viewtype='RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [547, 792]

# set active view
SetActiveView(renderView1)

# set active source
SetActiveSource(tableToPoints1)

# show data in view
tableToPoints1Display_1 = Show(tableToPoints1, renderView1)

# trace defaults for the display properties.
tableToPoints1Display_1.Representation = 'Surface'
tableToPoints1Display_1.AmbientColor = [1.0, 1.0, 1.0]
tableToPoints1Display_1.ColorArrayName = [None, '']
tableToPoints1Display_1.DiffuseColor = [1.0, 1.0, 1.0]
tableToPoints1Display_1.LookupTable = None
tableToPoints1Display_1.MapScalars = 1
tableToPoints1Display_1.MultiComponentsMapping = 0
tableToPoints1Display_1.InterpolateScalarsBeforeMapping = 1
tableToPoints1Display_1.Opacity = 1.0
tableToPoints1Display_1.PointSize = 2.0
tableToPoints1Display_1.LineWidth = 1.0
tableToPoints1Display_1.RenderLinesAsTubes = 0
tableToPoints1Display_1.RenderPointsAsSpheres = 0
tableToPoints1Display_1.Interpolation = 'Gouraud'
tableToPoints1Display_1.Specular = 0.0
tableToPoints1Display_1.SpecularColor = [1.0, 1.0, 1.0]
tableToPoints1Display_1.SpecularPower = 100.0
tableToPoints1Display_1.Luminosity = 0.0
tableToPoints1Display_1.Ambient = 0.0
tableToPoints1Display_1.Diffuse = 1.0
tableToPoints1Display_1.EdgeColor = [0.0, 0.0, 0.5]
tableToPoints1Display_1.BackfaceRepresentation = 'Follow Frontface'
tableToPoints1Display_1.BackfaceAmbientColor = [1.0, 1.0, 1.0]
tableToPoints1Display_1.BackfaceDiffuseColor = [1.0, 1.0, 1.0]
tableToPoints1Display_1.BackfaceOpacity = 1.0
tableToPoints1Display_1.Position = [0.0, 0.0, 0.0]
tableToPoints1Display_1.Scale = [1.0, 1.0, 1.0]
tableToPoints1Display_1.Orientation = [0.0, 0.0, 0.0]
tableToPoints1Display_1.Origin = [0.0, 0.0, 0.0]
tableToPoints1Display_1.Pickable = 1
tableToPoints1Display_1.Texture = None
tableToPoints1Display_1.Triangulate = 0
tableToPoints1Display_1.UseShaderReplacements = 0
tableToPoints1Display_1.ShaderReplacements = ''
tableToPoints1Display_1.NonlinearSubdivisionLevel = 1
tableToPoints1Display_1.UseDataPartitions = 0
tableToPoints1Display_1.OSPRayUseScaleArray = 0
tableToPoints1Display_1.OSPRayScaleArray = 'CLV3Sig1'
tableToPoints1Display_1.OSPRayScaleFunction = 'PiecewiseFunction'
tableToPoints1Display_1.OSPRayMaterial = 'None'
tableToPoints1Display_1.Orient = 0
tableToPoints1Display_1.OrientationMode = 'Direction'
tableToPoints1Display_1.SelectOrientationVectors = 'CLV3Sig1'
tableToPoints1Display_1.Scaling = 0
tableToPoints1Display_1.ScaleMode = 'No Data Scaling Off'
tableToPoints1Display_1.ScaleFactor = 1.82416
tableToPoints1Display_1.SelectScaleArray = 'CLV3Sig1'
tableToPoints1Display_1.GlyphType = 'Arrow'
tableToPoints1Display_1.UseGlyphTable = 0
tableToPoints1Display_1.GlyphTableIndexArray = 'CLV3Sig1'
tableToPoints1Display_1.UseCompositeGlyphTable = 0
tableToPoints1Display_1.UseGlyphCullingAndLOD = 0
tableToPoints1Display_1.LODValues = []
tableToPoints1Display_1.ColorByLODIndex = 0
tableToPoints1Display_1.GaussianRadius = 0.091208
tableToPoints1Display_1.ShaderPreset = 'Sphere'
tableToPoints1Display_1.CustomTriangleScale = 3
tableToPoints1Display_1.CustomShader = """ // This custom shader code define a gaussian blur
 // Please take a look into vtkSMPointGaussianRepresentation.cxx
 // for other custom shader examples
 //VTK::Color::Impl
   float dist2 = dot(offsetVCVSOutput.xy,offsetVCVSOutput.xy);
   float gaussian = exp(-0.5*dist2);
   opacity = opacity*gaussian;
"""
tableToPoints1Display_1.Emissive = 0
tableToPoints1Display_1.ScaleByArray = 0
tableToPoints1Display_1.SetScaleArray = ['POINTS', 'CLV3Sig1']
tableToPoints1Display_1.ScaleArrayComponent = ''
tableToPoints1Display_1.UseScaleFunction = 1
tableToPoints1Display_1.ScaleTransferFunction = 'PiecewiseFunction'
tableToPoints1Display_1.OpacityByArray = 0
tableToPoints1Display_1.OpacityArray = ['POINTS', 'CLV3Sig1']
tableToPoints1Display_1.OpacityArrayComponent = ''
tableToPoints1Display_1.OpacityTransferFunction = 'PiecewiseFunction'
tableToPoints1Display_1.DataAxesGrid = 'GridAxesRepresentation'
tableToPoints1Display_1.SelectionCellLabelBold = 0
tableToPoints1Display_1.SelectionCellLabelColor = [0.0, 1.0, 0.0]
tableToPoints1Display_1.SelectionCellLabelFontFamily = 'Arial'
tableToPoints1Display_1.SelectionCellLabelFontFile = ''
tableToPoints1Display_1.SelectionCellLabelFontSize = 18
tableToPoints1Display_1.SelectionCellLabelItalic = 0
tableToPoints1Display_1.SelectionCellLabelJustification = 'Left'
tableToPoints1Display_1.SelectionCellLabelOpacity = 1.0
tableToPoints1Display_1.SelectionCellLabelShadow = 0
tableToPoints1Display_1.SelectionPointLabelBold = 0
tableToPoints1Display_1.SelectionPointLabelColor = [1.0, 1.0, 0.0]
tableToPoints1Display_1.SelectionPointLabelFontFamily = 'Arial'
tableToPoints1Display_1.SelectionPointLabelFontFile = ''
tableToPoints1Display_1.SelectionPointLabelFontSize = 18
tableToPoints1Display_1.SelectionPointLabelItalic = 0
tableToPoints1Display_1.SelectionPointLabelJustification = 'Left'
tableToPoints1Display_1.SelectionPointLabelOpacity = 1.0
tableToPoints1Display_1.SelectionPointLabelShadow = 0
tableToPoints1Display_1.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
tableToPoints1Display_1.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
tableToPoints1Display_1.OSPRayScaleFunction.UseLogScale = 0

# init the 'Arrow' selected for 'GlyphType'
tableToPoints1Display_1.GlyphType.TipResolution = 6
tableToPoints1Display_1.GlyphType.TipRadius = 0.1
tableToPoints1Display_1.GlyphType.TipLength = 0.35
tableToPoints1Display_1.GlyphType.ShaftResolution = 6
tableToPoints1Display_1.GlyphType.ShaftRadius = 0.03
tableToPoints1Display_1.GlyphType.Invert = 0

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
tableToPoints1Display_1.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]
tableToPoints1Display_1.ScaleTransferFunction.UseLogScale = 0

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
tableToPoints1Display_1.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]
tableToPoints1Display_1.OpacityTransferFunction.UseLogScale = 0

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
tableToPoints1Display_1.DataAxesGrid.XTitle = 'X Axis'
tableToPoints1Display_1.DataAxesGrid.YTitle = 'Y Axis'
tableToPoints1Display_1.DataAxesGrid.ZTitle = 'Z Axis'
tableToPoints1Display_1.DataAxesGrid.XTitleColor = [1.0, 1.0, 1.0]
tableToPoints1Display_1.DataAxesGrid.XTitleFontFamily = 'Arial'
tableToPoints1Display_1.DataAxesGrid.XTitleFontFile = ''
tableToPoints1Display_1.DataAxesGrid.XTitleBold = 0
tableToPoints1Display_1.DataAxesGrid.XTitleItalic = 0
tableToPoints1Display_1.DataAxesGrid.XTitleFontSize = 12
tableToPoints1Display_1.DataAxesGrid.XTitleShadow = 0
tableToPoints1Display_1.DataAxesGrid.XTitleOpacity = 1.0
tableToPoints1Display_1.DataAxesGrid.YTitleColor = [1.0, 1.0, 1.0]
tableToPoints1Display_1.DataAxesGrid.YTitleFontFamily = 'Arial'
tableToPoints1Display_1.DataAxesGrid.YTitleFontFile = ''
tableToPoints1Display_1.DataAxesGrid.YTitleBold = 0
tableToPoints1Display_1.DataAxesGrid.YTitleItalic = 0
tableToPoints1Display_1.DataAxesGrid.YTitleFontSize = 12
tableToPoints1Display_1.DataAxesGrid.YTitleShadow = 0
tableToPoints1Display_1.DataAxesGrid.YTitleOpacity = 1.0
tableToPoints1Display_1.DataAxesGrid.ZTitleColor = [1.0, 1.0, 1.0]
tableToPoints1Display_1.DataAxesGrid.ZTitleFontFamily = 'Arial'
tableToPoints1Display_1.DataAxesGrid.ZTitleFontFile = ''
tableToPoints1Display_1.DataAxesGrid.ZTitleBold = 0
tableToPoints1Display_1.DataAxesGrid.ZTitleItalic = 0
tableToPoints1Display_1.DataAxesGrid.ZTitleFontSize = 12
tableToPoints1Display_1.DataAxesGrid.ZTitleShadow = 0
tableToPoints1Display_1.DataAxesGrid.ZTitleOpacity = 1.0
tableToPoints1Display_1.DataAxesGrid.FacesToRender = 63
tableToPoints1Display_1.DataAxesGrid.CullBackface = 0
tableToPoints1Display_1.DataAxesGrid.CullFrontface = 1
tableToPoints1Display_1.DataAxesGrid.GridColor = [1.0, 1.0, 1.0]
tableToPoints1Display_1.DataAxesGrid.ShowGrid = 0
tableToPoints1Display_1.DataAxesGrid.ShowEdges = 1
tableToPoints1Display_1.DataAxesGrid.ShowTicks = 1
tableToPoints1Display_1.DataAxesGrid.LabelUniqueEdgesOnly = 1
tableToPoints1Display_1.DataAxesGrid.AxesToLabel = 63
tableToPoints1Display_1.DataAxesGrid.XLabelColor = [1.0, 1.0, 1.0]
tableToPoints1Display_1.DataAxesGrid.XLabelFontFamily = 'Arial'
tableToPoints1Display_1.DataAxesGrid.XLabelFontFile = ''
tableToPoints1Display_1.DataAxesGrid.XLabelBold = 0
tableToPoints1Display_1.DataAxesGrid.XLabelItalic = 0
tableToPoints1Display_1.DataAxesGrid.XLabelFontSize = 12
tableToPoints1Display_1.DataAxesGrid.XLabelShadow = 0
tableToPoints1Display_1.DataAxesGrid.XLabelOpacity = 1.0
tableToPoints1Display_1.DataAxesGrid.YLabelColor = [1.0, 1.0, 1.0]
tableToPoints1Display_1.DataAxesGrid.YLabelFontFamily = 'Arial'
tableToPoints1Display_1.DataAxesGrid.YLabelFontFile = ''
tableToPoints1Display_1.DataAxesGrid.YLabelBold = 0
tableToPoints1Display_1.DataAxesGrid.YLabelItalic = 0
tableToPoints1Display_1.DataAxesGrid.YLabelFontSize = 12
tableToPoints1Display_1.DataAxesGrid.YLabelShadow = 0
tableToPoints1Display_1.DataAxesGrid.YLabelOpacity = 1.0
tableToPoints1Display_1.DataAxesGrid.ZLabelColor = [1.0, 1.0, 1.0]
tableToPoints1Display_1.DataAxesGrid.ZLabelFontFamily = 'Arial'
tableToPoints1Display_1.DataAxesGrid.ZLabelFontFile = ''
tableToPoints1Display_1.DataAxesGrid.ZLabelBold = 0
tableToPoints1Display_1.DataAxesGrid.ZLabelItalic = 0
tableToPoints1Display_1.DataAxesGrid.ZLabelFontSize = 12
tableToPoints1Display_1.DataAxesGrid.ZLabelShadow = 0
tableToPoints1Display_1.DataAxesGrid.ZLabelOpacity = 1.0
tableToPoints1Display_1.DataAxesGrid.XAxisNotation = 'Mixed'
tableToPoints1Display_1.DataAxesGrid.XAxisPrecision = 2
tableToPoints1Display_1.DataAxesGrid.XAxisUseCustomLabels = 0
tableToPoints1Display_1.DataAxesGrid.XAxisLabels = []
tableToPoints1Display_1.DataAxesGrid.YAxisNotation = 'Mixed'
tableToPoints1Display_1.DataAxesGrid.YAxisPrecision = 2
tableToPoints1Display_1.DataAxesGrid.YAxisUseCustomLabels = 0
tableToPoints1Display_1.DataAxesGrid.YAxisLabels = []
tableToPoints1Display_1.DataAxesGrid.ZAxisNotation = 'Mixed'
tableToPoints1Display_1.DataAxesGrid.ZAxisPrecision = 2
tableToPoints1Display_1.DataAxesGrid.ZAxisUseCustomLabels = 0
tableToPoints1Display_1.DataAxesGrid.ZAxisLabels = []
tableToPoints1Display_1.DataAxesGrid.UseCustomBounds = 0
tableToPoints1Display_1.DataAxesGrid.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
tableToPoints1Display_1.PolarAxes.Visibility = 0
tableToPoints1Display_1.PolarAxes.Translation = [0.0, 0.0, 0.0]
tableToPoints1Display_1.PolarAxes.Scale = [1.0, 1.0, 1.0]
tableToPoints1Display_1.PolarAxes.Orientation = [0.0, 0.0, 0.0]
tableToPoints1Display_1.PolarAxes.EnableCustomBounds = [0, 0, 0]
tableToPoints1Display_1.PolarAxes.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
tableToPoints1Display_1.PolarAxes.EnableCustomRange = 0
tableToPoints1Display_1.PolarAxes.CustomRange = [0.0, 1.0]
tableToPoints1Display_1.PolarAxes.PolarAxisVisibility = 1
tableToPoints1Display_1.PolarAxes.RadialAxesVisibility = 1
tableToPoints1Display_1.PolarAxes.DrawRadialGridlines = 1
tableToPoints1Display_1.PolarAxes.PolarArcsVisibility = 1
tableToPoints1Display_1.PolarAxes.DrawPolarArcsGridlines = 1
tableToPoints1Display_1.PolarAxes.NumberOfRadialAxes = 0
tableToPoints1Display_1.PolarAxes.AutoSubdividePolarAxis = 1
tableToPoints1Display_1.PolarAxes.NumberOfPolarAxis = 0
tableToPoints1Display_1.PolarAxes.MinimumRadius = 0.0
tableToPoints1Display_1.PolarAxes.MinimumAngle = 0.0
tableToPoints1Display_1.PolarAxes.MaximumAngle = 90.0
tableToPoints1Display_1.PolarAxes.RadialAxesOriginToPolarAxis = 1
tableToPoints1Display_1.PolarAxes.Ratio = 1.0
tableToPoints1Display_1.PolarAxes.PolarAxisColor = [1.0, 1.0, 1.0]
tableToPoints1Display_1.PolarAxes.PolarArcsColor = [1.0, 1.0, 1.0]
tableToPoints1Display_1.PolarAxes.LastRadialAxisColor = [1.0, 1.0, 1.0]
tableToPoints1Display_1.PolarAxes.SecondaryPolarArcsColor = [1.0, 1.0, 1.0]
tableToPoints1Display_1.PolarAxes.SecondaryRadialAxesColor = [1.0, 1.0, 1.0]
tableToPoints1Display_1.PolarAxes.PolarAxisTitleVisibility = 1
tableToPoints1Display_1.PolarAxes.PolarAxisTitle = 'Radial Distance'
tableToPoints1Display_1.PolarAxes.PolarAxisTitleLocation = 'Bottom'
tableToPoints1Display_1.PolarAxes.PolarLabelVisibility = 1
tableToPoints1Display_1.PolarAxes.PolarLabelFormat = '%-#6.3g'
tableToPoints1Display_1.PolarAxes.PolarLabelExponentLocation = 'Labels'
tableToPoints1Display_1.PolarAxes.RadialLabelVisibility = 1
tableToPoints1Display_1.PolarAxes.RadialLabelFormat = '%-#3.1f'
tableToPoints1Display_1.PolarAxes.RadialLabelLocation = 'Bottom'
tableToPoints1Display_1.PolarAxes.RadialUnitsVisibility = 1
tableToPoints1Display_1.PolarAxes.ScreenSize = 10.0
tableToPoints1Display_1.PolarAxes.PolarAxisTitleColor = [1.0, 1.0, 1.0]
tableToPoints1Display_1.PolarAxes.PolarAxisTitleOpacity = 1.0
tableToPoints1Display_1.PolarAxes.PolarAxisTitleFontFamily = 'Arial'
tableToPoints1Display_1.PolarAxes.PolarAxisTitleFontFile = ''
tableToPoints1Display_1.PolarAxes.PolarAxisTitleBold = 0
tableToPoints1Display_1.PolarAxes.PolarAxisTitleItalic = 0
tableToPoints1Display_1.PolarAxes.PolarAxisTitleShadow = 0
tableToPoints1Display_1.PolarAxes.PolarAxisTitleFontSize = 12
tableToPoints1Display_1.PolarAxes.PolarAxisLabelColor = [1.0, 1.0, 1.0]
tableToPoints1Display_1.PolarAxes.PolarAxisLabelOpacity = 1.0
tableToPoints1Display_1.PolarAxes.PolarAxisLabelFontFamily = 'Arial'
tableToPoints1Display_1.PolarAxes.PolarAxisLabelFontFile = ''
tableToPoints1Display_1.PolarAxes.PolarAxisLabelBold = 0
tableToPoints1Display_1.PolarAxes.PolarAxisLabelItalic = 0
tableToPoints1Display_1.PolarAxes.PolarAxisLabelShadow = 0
tableToPoints1Display_1.PolarAxes.PolarAxisLabelFontSize = 12
tableToPoints1Display_1.PolarAxes.LastRadialAxisTextColor = [1.0, 1.0, 1.0]
tableToPoints1Display_1.PolarAxes.LastRadialAxisTextOpacity = 1.0
tableToPoints1Display_1.PolarAxes.LastRadialAxisTextFontFamily = 'Arial'
tableToPoints1Display_1.PolarAxes.LastRadialAxisTextFontFile = ''
tableToPoints1Display_1.PolarAxes.LastRadialAxisTextBold = 0
tableToPoints1Display_1.PolarAxes.LastRadialAxisTextItalic = 0
tableToPoints1Display_1.PolarAxes.LastRadialAxisTextShadow = 0
tableToPoints1Display_1.PolarAxes.LastRadialAxisTextFontSize = 12
tableToPoints1Display_1.PolarAxes.SecondaryRadialAxesTextColor = [1.0, 1.0, 1.0]
tableToPoints1Display_1.PolarAxes.SecondaryRadialAxesTextOpacity = 1.0
tableToPoints1Display_1.PolarAxes.SecondaryRadialAxesTextFontFamily = 'Arial'
tableToPoints1Display_1.PolarAxes.SecondaryRadialAxesTextFontFile = ''
tableToPoints1Display_1.PolarAxes.SecondaryRadialAxesTextBold = 0
tableToPoints1Display_1.PolarAxes.SecondaryRadialAxesTextItalic = 0
tableToPoints1Display_1.PolarAxes.SecondaryRadialAxesTextShadow = 0
tableToPoints1Display_1.PolarAxes.SecondaryRadialAxesTextFontSize = 12
tableToPoints1Display_1.PolarAxes.EnableDistanceLOD = 1
tableToPoints1Display_1.PolarAxes.DistanceLODThreshold = 0.7
tableToPoints1Display_1.PolarAxes.EnableViewAngleLOD = 1
tableToPoints1Display_1.PolarAxes.ViewAngleLODThreshold = 0.7
tableToPoints1Display_1.PolarAxes.SmallestVisiblePolarAngle = 0.5
tableToPoints1Display_1.PolarAxes.PolarTicksVisibility = 1
tableToPoints1Display_1.PolarAxes.ArcTicksOriginToPolarAxis = 1
tableToPoints1Display_1.PolarAxes.TickLocation = 'Both'
tableToPoints1Display_1.PolarAxes.AxisTickVisibility = 1
tableToPoints1Display_1.PolarAxes.AxisMinorTickVisibility = 0
tableToPoints1Display_1.PolarAxes.ArcTickVisibility = 1
tableToPoints1Display_1.PolarAxes.ArcMinorTickVisibility = 0
tableToPoints1Display_1.PolarAxes.DeltaAngleMajor = 10.0
tableToPoints1Display_1.PolarAxes.DeltaAngleMinor = 5.0
tableToPoints1Display_1.PolarAxes.PolarAxisMajorTickSize = 0.0
tableToPoints1Display_1.PolarAxes.PolarAxisTickRatioSize = 0.3
tableToPoints1Display_1.PolarAxes.PolarAxisMajorTickThickness = 1.0
tableToPoints1Display_1.PolarAxes.PolarAxisTickRatioThickness = 0.5
tableToPoints1Display_1.PolarAxes.LastRadialAxisMajorTickSize = 0.0
tableToPoints1Display_1.PolarAxes.LastRadialAxisTickRatioSize = 0.3
tableToPoints1Display_1.PolarAxes.LastRadialAxisMajorTickThickness = 1.0
tableToPoints1Display_1.PolarAxes.LastRadialAxisTickRatioThickness = 0.5
tableToPoints1Display_1.PolarAxes.ArcMajorTickSize = 0.0
tableToPoints1Display_1.PolarAxes.ArcTickRatioSize = 0.3
tableToPoints1Display_1.PolarAxes.ArcMajorTickThickness = 1.0
tableToPoints1Display_1.PolarAxes.ArcTickRatioThickness = 0.5
tableToPoints1Display_1.PolarAxes.Use2DMode = 0
tableToPoints1Display_1.PolarAxes.UseLogAxis = 0

# reset view to fit data
renderView1.ResetCamera()

# create a new 'Clip'
clip1 = Clip(Input=tableToPoints1)
clip1.ClipType = 'Plane'
clip1.Scalars = ['POINTS', 'CLV3Sig1']
clip1.Value = 0.0
clip1.Invert = 1
clip1.Crinkleclip = 0
clip1.Exact = 0

# init the 'Plane' selected for 'ClipType'
clip1.ClipType.Origin = [-0.031060000000000088, -0.056819999999999204, 4.889417]
clip1.ClipType.Normal = [1.0, 0.0, 0.0]
clip1.ClipType.Offset = 0.0

# Properties modified on clip1.ClipType
clip1.ClipType.Origin = [0.0, 0.0, 0.0]

# show data in view
clip1Display = Show(clip1, renderView1)

# trace defaults for the display properties.
clip1Display.Representation = 'Surface'
clip1Display.AmbientColor = [1.0, 1.0, 1.0]
clip1Display.ColorArrayName = [None, '']
clip1Display.DiffuseColor = [1.0, 1.0, 1.0]
clip1Display.LookupTable = None
clip1Display.MapScalars = 1
clip1Display.MultiComponentsMapping = 0
clip1Display.InterpolateScalarsBeforeMapping = 1
clip1Display.Opacity = 1.0
clip1Display.PointSize = 2.0
clip1Display.LineWidth = 1.0
clip1Display.RenderLinesAsTubes = 0
clip1Display.RenderPointsAsSpheres = 0
clip1Display.Interpolation = 'Gouraud'
clip1Display.Specular = 0.0
clip1Display.SpecularColor = [1.0, 1.0, 1.0]
clip1Display.SpecularPower = 100.0
clip1Display.Luminosity = 0.0
clip1Display.Ambient = 0.0
clip1Display.Diffuse = 1.0
clip1Display.EdgeColor = [0.0, 0.0, 0.5]
clip1Display.BackfaceRepresentation = 'Follow Frontface'
clip1Display.BackfaceAmbientColor = [1.0, 1.0, 1.0]
clip1Display.BackfaceDiffuseColor = [1.0, 1.0, 1.0]
clip1Display.BackfaceOpacity = 1.0
clip1Display.Position = [0.0, 0.0, 0.0]
clip1Display.Scale = [1.0, 1.0, 1.0]
clip1Display.Orientation = [0.0, 0.0, 0.0]
clip1Display.Origin = [0.0, 0.0, 0.0]
clip1Display.Pickable = 1
clip1Display.Texture = None
clip1Display.Triangulate = 0
clip1Display.UseShaderReplacements = 0
clip1Display.ShaderReplacements = ''
clip1Display.NonlinearSubdivisionLevel = 1
clip1Display.UseDataPartitions = 0
clip1Display.OSPRayUseScaleArray = 0
clip1Display.OSPRayScaleArray = 'CLV3Sig1'
clip1Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip1Display.OSPRayMaterial = 'None'
clip1Display.Orient = 0
clip1Display.OrientationMode = 'Direction'
clip1Display.SelectOrientationVectors = 'CLV3Sig1'
clip1Display.Scaling = 0
clip1Display.ScaleMode = 'No Data Scaling Off'
clip1Display.ScaleFactor = 1.8161110000000003
clip1Display.SelectScaleArray = 'CLV3Sig1'
clip1Display.GlyphType = 'Arrow'
clip1Display.UseGlyphTable = 0
clip1Display.GlyphTableIndexArray = 'CLV3Sig1'
clip1Display.UseCompositeGlyphTable = 0
clip1Display.UseGlyphCullingAndLOD = 0
clip1Display.LODValues = []
clip1Display.ColorByLODIndex = 0
clip1Display.GaussianRadius = 0.09080555
clip1Display.ShaderPreset = 'Sphere'
clip1Display.CustomTriangleScale = 3
clip1Display.CustomShader = """ // This custom shader code define a gaussian blur
 // Please take a look into vtkSMPointGaussianRepresentation.cxx
 // for other custom shader examples
 //VTK::Color::Impl
   float dist2 = dot(offsetVCVSOutput.xy,offsetVCVSOutput.xy);
   float gaussian = exp(-0.5*dist2);
   opacity = opacity*gaussian;
"""
clip1Display.Emissive = 0
clip1Display.ScaleByArray = 0
clip1Display.SetScaleArray = ['POINTS', 'CLV3Sig1']
clip1Display.ScaleArrayComponent = ''
clip1Display.UseScaleFunction = 1
clip1Display.ScaleTransferFunction = 'PiecewiseFunction'
clip1Display.OpacityByArray = 0
clip1Display.OpacityArray = ['POINTS', 'CLV3Sig1']
clip1Display.OpacityArrayComponent = ''
clip1Display.OpacityTransferFunction = 'PiecewiseFunction'
clip1Display.DataAxesGrid = 'GridAxesRepresentation'
clip1Display.SelectionCellLabelBold = 0
clip1Display.SelectionCellLabelColor = [0.0, 1.0, 0.0]
clip1Display.SelectionCellLabelFontFamily = 'Arial'
clip1Display.SelectionCellLabelFontFile = ''
clip1Display.SelectionCellLabelFontSize = 18
clip1Display.SelectionCellLabelItalic = 0
clip1Display.SelectionCellLabelJustification = 'Left'
clip1Display.SelectionCellLabelOpacity = 1.0
clip1Display.SelectionCellLabelShadow = 0
clip1Display.SelectionPointLabelBold = 0
clip1Display.SelectionPointLabelColor = [1.0, 1.0, 0.0]
clip1Display.SelectionPointLabelFontFamily = 'Arial'
clip1Display.SelectionPointLabelFontFile = ''
clip1Display.SelectionPointLabelFontSize = 18
clip1Display.SelectionPointLabelItalic = 0
clip1Display.SelectionPointLabelJustification = 'Left'
clip1Display.SelectionPointLabelOpacity = 1.0
clip1Display.SelectionPointLabelShadow = 0
clip1Display.PolarAxes = 'PolarAxesRepresentation'
clip1Display.ScalarOpacityFunction = None
clip1Display.ScalarOpacityUnitDistance = 2.5332876999781373
clip1Display.ExtractedBlockIndex = 0
clip1Display.SelectMapper = 'Projected tetra'
clip1Display.SamplingDimensions = [128, 128, 128]
clip1Display.UseFloatingPointFrameBuffer = 1

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
clip1Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
clip1Display.OSPRayScaleFunction.UseLogScale = 0

# init the 'Arrow' selected for 'GlyphType'
clip1Display.GlyphType.TipResolution = 6
clip1Display.GlyphType.TipRadius = 0.1
clip1Display.GlyphType.TipLength = 0.35
clip1Display.GlyphType.ShaftResolution = 6
clip1Display.GlyphType.ShaftRadius = 0.03
clip1Display.GlyphType.Invert = 0

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]
clip1Display.ScaleTransferFunction.UseLogScale = 0

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]
clip1Display.OpacityTransferFunction.UseLogScale = 0

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
clip1Display.DataAxesGrid.XTitle = 'X Axis'
clip1Display.DataAxesGrid.YTitle = 'Y Axis'
clip1Display.DataAxesGrid.ZTitle = 'Z Axis'
clip1Display.DataAxesGrid.XTitleColor = [1.0, 1.0, 1.0]
clip1Display.DataAxesGrid.XTitleFontFamily = 'Arial'
clip1Display.DataAxesGrid.XTitleFontFile = ''
clip1Display.DataAxesGrid.XTitleBold = 0
clip1Display.DataAxesGrid.XTitleItalic = 0
clip1Display.DataAxesGrid.XTitleFontSize = 12
clip1Display.DataAxesGrid.XTitleShadow = 0
clip1Display.DataAxesGrid.XTitleOpacity = 1.0
clip1Display.DataAxesGrid.YTitleColor = [1.0, 1.0, 1.0]
clip1Display.DataAxesGrid.YTitleFontFamily = 'Arial'
clip1Display.DataAxesGrid.YTitleFontFile = ''
clip1Display.DataAxesGrid.YTitleBold = 0
clip1Display.DataAxesGrid.YTitleItalic = 0
clip1Display.DataAxesGrid.YTitleFontSize = 12
clip1Display.DataAxesGrid.YTitleShadow = 0
clip1Display.DataAxesGrid.YTitleOpacity = 1.0
clip1Display.DataAxesGrid.ZTitleColor = [1.0, 1.0, 1.0]
clip1Display.DataAxesGrid.ZTitleFontFamily = 'Arial'
clip1Display.DataAxesGrid.ZTitleFontFile = ''
clip1Display.DataAxesGrid.ZTitleBold = 0
clip1Display.DataAxesGrid.ZTitleItalic = 0
clip1Display.DataAxesGrid.ZTitleFontSize = 12
clip1Display.DataAxesGrid.ZTitleShadow = 0
clip1Display.DataAxesGrid.ZTitleOpacity = 1.0
clip1Display.DataAxesGrid.FacesToRender = 63
clip1Display.DataAxesGrid.CullBackface = 0
clip1Display.DataAxesGrid.CullFrontface = 1
clip1Display.DataAxesGrid.GridColor = [1.0, 1.0, 1.0]
clip1Display.DataAxesGrid.ShowGrid = 0
clip1Display.DataAxesGrid.ShowEdges = 1
clip1Display.DataAxesGrid.ShowTicks = 1
clip1Display.DataAxesGrid.LabelUniqueEdgesOnly = 1
clip1Display.DataAxesGrid.AxesToLabel = 63
clip1Display.DataAxesGrid.XLabelColor = [1.0, 1.0, 1.0]
clip1Display.DataAxesGrid.XLabelFontFamily = 'Arial'
clip1Display.DataAxesGrid.XLabelFontFile = ''
clip1Display.DataAxesGrid.XLabelBold = 0
clip1Display.DataAxesGrid.XLabelItalic = 0
clip1Display.DataAxesGrid.XLabelFontSize = 12
clip1Display.DataAxesGrid.XLabelShadow = 0
clip1Display.DataAxesGrid.XLabelOpacity = 1.0
clip1Display.DataAxesGrid.YLabelColor = [1.0, 1.0, 1.0]
clip1Display.DataAxesGrid.YLabelFontFamily = 'Arial'
clip1Display.DataAxesGrid.YLabelFontFile = ''
clip1Display.DataAxesGrid.YLabelBold = 0
clip1Display.DataAxesGrid.YLabelItalic = 0
clip1Display.DataAxesGrid.YLabelFontSize = 12
clip1Display.DataAxesGrid.YLabelShadow = 0
clip1Display.DataAxesGrid.YLabelOpacity = 1.0
clip1Display.DataAxesGrid.ZLabelColor = [1.0, 1.0, 1.0]
clip1Display.DataAxesGrid.ZLabelFontFamily = 'Arial'
clip1Display.DataAxesGrid.ZLabelFontFile = ''
clip1Display.DataAxesGrid.ZLabelBold = 0
clip1Display.DataAxesGrid.ZLabelItalic = 0
clip1Display.DataAxesGrid.ZLabelFontSize = 12
clip1Display.DataAxesGrid.ZLabelShadow = 0
clip1Display.DataAxesGrid.ZLabelOpacity = 1.0
clip1Display.DataAxesGrid.XAxisNotation = 'Mixed'
clip1Display.DataAxesGrid.XAxisPrecision = 2
clip1Display.DataAxesGrid.XAxisUseCustomLabels = 0
clip1Display.DataAxesGrid.XAxisLabels = []
clip1Display.DataAxesGrid.YAxisNotation = 'Mixed'
clip1Display.DataAxesGrid.YAxisPrecision = 2
clip1Display.DataAxesGrid.YAxisUseCustomLabels = 0
clip1Display.DataAxesGrid.YAxisLabels = []
clip1Display.DataAxesGrid.ZAxisNotation = 'Mixed'
clip1Display.DataAxesGrid.ZAxisPrecision = 2
clip1Display.DataAxesGrid.ZAxisUseCustomLabels = 0
clip1Display.DataAxesGrid.ZAxisLabels = []
clip1Display.DataAxesGrid.UseCustomBounds = 0
clip1Display.DataAxesGrid.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
clip1Display.PolarAxes.Visibility = 0
clip1Display.PolarAxes.Translation = [0.0, 0.0, 0.0]
clip1Display.PolarAxes.Scale = [1.0, 1.0, 1.0]
clip1Display.PolarAxes.Orientation = [0.0, 0.0, 0.0]
clip1Display.PolarAxes.EnableCustomBounds = [0, 0, 0]
clip1Display.PolarAxes.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
clip1Display.PolarAxes.EnableCustomRange = 0
clip1Display.PolarAxes.CustomRange = [0.0, 1.0]
clip1Display.PolarAxes.PolarAxisVisibility = 1
clip1Display.PolarAxes.RadialAxesVisibility = 1
clip1Display.PolarAxes.DrawRadialGridlines = 1
clip1Display.PolarAxes.PolarArcsVisibility = 1
clip1Display.PolarAxes.DrawPolarArcsGridlines = 1
clip1Display.PolarAxes.NumberOfRadialAxes = 0
clip1Display.PolarAxes.AutoSubdividePolarAxis = 1
clip1Display.PolarAxes.NumberOfPolarAxis = 0
clip1Display.PolarAxes.MinimumRadius = 0.0
clip1Display.PolarAxes.MinimumAngle = 0.0
clip1Display.PolarAxes.MaximumAngle = 90.0
clip1Display.PolarAxes.RadialAxesOriginToPolarAxis = 1
clip1Display.PolarAxes.Ratio = 1.0
clip1Display.PolarAxes.PolarAxisColor = [1.0, 1.0, 1.0]
clip1Display.PolarAxes.PolarArcsColor = [1.0, 1.0, 1.0]
clip1Display.PolarAxes.LastRadialAxisColor = [1.0, 1.0, 1.0]
clip1Display.PolarAxes.SecondaryPolarArcsColor = [1.0, 1.0, 1.0]
clip1Display.PolarAxes.SecondaryRadialAxesColor = [1.0, 1.0, 1.0]
clip1Display.PolarAxes.PolarAxisTitleVisibility = 1
clip1Display.PolarAxes.PolarAxisTitle = 'Radial Distance'
clip1Display.PolarAxes.PolarAxisTitleLocation = 'Bottom'
clip1Display.PolarAxes.PolarLabelVisibility = 1
clip1Display.PolarAxes.PolarLabelFormat = '%-#6.3g'
clip1Display.PolarAxes.PolarLabelExponentLocation = 'Labels'
clip1Display.PolarAxes.RadialLabelVisibility = 1
clip1Display.PolarAxes.RadialLabelFormat = '%-#3.1f'
clip1Display.PolarAxes.RadialLabelLocation = 'Bottom'
clip1Display.PolarAxes.RadialUnitsVisibility = 1
clip1Display.PolarAxes.ScreenSize = 10.0
clip1Display.PolarAxes.PolarAxisTitleColor = [1.0, 1.0, 1.0]
clip1Display.PolarAxes.PolarAxisTitleOpacity = 1.0
clip1Display.PolarAxes.PolarAxisTitleFontFamily = 'Arial'
clip1Display.PolarAxes.PolarAxisTitleFontFile = ''
clip1Display.PolarAxes.PolarAxisTitleBold = 0
clip1Display.PolarAxes.PolarAxisTitleItalic = 0
clip1Display.PolarAxes.PolarAxisTitleShadow = 0
clip1Display.PolarAxes.PolarAxisTitleFontSize = 12
clip1Display.PolarAxes.PolarAxisLabelColor = [1.0, 1.0, 1.0]
clip1Display.PolarAxes.PolarAxisLabelOpacity = 1.0
clip1Display.PolarAxes.PolarAxisLabelFontFamily = 'Arial'
clip1Display.PolarAxes.PolarAxisLabelFontFile = ''
clip1Display.PolarAxes.PolarAxisLabelBold = 0
clip1Display.PolarAxes.PolarAxisLabelItalic = 0
clip1Display.PolarAxes.PolarAxisLabelShadow = 0
clip1Display.PolarAxes.PolarAxisLabelFontSize = 12
clip1Display.PolarAxes.LastRadialAxisTextColor = [1.0, 1.0, 1.0]
clip1Display.PolarAxes.LastRadialAxisTextOpacity = 1.0
clip1Display.PolarAxes.LastRadialAxisTextFontFamily = 'Arial'
clip1Display.PolarAxes.LastRadialAxisTextFontFile = ''
clip1Display.PolarAxes.LastRadialAxisTextBold = 0
clip1Display.PolarAxes.LastRadialAxisTextItalic = 0
clip1Display.PolarAxes.LastRadialAxisTextShadow = 0
clip1Display.PolarAxes.LastRadialAxisTextFontSize = 12
clip1Display.PolarAxes.SecondaryRadialAxesTextColor = [1.0, 1.0, 1.0]
clip1Display.PolarAxes.SecondaryRadialAxesTextOpacity = 1.0
clip1Display.PolarAxes.SecondaryRadialAxesTextFontFamily = 'Arial'
clip1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''
clip1Display.PolarAxes.SecondaryRadialAxesTextBold = 0
clip1Display.PolarAxes.SecondaryRadialAxesTextItalic = 0
clip1Display.PolarAxes.SecondaryRadialAxesTextShadow = 0
clip1Display.PolarAxes.SecondaryRadialAxesTextFontSize = 12
clip1Display.PolarAxes.EnableDistanceLOD = 1
clip1Display.PolarAxes.DistanceLODThreshold = 0.7
clip1Display.PolarAxes.EnableViewAngleLOD = 1
clip1Display.PolarAxes.ViewAngleLODThreshold = 0.7
clip1Display.PolarAxes.SmallestVisiblePolarAngle = 0.5
clip1Display.PolarAxes.PolarTicksVisibility = 1
clip1Display.PolarAxes.ArcTicksOriginToPolarAxis = 1
clip1Display.PolarAxes.TickLocation = 'Both'
clip1Display.PolarAxes.AxisTickVisibility = 1
clip1Display.PolarAxes.AxisMinorTickVisibility = 0
clip1Display.PolarAxes.ArcTickVisibility = 1
clip1Display.PolarAxes.ArcMinorTickVisibility = 0
clip1Display.PolarAxes.DeltaAngleMajor = 10.0
clip1Display.PolarAxes.DeltaAngleMinor = 5.0
clip1Display.PolarAxes.PolarAxisMajorTickSize = 0.0
clip1Display.PolarAxes.PolarAxisTickRatioSize = 0.3
clip1Display.PolarAxes.PolarAxisMajorTickThickness = 1.0
clip1Display.PolarAxes.PolarAxisTickRatioThickness = 0.5
clip1Display.PolarAxes.LastRadialAxisMajorTickSize = 0.0
clip1Display.PolarAxes.LastRadialAxisTickRatioSize = 0.3
clip1Display.PolarAxes.LastRadialAxisMajorTickThickness = 1.0
clip1Display.PolarAxes.LastRadialAxisTickRatioThickness = 0.5
clip1Display.PolarAxes.ArcMajorTickSize = 0.0
clip1Display.PolarAxes.ArcTickRatioSize = 0.3
clip1Display.PolarAxes.ArcMajorTickThickness = 1.0
clip1Display.PolarAxes.ArcTickRatioThickness = 0.5
clip1Display.PolarAxes.Use2DMode = 0
clip1Display.PolarAxes.UseLogAxis = 0

# hide data in view
Hide(tableToPoints1, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# reset view to fit data
renderView1.ResetCamera()

# change representation type
clip1Display.SetRepresentationType('Point Gaussian')

# Properties modified on clip1Display
clip1Display.GaussianRadius = 0.6

# set scalar coloring
ColorBy(clip1Display, ('POINTS', 'CLV3Sig1'))

# rescale color and/or opacity maps used to include current data range
clip1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
clip1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'CLV3Sig1'
cLV3Sig1LUT = GetColorTransferFunction('CLV3Sig1')
cLV3Sig1LUT.AutomaticRescaleRangeMode = "Grow and update on 'Apply'"
cLV3Sig1LUT.InterpretValuesAsCategories = 0
cLV3Sig1LUT.AnnotationsInitialized = 0
cLV3Sig1LUT.ShowCategoricalColorsinDataRangeOnly = 0
cLV3Sig1LUT.RescaleOnVisibilityChange = 0
cLV3Sig1LUT.EnableOpacityMapping = 0
cLV3Sig1LUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 5.878906683738906e-39, 0.865003, 0.865003, 0.865003, 1.1757813367477812e-38, 0.705882, 0.0156863, 0.14902]
cLV3Sig1LUT.UseLogScale = 0
cLV3Sig1LUT.ColorSpace = 'Diverging'
cLV3Sig1LUT.UseBelowRangeColor = 0
cLV3Sig1LUT.BelowRangeColor = [0.0, 0.0, 0.0]
cLV3Sig1LUT.UseAboveRangeColor = 0
cLV3Sig1LUT.AboveRangeColor = [0.5, 0.5, 0.5]
cLV3Sig1LUT.NanColor = [1.0, 1.0, 0.0]
cLV3Sig1LUT.NanOpacity = 1.0
cLV3Sig1LUT.Discretize = 1
cLV3Sig1LUT.NumberOfTableValues = 256
cLV3Sig1LUT.ScalarRangeInitialized = 1.0
cLV3Sig1LUT.HSVWrap = 0
cLV3Sig1LUT.VectorComponent = 0
cLV3Sig1LUT.VectorMode = 'Magnitude'
cLV3Sig1LUT.AllowDuplicateScalars = 1
cLV3Sig1LUT.Annotations = []
cLV3Sig1LUT.ActiveAnnotatedValues = []
cLV3Sig1LUT.IndexedColors = []
cLV3Sig1LUT.IndexedOpacities = []

# get opacity transfer function/opacity map for 'CLV3Sig1'
cLV3Sig1PWF = GetOpacityTransferFunction('CLV3Sig1')
cLV3Sig1PWF.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]
cLV3Sig1PWF.AllowDuplicateScalars = 1
cLV3Sig1PWF.UseLogScale = 0
cLV3Sig1PWF.ScalarRangeInitialized = 1

# Rescale transfer function
cLV3Sig1LUT.RescaleTransferFunction(0.0, 100.0)

# Rescale transfer function
cLV3Sig1PWF.RescaleTransferFunction(0.0, 100.0)

animationScene1.Play()

# set active source
SetActiveSource(tableToPoints1)

# current camera placement for renderView1
renderView1.CameraPosition = [90.3652515703269, -0.056819915771484375, 4.81049570441246]
renderView1.CameraFocalPoint = [-0.031060000000000088, -0.056819915771484375, 4.81049570441246]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraParallelScale = 23.817211739818177

# save screenshot
SaveScreenshot('/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Pics/CLV3AccScreenShot_wtresults021621DcTrWt.png', renderView1, SaveAllViews=0,
    ImageResolution=[1896, 3168],
    FontScaling='Scale fonts proportionally',
    SeparatorWidth=1,
    SeparatorColor=[0.937, 0.922, 0.906],
    OverrideColorPalette='',
    StereoMode='No change',
    TransparentBackground=0, 
    # PNG options
    CompressionLevel='5')

# set active source
SetActiveSource(clip1)

# set scalar coloring
ColorBy(clip1Display, ('POINTS', 'CLV3Sig2'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(cLV3Sig1LUT, renderView1)

# rescale color and/or opacity maps used to include current data range
clip1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
clip1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'CLV3Sig2'
cLV3Sig2LUT = GetColorTransferFunction('CLV3Sig2')
cLV3Sig2LUT.AutomaticRescaleRangeMode = "Grow and update on 'Apply'"
cLV3Sig2LUT.InterpretValuesAsCategories = 0
cLV3Sig2LUT.AnnotationsInitialized = 0
cLV3Sig2LUT.ShowCategoricalColorsinDataRangeOnly = 0
cLV3Sig2LUT.RescaleOnVisibilityChange = 0
cLV3Sig2LUT.EnableOpacityMapping = 0
cLV3Sig2LUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 150.575, 0.865003, 0.865003, 0.865003, 301.15, 0.705882, 0.0156863, 0.14902]
cLV3Sig2LUT.UseLogScale = 0
cLV3Sig2LUT.ColorSpace = 'Diverging'
cLV3Sig2LUT.UseBelowRangeColor = 0
cLV3Sig2LUT.BelowRangeColor = [0.0, 0.0, 0.0]
cLV3Sig2LUT.UseAboveRangeColor = 0
cLV3Sig2LUT.AboveRangeColor = [0.5, 0.5, 0.5]
cLV3Sig2LUT.NanColor = [1.0, 1.0, 0.0]
cLV3Sig2LUT.NanOpacity = 1.0
cLV3Sig2LUT.Discretize = 1
cLV3Sig2LUT.NumberOfTableValues = 256
cLV3Sig2LUT.ScalarRangeInitialized = 1.0
cLV3Sig2LUT.HSVWrap = 0
cLV3Sig2LUT.VectorComponent = 0
cLV3Sig2LUT.VectorMode = 'Magnitude'
cLV3Sig2LUT.AllowDuplicateScalars = 1
cLV3Sig2LUT.Annotations = []
cLV3Sig2LUT.ActiveAnnotatedValues = []
cLV3Sig2LUT.IndexedColors = []
cLV3Sig2LUT.IndexedOpacities = []

# get opacity transfer function/opacity map for 'CLV3Sig2'
cLV3Sig2PWF = GetOpacityTransferFunction('CLV3Sig2')
cLV3Sig2PWF.Points = [0.0, 0.0, 0.5, 0.0, 301.15, 1.0, 0.5, 0.0]
cLV3Sig2PWF.AllowDuplicateScalars = 1
cLV3Sig2PWF.UseLogScale = 0
cLV3Sig2PWF.ScalarRangeInitialized = 1

# Rescale transfer function
cLV3Sig2LUT.RescaleTransferFunction(0.0, 100.0)

# Rescale transfer function
cLV3Sig2PWF.RescaleTransferFunction(0.0, 100.0)

animationScene1.Play()

# set active source
SetActiveSource(tableToPoints1)

# current camera placement for renderView1
renderView1.CameraPosition = [90.3652515703269, -0.056819915771484375, 4.81049570441246]
renderView1.CameraFocalPoint = [-0.031060000000000088, -0.056819915771484375, 4.81049570441246]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraParallelScale = 23.817211739818177

# save screenshot
SaveScreenshot('/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Pics/CLV3Ac2ScreenShot_wtresults021621DcTrWt.png', renderView1, SaveAllViews=0,
    ImageResolution=[1896, 3168],
    FontScaling='Scale fonts proportionally',
    SeparatorWidth=1,
    SeparatorColor=[0.937, 0.922, 0.906],
    OverrideColorPalette='',
    StereoMode='No change',
    TransparentBackground=0, 
    # PNG options
    CompressionLevel='5')

# set active source
SetActiveSource(clip1)

# set scalar coloring
ColorBy(clip1Display, ('POINTS', 'CLV3_Peptide'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(cLV3Sig2LUT, renderView1)

# rescale color and/or opacity maps used to include current data range
clip1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
clip1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'CLV3_Peptide'
cLV3_PeptideLUT = GetColorTransferFunction('CLV3_Peptide')
cLV3_PeptideLUT.AutomaticRescaleRangeMode = "Grow and update on 'Apply'"
cLV3_PeptideLUT.InterpretValuesAsCategories = 0
cLV3_PeptideLUT.AnnotationsInitialized = 0
cLV3_PeptideLUT.ShowCategoricalColorsinDataRangeOnly = 0
cLV3_PeptideLUT.RescaleOnVisibilityChange = 0
cLV3_PeptideLUT.EnableOpacityMapping = 0
cLV3_PeptideLUT.RGBPoints = [0.0886145, 0.231373, 0.298039, 0.752941, 13.06280725, 0.865003, 0.865003, 0.865003, 26.037, 0.705882, 0.0156863, 0.14902]
cLV3_PeptideLUT.UseLogScale = 0
cLV3_PeptideLUT.ColorSpace = 'Diverging'
cLV3_PeptideLUT.UseBelowRangeColor = 0
cLV3_PeptideLUT.BelowRangeColor = [0.0, 0.0, 0.0]
cLV3_PeptideLUT.UseAboveRangeColor = 0
cLV3_PeptideLUT.AboveRangeColor = [0.5, 0.5, 0.5]
cLV3_PeptideLUT.NanColor = [1.0, 1.0, 0.0]
cLV3_PeptideLUT.NanOpacity = 1.0
cLV3_PeptideLUT.Discretize = 1
cLV3_PeptideLUT.NumberOfTableValues = 256
cLV3_PeptideLUT.ScalarRangeInitialized = 1.0
cLV3_PeptideLUT.HSVWrap = 0
cLV3_PeptideLUT.VectorComponent = 0
cLV3_PeptideLUT.VectorMode = 'Magnitude'
cLV3_PeptideLUT.AllowDuplicateScalars = 1
cLV3_PeptideLUT.Annotations = []
cLV3_PeptideLUT.ActiveAnnotatedValues = []
cLV3_PeptideLUT.IndexedColors = []
cLV3_PeptideLUT.IndexedOpacities = []

# get opacity transfer function/opacity map for 'CLV3_Peptide'
cLV3_PeptidePWF = GetOpacityTransferFunction('CLV3_Peptide')
cLV3_PeptidePWF.Points = [0.0886145, 0.0, 0.5, 0.0, 26.037, 1.0, 0.5, 0.0]
cLV3_PeptidePWF.AllowDuplicateScalars = 1
cLV3_PeptidePWF.UseLogScale = 0
cLV3_PeptidePWF.ScalarRangeInitialized = 1

# Rescale transfer function
cLV3_PeptideLUT.RescaleTransferFunction(0, 250.0)

# Rescale transfer function
cLV3_PeptidePWF.RescaleTransferFunction(0, 250.0)

# set active source
SetActiveSource(tableToPoints1)

animationScene1.Play()

# current camera placement for renderView1
renderView1.CameraPosition = [90.3652515703269, -0.056819915771484375, 4.81049570441246]
renderView1.CameraFocalPoint = [-0.031060000000000088, -0.056819915771484375, 4.81049570441246]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraParallelScale = 23.817211739818177

# save screenshot
SaveScreenshot('/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Pics/CLV3PScreenShot_wtresults021621DcTrWt.png', renderView1, SaveAllViews=0,
    ImageResolution=[1896, 3168],
    FontScaling='Scale fonts proportionally',
    SeparatorWidth=1,
    SeparatorColor=[0.937, 0.922, 0.906],
    OverrideColorPalette='',
    StereoMode='No change',
    TransparentBackground=0, 
    # PNG options
    CompressionLevel='5')

# set active source
SetActiveSource(clip1)

# set scalar coloring
ColorBy(clip1Display, ('POINTS', 'Dimer'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(cLV3_PeptideLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
clip1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
clip1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'Dimer'
dimerLUT = GetColorTransferFunction('Dimer')
dimerLUT.AutomaticRescaleRangeMode = "Grow and update on 'Apply'"
dimerLUT.InterpretValuesAsCategories = 0
dimerLUT.AnnotationsInitialized = 0
dimerLUT.ShowCategoricalColorsinDataRangeOnly = 0
dimerLUT.RescaleOnVisibilityChange = 0
dimerLUT.EnableOpacityMapping = 0
dimerLUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 4.0, 0.865003, 0.865003, 0.865003, 8.0, 0.705882, 0.0156863, 0.14902]
dimerLUT.UseLogScale = 0
dimerLUT.ColorSpace = 'Diverging'
dimerLUT.UseBelowRangeColor = 0
dimerLUT.BelowRangeColor = [0.0, 0.0, 0.0]
dimerLUT.UseAboveRangeColor = 0
dimerLUT.AboveRangeColor = [0.5, 0.5, 0.5]
dimerLUT.NanColor = [1.0, 1.0, 0.0]
dimerLUT.NanOpacity = 1.0
dimerLUT.Discretize = 1
dimerLUT.NumberOfTableValues = 256
dimerLUT.ScalarRangeInitialized = 1.0
dimerLUT.HSVWrap = 0
dimerLUT.VectorComponent = 0
dimerLUT.VectorMode = 'Magnitude'
dimerLUT.AllowDuplicateScalars = 1
dimerLUT.Annotations = []
dimerLUT.ActiveAnnotatedValues = []
dimerLUT.IndexedColors = []
dimerLUT.IndexedOpacities = []

# get opacity transfer function/opacity map for 'Dimer'
dimerPWF = GetOpacityTransferFunction('Dimer')
dimerPWF.Points = [0.0, 0.0, 0.5, 0.0, 8.0, 1.0, 0.5, 0.0]
dimerPWF.AllowDuplicateScalars = 1
dimerPWF.UseLogScale = 0
dimerPWF.ScalarRangeInitialized = 1

# Rescale transfer function
dimerLUT.RescaleTransferFunction(0.0, 10.0)

# Rescale transfer function
dimerPWF.RescaleTransferFunction(0.0, 10.0)

animationScene1.Play()

# set active source
SetActiveSource(tableToPoints1)

# current camera placement for renderView1
renderView1.CameraPosition = [90.3652515703269, -0.056819915771484375, 4.81049570441246]
renderView1.CameraFocalPoint = [-0.031060000000000088, -0.056819915771484375, 4.81049570441246]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraParallelScale = 23.817211739818177

# save screenshot
SaveScreenshot('/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Pics/DimerScreenShot_wtresults021621DcTrWt.png', renderView1, SaveAllViews=0,
    ImageResolution=[1896, 3168],
    FontScaling='Scale fonts proportionally',
    SeparatorWidth=1,
    SeparatorColor=[0.937, 0.922, 0.906],
    OverrideColorPalette='',
    StereoMode='No change',
    TransparentBackground=0, 
    # PNG options
    CompressionLevel='5')

# set active source
SetActiveSource(clip1)

# set scalar coloring
ColorBy(clip1Display, ('POINTS', 'Monomer'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(dimerLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
clip1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
clip1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'Monomer'
monomerLUT = GetColorTransferFunction('Monomer')
monomerLUT.AutomaticRescaleRangeMode = "Grow and update on 'Apply'"
monomerLUT.InterpretValuesAsCategories = 0
monomerLUT.AnnotationsInitialized = 0
monomerLUT.ShowCategoricalColorsinDataRangeOnly = 0
monomerLUT.RescaleOnVisibilityChange = 0
monomerLUT.EnableOpacityMapping = 0
monomerLUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 3.5, 0.865003, 0.865003, 0.865003, 7.0, 0.705882, 0.0156863, 0.14902]
monomerLUT.UseLogScale = 0
monomerLUT.ColorSpace = 'Diverging'
monomerLUT.UseBelowRangeColor = 0
monomerLUT.BelowRangeColor = [0.0, 0.0, 0.0]
monomerLUT.UseAboveRangeColor = 0
monomerLUT.AboveRangeColor = [0.5, 0.5, 0.5]
monomerLUT.NanColor = [1.0, 1.0, 0.0]
monomerLUT.NanOpacity = 1.0
monomerLUT.Discretize = 1
monomerLUT.NumberOfTableValues = 256
monomerLUT.ScalarRangeInitialized = 1.0
monomerLUT.HSVWrap = 0
monomerLUT.VectorComponent = 0
monomerLUT.VectorMode = 'Magnitude'
monomerLUT.AllowDuplicateScalars = 1
monomerLUT.Annotations = []
monomerLUT.ActiveAnnotatedValues = []
monomerLUT.IndexedColors = []
monomerLUT.IndexedOpacities = []

# get opacity transfer function/opacity map for 'Monomer'
monomerPWF = GetOpacityTransferFunction('Monomer')
monomerPWF.Points = [0.0, 0.0, 0.5, 0.0, 7.0, 1.0, 0.5, 0.0]
monomerPWF.AllowDuplicateScalars = 1
monomerPWF.UseLogScale = 0
monomerPWF.ScalarRangeInitialized = 1

# Rescale transfer function
monomerLUT.RescaleTransferFunction(0.0, 10.0)

# Rescale transfer function
monomerPWF.RescaleTransferFunction(0.0, 10.0)

animationScene1.Play()

# set active source
SetActiveSource(tableToPoints1)

# current camera placement for renderView1
renderView1.CameraPosition = [90.3652515703269, -0.056819915771484375, 4.81049570441246]
renderView1.CameraFocalPoint = [-0.031060000000000088, -0.056819915771484375, 4.81049570441246]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraParallelScale = 23.817211739818177

# save screenshot
SaveScreenshot('/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Pics/MonomerScreenShot_wtresults021621DcTrWt.png', renderView1, SaveAllViews=0,
    ImageResolution=[1896, 3168],
    FontScaling='Scale fonts proportionally',
    SeparatorWidth=1,
    SeparatorColor=[0.937, 0.922, 0.906],
    OverrideColorPalette='',
    StereoMode='No change',
    TransparentBackground=0, 
    # PNG options
    CompressionLevel='5')

# set active source
SetActiveSource(clip1)

# set scalar coloring
ColorBy(clip1Display, ('POINTS', 'WUSCyto'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(monomerLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
clip1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
clip1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'WUSCyto'
wUSCytoLUT = GetColorTransferFunction('WUSCyto')
wUSCytoLUT.AutomaticRescaleRangeMode = "Grow and update on 'Apply'"
wUSCytoLUT.InterpretValuesAsCategories = 0
wUSCytoLUT.AnnotationsInitialized = 0
wUSCytoLUT.ShowCategoricalColorsinDataRangeOnly = 0
wUSCytoLUT.RescaleOnVisibilityChange = 0
wUSCytoLUT.EnableOpacityMapping = 0
wUSCytoLUT.RGBPoints = [0.0558413, 0.231373, 0.298039, 0.752941, 27.60207065, 0.865003, 0.865003, 0.865003, 55.1483, 0.705882, 0.0156863, 0.14902]
wUSCytoLUT.UseLogScale = 0
wUSCytoLUT.ColorSpace = 'Diverging'
wUSCytoLUT.UseBelowRangeColor = 0
wUSCytoLUT.BelowRangeColor = [0.0, 0.0, 0.0]
wUSCytoLUT.UseAboveRangeColor = 0
wUSCytoLUT.AboveRangeColor = [0.5, 0.5, 0.5]
wUSCytoLUT.NanColor = [1.0, 1.0, 0.0]
wUSCytoLUT.NanOpacity = 1.0
wUSCytoLUT.Discretize = 1
wUSCytoLUT.NumberOfTableValues = 256
wUSCytoLUT.ScalarRangeInitialized = 1.0
wUSCytoLUT.HSVWrap = 0
wUSCytoLUT.VectorComponent = 0
wUSCytoLUT.VectorMode = 'Magnitude'
wUSCytoLUT.AllowDuplicateScalars = 1
wUSCytoLUT.Annotations = []
wUSCytoLUT.ActiveAnnotatedValues = []
wUSCytoLUT.IndexedColors = []
wUSCytoLUT.IndexedOpacities = []

# get opacity transfer function/opacity map for 'WUSCyto'
wUSCytoPWF = GetOpacityTransferFunction('WUSCyto')
wUSCytoPWF.Points = [0.0558413, 0.0, 0.5, 0.0, 55.1483, 1.0, 0.5, 0.0]
wUSCytoPWF.AllowDuplicateScalars = 1
wUSCytoPWF.UseLogScale = 0
wUSCytoPWF.ScalarRangeInitialized = 1

# Rescale transfer function
wUSCytoLUT.RescaleTransferFunction(0, 150.0)

# Rescale transfer function
wUSCytoPWF.RescaleTransferFunction(0, 150.0)

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
wUSCytoLUT.ApplyPreset('erdc_purple2green_dark', True)

animationScene1.Play()

# set active source
SetActiveSource(tableToPoints1)

# current camera placement for renderView1
renderView1.CameraPosition = [90.3652515703269, -0.056819915771484375, 4.81049570441246]
renderView1.CameraFocalPoint = [-0.031060000000000088, -0.056819915771484375, 4.81049570441246]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraParallelScale = 23.817211739818177

# save screenshot
SaveScreenshot('/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Pics/WUSCScreenShot_wtresults021621DcTrWt.png', renderView1, SaveAllViews=0,
    ImageResolution=[1896, 3168],
    FontScaling='Scale fonts proportionally',
    SeparatorWidth=1,
    SeparatorColor=[0.937, 0.922, 0.906],
    OverrideColorPalette='',
    StereoMode='No change',
    TransparentBackground=0, 
    # PNG options
    CompressionLevel='5')

# set active source
SetActiveSource(clip1)

# set scalar coloring
ColorBy(clip1Display, ('POINTS', 'WUSNuc_WS'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(wUSCytoLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
clip1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
clip1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'WUSNuc_WS'
wUSNuc_WSLUT = GetColorTransferFunction('WUSNuc_WS')
wUSNuc_WSLUT.AutomaticRescaleRangeMode = "Grow and update on 'Apply'"
wUSNuc_WSLUT.InterpretValuesAsCategories = 0
wUSNuc_WSLUT.AnnotationsInitialized = 0
wUSNuc_WSLUT.ShowCategoricalColorsinDataRangeOnly = 0
wUSNuc_WSLUT.RescaleOnVisibilityChange = 0
wUSNuc_WSLUT.EnableOpacityMapping = 0
wUSNuc_WSLUT.RGBPoints = [0.0389217, 0.231373, 0.298039, 0.752941, 86.28796085, 0.865003, 0.865003, 0.865003, 172.537, 0.705882, 0.0156863, 0.14902]
wUSNuc_WSLUT.UseLogScale = 0
wUSNuc_WSLUT.ColorSpace = 'Diverging'
wUSNuc_WSLUT.UseBelowRangeColor = 0
wUSNuc_WSLUT.BelowRangeColor = [0.0, 0.0, 0.0]
wUSNuc_WSLUT.UseAboveRangeColor = 0
wUSNuc_WSLUT.AboveRangeColor = [0.5, 0.5, 0.5]
wUSNuc_WSLUT.NanColor = [1.0, 1.0, 0.0]
wUSNuc_WSLUT.NanOpacity = 1.0
wUSNuc_WSLUT.Discretize = 1
wUSNuc_WSLUT.NumberOfTableValues = 256
wUSNuc_WSLUT.ScalarRangeInitialized = 1.0
wUSNuc_WSLUT.HSVWrap = 0
wUSNuc_WSLUT.VectorComponent = 0
wUSNuc_WSLUT.VectorMode = 'Magnitude'
wUSNuc_WSLUT.AllowDuplicateScalars = 1
wUSNuc_WSLUT.Annotations = []
wUSNuc_WSLUT.ActiveAnnotatedValues = []
wUSNuc_WSLUT.IndexedColors = []
wUSNuc_WSLUT.IndexedOpacities = []

# get opacity transfer function/opacity map for 'WUSNuc_WS'
wUSNuc_WSPWF = GetOpacityTransferFunction('WUSNuc_WS')
wUSNuc_WSPWF.Points = [0.0389217, 0.0, 0.5, 0.0, 172.537, 1.0, 0.5, 0.0]
wUSNuc_WSPWF.AllowDuplicateScalars = 1
wUSNuc_WSPWF.UseLogScale = 0
wUSNuc_WSPWF.ScalarRangeInitialized = 1

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
wUSNuc_WSLUT.ApplyPreset('erdc_purple2green_dark', True)

# Rescale transfer function
wUSNuc_WSLUT.RescaleTransferFunction(0.0, 250.0)

# Rescale transfer function
wUSNuc_WSPWF.RescaleTransferFunction(0.0, 250.0)

animationScene1.Play()

# set active source
SetActiveSource(tableToPoints1)

# current camera placement for renderView1
renderView1.CameraPosition = [90.3652515703269, -0.056819915771484375, 4.81049570441246]
renderView1.CameraFocalPoint = [-0.031060000000000088, -0.056819915771484375, 4.81049570441246]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraParallelScale = 23.817211739818177

# save screenshot
SaveScreenshot('/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Pics/WUSNScreenShot_wtresults021621DcTrWt.png', renderView1, SaveAllViews=0,
    ImageResolution=[1896, 3168],
    FontScaling='Scale fonts proportionally',
    SeparatorWidth=1,
    SeparatorColor=[0.937, 0.922, 0.906],
    OverrideColorPalette='',
    StereoMode='No change',
    TransparentBackground=0, 
    # PNG options
    CompressionLevel='5')

# set active source
SetActiveSource(clip1)

# set scalar coloring
ColorBy(clip1Display, ('POINTS', 'WUSRNA'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(wUSNuc_WSLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
clip1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
clip1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'WUSRNA'
wUSRNALUT = GetColorTransferFunction('WUSRNA')
wUSRNALUT.AutomaticRescaleRangeMode = "Grow and update on 'Apply'"
wUSRNALUT.InterpretValuesAsCategories = 0
wUSRNALUT.AnnotationsInitialized = 0
wUSRNALUT.ShowCategoricalColorsinDataRangeOnly = 0
wUSRNALUT.RescaleOnVisibilityChange = 0
wUSRNALUT.EnableOpacityMapping = 0
wUSRNALUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 12.4997, 0.865003, 0.865003, 0.865003, 24.9994, 0.705882, 0.0156863, 0.14902]
wUSRNALUT.UseLogScale = 0
wUSRNALUT.ColorSpace = 'Diverging'
wUSRNALUT.UseBelowRangeColor = 0
wUSRNALUT.BelowRangeColor = [0.0, 0.0, 0.0]
wUSRNALUT.UseAboveRangeColor = 0
wUSRNALUT.AboveRangeColor = [0.5, 0.5, 0.5]
wUSRNALUT.NanColor = [1.0, 1.0, 0.0]
wUSRNALUT.NanOpacity = 1.0
wUSRNALUT.Discretize = 1
wUSRNALUT.NumberOfTableValues = 256
wUSRNALUT.ScalarRangeInitialized = 1.0
wUSRNALUT.HSVWrap = 0
wUSRNALUT.VectorComponent = 0
wUSRNALUT.VectorMode = 'Magnitude'
wUSRNALUT.AllowDuplicateScalars = 1
wUSRNALUT.Annotations = []
wUSRNALUT.ActiveAnnotatedValues = []
wUSRNALUT.IndexedColors = []
wUSRNALUT.IndexedOpacities = []

# get opacity transfer function/opacity map for 'WUSRNA'
wUSRNAPWF = GetOpacityTransferFunction('WUSRNA')
wUSRNAPWF.Points = [0.0, 0.0, 0.5, 0.0, 24.9994, 1.0, 0.5, 0.0]
wUSRNAPWF.AllowDuplicateScalars = 1
wUSRNAPWF.UseLogScale = 0
wUSRNAPWF.ScalarRangeInitialized = 1

# Rescale transfer function
wUSRNALUT.RescaleTransferFunction(0.0, 80.0)

# Rescale transfer function
wUSRNAPWF.RescaleTransferFunction(0.0, 80.0)

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
wUSRNALUT.ApplyPreset('erdc_purple2green_dark', True)

# set active source
SetActiveSource(tableToPoints1)

animationScene1.Play()

# current camera placement for renderView1
renderView1.CameraPosition = [90.3652515703269, -0.056819915771484375, 4.81049570441246]
renderView1.CameraFocalPoint = [-0.031060000000000088, -0.056819915771484375, 4.81049570441246]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraParallelScale = 23.817211739818177

# save screenshot
SaveScreenshot('/home/ace/Documents/Software_Projects/OrgSimMerc_IDE/Results/Pics/WUSRNAScreenShot_wtresults021621DcTrWt.png', renderView1, SaveAllViews=0,
    ImageResolution=[1896, 3168],
    FontScaling='Scale fonts proportionally',
    SeparatorWidth=1,
    SeparatorColor=[0.937, 0.922, 0.906],
    OverrideColorPalette='',
    StereoMode='No change',
    TransparentBackground=0, 
    # PNG options
    CompressionLevel='5')

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [90.3652515703269, -0.056819915771484375, 4.81049570441246]
renderView1.CameraFocalPoint = [-0.031060000000000088, -0.056819915771484375, 4.81049570441246]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraParallelScale = 23.817211739818177


#save state
SaveState("./Results/Paraview_State_Files/wtresults021621DcTrWt.pvsm")

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
