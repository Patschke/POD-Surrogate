import sys
import os
from optparse import OptionParser
import re

from paraview.simple import *

def main(): 

    parser = OptionParser()
    parser.add_option("-i", "--inputvector", dest="input",
                      help="read ascii snapshot vector FILE", metavar="FILE")
    parser.add_option("-t", "--template", dest="template",
                      help="use template FILE", metavar="FILE")
    parser.add_option("-o", "--output", dest="output",
                      help="save output flowfile as FILE", metavar="FILE")
    parser.add_option("-r", "--render", dest="render", default="",
                      help="optional: outputfile for rendered view", metavar="FILE")

    (options, args) = parser.parse_args()

    reconstructFlowFile(options.input, options.template, options.output)
    if not options.render == "":
        visualizeFlow(options.output, options.render)

# Note: Inputvector must contain variables in following order to be matched correctly: Coordinates, Density, MomentumX, MomentumY, Pressure.
def reconstructFlowFile(inputfile, templatefile, outfile):
    flow_template_file = open(templatefile, 'r')
    flow_template = flow_template_file.read()
    flow_template_file.close()

    number_of_properties = 0 # parse template to figure out how many properties are included in the inputvector

    use_coordinates = False
    use_density = False
    use_momentumx = False
    use_momentumy = False
    use_pressure = False
    if "%COORDINATES%" in flow_template:
        number_of_properties +=2 # assuming we are 2D
        use_coordinates = True
    if "%DENSITY%" in flow_template:
        number_of_properties += 1
        use_density = True
    if "%MOMENTUMX%" in flow_template:
        number_of_properties += 1
        use_momentumx = True
    if "%MOMENTUMY%" in flow_template:
        number_of_properties += 1
        use_momentumy = True
    if "%PRESSURE%" in flow_template:
        number_of_properties += 1
        use_pressure = True
    
    data_file = open(inputfile, 'r')
    data = data_file.read()
    data_file.close()

    # some preprocessing to make armadillo output compatible with vtk notation
    data = data.replace(' ','').replace('\t','').split('\n')

    counter = 0
    valuecount = len(data)/number_of_properties # how many values per propertie?
    
    if use_coordinates: 
        coordinates = ""
        for i in range(valuecount*2):
            coordinates += data[i] + "\t"
            if i % 2 == 1:
                coordinates += "0.0\t" # vtk requires 3D coordinats, so we add 0s in 
        coordinates = coordinates[:-1]
        flow_template = flow_template.replace("%COORDINATES%", coordinates)
        counter += 2
        
    if use_density:
        density = ""
        for i in range(counter*valuecount, (counter+1)*valuecount):
            density += data[i] + "\t"
        density = density[:-1]
        flow_template = flow_template.replace("%DENSITY%", density)
        counter += 1
        
    if use_momentumx:
        momentum_x = ""
        for i in range(counter*valuecount, (counter+1)*valuecount):
            momentum_x += data[i] + "\t"
        momentum_x = momentum_x[:-1]
        flow_template = flow_template.replace("%MOMENTUMX%", momentum_x)
        counter += 1

    if use_momentumy:
        momentum_y = ""
        for i in range(counter*valuecount, (counter+1)*valuecount):
            momentum_y += data[i] + "\t"
        momentum_y = momentum_y[:-1]
        flow_template = flow_template.replace("%MOMENTUMY%", momentum_y)
        counter += 1

    if use_pressure:
        pressure = ""
        for i in range(counter*valuecount, (counter+1)*valuecount):
            pressure += data[i] + "\t"
        pressure = pressure[:-1]
        flow_template = flow_template.replace("%PRESSURE%", pressure)
        counter += 1

    reconstructedFile = open(outfile,'w')
    reconstructedFile.write(flow_template)
    reconstructedFile.close()

# uses paraview to create an image of the flow file
def visualizeFlow(flow_file, output_file):
    Connect()
    #### disable automatic camera reset on 'Show'
    paraview.simple._DisableFirstRenderCameraReset()

    # create a new 'Legacy VTK Reader'
    data = LegacyVTKReader(FileNames=[flow_file])

    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')

    xc = 0.5 # tip of airfoil is at 0/0, so shifting in .5 in x should center on airfoil instead of airfoil tip
    yc = 0
    zc = 0
    width = 4.0
    height = 3.0
    scalar = None
    colormap_min = 0.5
    colormap_max = 1.5

    ratio = height / width
    magnification = 1
    height_p = 855 * magnification
    width_p = int(height_p * 1.0 / ratio / magnification)  

    renderView1.ViewSize = [width_p , height_p]

    # show data in view
    dataDisplay = Show(data, renderView1)
    # trace defaults for the display properties.
    dataDisplay.ColorArrayName = ['POINTS', 'Density']
    # set scalar coloring
    ColorBy(dataDisplay, ('POINTS', 'Density'))
    # rescale color and/or opacity maps used to include current data range
    dataDisplay.RescaleTransferFunctionToDataRange(True)
    dataDisplay.SetScalarBarVisibility(renderView1, True)

    # get color transfer function/color map for 'irradiation'
    irradiationLUT = GetColorTransferFunction(scalar)
    # Rescale transfer function
    irradiationLUT.RescaleTransferFunction(colormap_min, colormap_max)
    
    irradiationLUT.ColorSpace = 'RGB'
    irradiationLUT.NanColor = [0.498039, 0.0, 0.0]

    renderView1.InteractionMode = '2D'
    renderView1.CameraPosition = [xc, yc, 10000.0 + zc]
    renderView1.CameraFocalPoint = [xc, yc, zc]
    renderView1.CameraParallelScale = (height / 2.0)

    # save screenshot
    SaveScreenshot(output_file, magnification=magnification, quality=100, view=renderView1)

    Disconnect()




# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()
