import sys
import os
sys.path.append(os.environ['SU2_RUN']) 
import SU2
from optparse import OptionParser
from pyDOE import lhs
import shutil
import numpy
import matplotlib.pyplot as plt
import re

from mesh_deformation import mesh_deformation
from parallel_computation import parallel_computation

# from paraview.simple import *

def data_aquisition(config_file, scaling = 1, sample_dir = 'InitialSamples', partitions=0, nzones=1, quiet=False, samples = 400, use_coordinates = False, use_pressure = True, use_density = True, use_momentum = True):
    if not config_file:
        sys.stdout.write('No config file provided. Aborting.\n')
        sys.exit(1)

    config = SU2.io.Config(config_file)
    
    if config.OUTPUT_FORMAT != 'PARAVIEW':
        sys.stdout.write('Format has to be PARAVIEW. Aborting. \n')
        sys.exit(1)

    config.NUMBER_PART = partitions
    config.NZONES = int(nzones)
    if quiet: config.CONSOLE = 'CONCISE'
    
    # load the DV_DEFINITION, so we can write the proper deformation settings in the configs for the samples
    def_dv = config.DEFINITION_DV                               # complete definition of the desing variable

    n_dv = sum(def_dv['SIZE'])                                # number of design variables
    bound_upper      = float ( config.OPT_BOUND_UPPER )*scaling                   # variable bound to be scaled by the line search
    bound_lower      = float ( config.OPT_BOUND_LOWER )*scaling                   # variable bound to be scaled by the line search
    x0 = [0.0] * n_dv # initial design
    xb_low           = [float(bound_lower)]*n_dv      # lower dv bound it includes the line search acceleration factor
    xb_up            = [float(bound_upper)]*n_dv      # upper dv bound it includes the line search acceleration fa
    xb          = list(zip(xb_low, xb_up)) # design bounds
    mesh_in = config.MESH_FILENAME
    mesh_out = config.MESH_OUT_FILENAME

    flowfilename = config.VOLUME_FLOW_FILENAME + ".vtk"

    kindstring = ""
    for kind in def_dv['KIND']:
        if kindstring:
            kindstring += ", "
        kindstring += kind

    dvparamstring = ""
    for param in def_dv['PARAM']:
        if dvparamstring:
            dvparamstring += ' ; '
        dvparamstring += "( %f, %f)"%(param[0],param[1])

    # State
    state = SU2.io.State()
    state.find_files(config)
    
    sys.stdout.write('Generating LHS Points Distribution.\n')

    # compute a lhs sample
    lhs_sample = lhs(n_dv, samples)

    # transform lhs sample to design space
    dv_values_array = lhs_sample*(bound_upper - bound_lower) + bound_lower

    # write config files containing the deformations
    if os.path.exists(sample_dir):
        if query_yes_no('The given sample dir already exists. Any content will be deleted. Do you want to continue?'):
            shutil.rmtree(sample_dir)
        else: 
            sys.stdout.write('Aborting.\n')
            sys.exit(1)
            
    cfg_input_file = open(config_file,'r')
    cfg_input = cfg_input_file.read()
    cfg_input_file.close()

    sys.stdout.write('Writing Config Files.\n')

    dvzerostring = ''
    for _ in dv_values_array[0]:
        if dvzerostring:
            dvzerostring += ', '
        dvzerostring +='0.0'

    for i in range(samples):
        os.makedirs(os.path.join(sample_dir,'Sample'+str(i)))
        config_DEF = open(os.path.join(sample_dir,'Sample'+str(i),'config_DEF.cfg'), 'w')
        first = True
        dvstring = ''
        for dv_value in dv_values_array[i]:
            if not first:
                dvstring += ', '
            first = False
            dvstring += str(dv_value)
        config_DEF.write(cfg_input.replace('DV_VALUE','DV_VALUE='+dvstring+'\n%').replace('DV_PARAM','DV_PARAM='+dvparamstring+'\n%').replace('DV_KIND','DV_KIND='+kindstring+'\n%')) # comment out stuff we want to change, in case it exists
        # write deformation in config file for sample
        # config_DEF.write('\nDV_VALUE=')
        # config_DEF.write(dvstring)
        config_DEF.write('\nDV_VALUE_NEW=')
        config_DEF.write(dvstring)
        config_DEF.write('\nDV_VALUE_OLD=')
        config_DEF.write(dvzerostring)
        #config_DEF.write('\nDV_KIND=')
        #config_DEF.write(kindstring)
        #config_DEF.write('\nDV_PARAM=')
        #config_DEF.write(dvparamstring)
        config_DEF.write('\n')
        config_DEF.close()
        
    sys.stdout.write('Save Config Settings.\n')

    # save the sample values for each sample
    for i in range(samples):
        numpy.savetxt(os.path.join(sample_dir,'Sample'+str(i),'snapshot_settings.txt'),dv_values_array[i])
               
    sys.stdout.write('Performing Mesh Deformations.\n')

    # do the mesh deform
    for i in range(samples):
        sys.stdout.write('Mesh Deformation %d/%d.\n'%(i,samples))
        mesh_deformation(os.path.join(sample_dir,'Sample'+str(i),'config_DEF.cfg'))
        os.rename(mesh_out,os.path.join(sample_dir,'Sample'+str(i),'deformed_mesh.su2'))
                
    sys.stdout.write('Start Simulations.\n')

    # run simulations
    os.rename(mesh_in,mesh_in+'.orig')
    for i in range(samples):
        sys.stdout.write('Starting Simulation %d/%d.\n'%(i,samples))
        os.rename(os.path.join(sample_dir,'Sample'+str(i),'deformed_mesh.su2'), mesh_in)
        parallel_computation(config_file,partitions=partitions)
        os.rename(flowfilename,os.path.join(sample_dir,'Sample'+str(i),flowfilename))
        os.rename(mesh_in, os.path.join(sample_dir,'Sample'+str(i),'deformed_mesh.su2'))
    os.rename(mesh_in+'.orig',mesh_in)

    # run simulation for undeformed mesh
    parallel_computation(config_file, partitions=partitions)

    # generate a template (containing basic information like connection between mesh points. Used to create flow files from surrogate output)
    generate_template(flowfilename, sample_dir, use_coordinates = use_coordinates, use_pressure = use_pressure, use_density = use_density, use_momentum = use_momentum)

    generate_snapshot_matrix(samples, flowfilename = flowfilename, sample_dir = sample_dir, 
                             use_coordinates = use_coordinates, use_pressure = use_pressure, use_density = use_density, use_momentum = use_momentum)

    
def generate_template(flowfile, sample_dir, use_coordinates = False, use_pressure = True, use_density = True, use_momentum = True):
    
    flow_file = open(flowfile,'r')
    template_text = flow_file.read()
    flow_file.close()
    
    template_data = re.search('(?P<fileheader># vtk DataFile .*\n)' + # capture header
                            '(?P<coordinatesheader>POINTS[^\n]*\n)' + # capture coordinates header
                            '(?P<coordinatesdata>.*\n)' +  # capture coordinates (in case they are not part of the surrogate, we'll provide them as template)
                            '(?P<meshstuff>CELLS.*CELL_TYPES.*POINT_DATA[^\n]*\n).*',  # capture mesh data (connection between points)
                            template_text, flags=re.DOTALL)
    
    out_file = open(os.path.join(sample_dir,'flow_file_template.vtk'),'w')
    out_file.write(template_data.group('fileheader'))
    out_file.write(template_data.group('coordinatesheader'))
    if use_coordinates:
        out_file.write('%COORDINATES%') # if coordinates are part of the surrogate, they will be inserted as surrogate result
    else: 
        out_file.write(template_data.group('coordinatesdata')) # otherwise we use an undeformed mesh
    out_file.write(template_data.group('meshstuff'))
    if use_density:
        out_file.write('SCALARS Density float 1\nLOOKUP_TABLE default\n%DENSITY%')
    if use_momentum:
        out_file.write('SCALARS X-Momentum float 1\nLOOKUP_TABLE default\n%MOMENTUMX%')
        out_file.write('SCALARS Y-Momentum float 1\nLOOKUP_TABLE default\n%MOMENTUMY%')
    if use_pressure:
        out_file.write('SCALARS Pressure float 1\nLOOKUP_TABLE default\n%PRESSURE%')
    out_file.close()
    
def generate_snapshot_matrix(samples, flowfilename = 'flow.vtk', sample_dir = 'InitialSamples', use_coordinates = False, use_pressure = True, use_density = True, use_momentum = True):
    sys.stdout.write('Generating snapshot matrix...\n')
    out_file = open(os.path.join(sample_dir,'snapshot_matrix.txt'),'w')
    for sample in range(samples):
        sys.stdout.write(str(sample+1)+'/'+str(samples)+'\n')
        sample_file = open(os.path.join(sample_dir,'Sample'+str(sample),flowfilename))
        sample_data = sample_file.read().replace('\t',' ')
        sample_file.close()
        if sample:
            out_file.write(';\n')
        sample_data = re.search('POINTS [0-9]+[^\n]*\n(?P<coordinates>[^\n]*) \n.*' +
                                'Density[^\n]*\nLOOKUP_TABLE[^\n]*\n(?P<density>[^\n]*) \n.*' + 
                                'X-Momentum[^\n]*\nLOOKUP_TABLE[^\n]*\n(?P<momentum_x>[^\n]*) \n.*' +
                                'Y-Momentum[^\n]*\nLOOKUP_TABLE[^\n]*\n(?P<momentum_y>[^\n]*) \n.*' + 
                                'Pressure [^\n]*\nLOOKUP_TABLE[^\n]*\n(?P<pressure>[^\n]*) [\n,\Z].*', sample_data, flags=re.DOTALL)
        if use_coordinates:
            # remove all z-coordinates, as we only are 2D
            coordinates = sample_data.group('coordinates').split(' ')
            coordinates_without_z=""
            for index, val in enumerate(coordinates):
                if index%3 == 2:
                    continue
                coordinates_without_z += val + ' '

            out_file.write(coordinates_without_z)
        
        if use_density:
            out_file.write(sample_data.group('density'))
            out_file.write(' ')
        if use_momentum:
            out_file.write(sample_data.group('momentum_x'))
            out_file.write(' ')
            out_file.write(sample_data.group('momentum_y'))
            out_file.write(' ')
        if use_pressure:
            out_file.write(sample_data.group('pressure'))
    out_file.close()

# Method from https://stackoverflow.com/questions/3041986/apt-command-line-interface-like-yes-no-input
def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is True for "yes" or False for "no".
    """
    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")


def main(): 

    parser = OptionParser()
    parser.add_option("-f", "--file", dest="filename",
                      help="read config from FILE", metavar="FILE")
    parser.add_option("-p", "--samplepath", dest="sampledir",
                      help="store samples in PATH (directory)", metavar="SAMPLE PATH")
    parser.add_option("-n", "--partitions", dest="partitions", default=1,
                      help="number of PARTITIONS", metavar="PARTITIONS")
    parser.add_option("-q", "--quiet", dest="quiet", default="True",
                      help="True/False Quiet SU2 output", metavar="QUIET")
    parser.add_option("-r", "--pressure", dest="pressure", default="True",
                      help="True/False include pressure", metavar="PRESSURE")
    parser.add_option("-m", "--momentum", dest="momentum", default="True",
                      help="True/False include momentum", metavar="MOMENTUM")
    parser.add_option("-c", "--coordinates", dest="coordinates", default="False",
                      help="True/False include coordinates", metavar="COORDINATES")
    parser.add_option("-d", "--density", dest="density", default="True",
                      help="True/False include density", metavar="DENSITY")
    parser.add_option("-z", "--zones", dest="nzones", default="1",
                      help="Number of Zones", metavar="ZONES")
    parser.add_option("-s", "--samples", dest="samples", default="100", 
                      help="Number of Initial Samples", metavar="SAMPLES")
    parser.add_option("-g", "--deformationscaling", dest="scaling", default="1", 
                      help="scaling applied to OPT_BOUND (config file) as limits for sampling region", metavar="SCALING")

    (options, args) = parser.parse_args()
    
    # process inputs
    options.partitions = int(options.partitions)
    options.quiet = options.quiet.upper() == 'TRUE'
    options.density = options.density.upper() == 'TRUE'
    options.coordinates = options.coordinates.upper() == 'TRUE'
    options.momentum = options.momentum.upper() == 'TRUE'
    options.pressure = options.pressure.upper() == 'TRUE'
    options.nzones = int(options.nzones)
    options.samples = int(options.samples)
    options.scaling = float(options.scaling)

    data_aquisition(options.filename, scaling = options.scaling, sample_dir=options.sampledir, partitions=options.partitions, nzones=options.nzones, quiet=options.quiet, samples = options.samples, use_coordinates = options.coordinates, use_pressure = options.pressure, use_density = options.density, use_momentum = options.density)

    generate_snapshot_matrix(options.samples, sample_dir = options.sampledir, use_coordinates = options.coordinates, use_pressure = options.pressure, use_density = options.density, use_momentum = options.density)



if __name__ == '__main__':
    main()
