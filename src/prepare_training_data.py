import sys
import os
from optparse import OptionParser

def createTrainingData(pod_coefficients_path, sample_file_path, outputfile):
    # Input: pod_coefficient_path (see calculatePODCoefficients, Matrix of podCoefficients as .txt, separated by lines and "  " (double whitespace))
    # Sample files: Folder containing Folders called SampleX, where X is the Sample Number. Inside those a snapshot_settings.txt, containing the deformation parameters separated by lines
    # Output: creates one .csv file per POD Mode, each containing one line per sample, formatted like "design_var_1, ..., design_var_n, pod_coefficient"
    
    sys.stdout.write('Reading POD coefficients...\n')

    # load coefficients as 2D array
    pod_coefficient_file = open(pod_coefficients_path)
    pod_coefficient_data = ([[j for j in i.lstrip().split('  ')] for i in pod_coefficient_file.read().splitlines()]) 
    pod_coefficient_file.close()

    # create csv output list (list containing strings, that will be saved as .csv)
    csv_out = []
    for _ in pod_coefficient_data[0]:
        csv_out.append('')

    sys.stdout.write('Reading samples\n')
    for sample in range(len(pod_coefficient_data)):
        sys.stdout.write('Processing Sample %d/%d\n' % (sample+1, len(pod_coefficient_data)))
        deformation_parameters_file = open(os.path.join(sample_file_path, 'Sample%d'%sample, 'snapshot_settings.txt'))
        deformation_parameters = deformation_parameters_file.read().replace('\n',', ')
        deformation_parameters_file.close()
        for mode in range(len(pod_coefficient_data[sample])):
            csv_out[mode] += deformation_parameters + pod_coefficient_data[sample][mode] + '\n'
    
    sys.stdout.write('Writing training data\n')
    if not os.path.exists(outputfile):
        sys.stdout.write('Output folder did not exist, creating it for you.\n')
        os.makedirs(outputfile)
    for mode in range(len(csv_out)):
        sys.stdout.write('Writing mode %d/%d\n'%(mode+1, len(csv_out)))
        csv_file = open(os.path.join(outputfile, 'TrainingDataMode%d.csv' % mode), 'w')
        csv_file.write(csv_out[mode])
        csv_file.close()


def main(): 
    
    parser = OptionParser()
    parser.add_option("-f", "--podfile", dest="podfile",
                      help="text file containing matrix of pod COEFFICIENTS FILE", metavar="COEFFICIENTS FILE")
    parser.add_option("-p", "--samplepath", dest="sampledir",
                      help="PATH (directory) containing the sample folder", metavar="SAMPLE PATH")
    parser.add_option("-o", "--output", dest="output",
                      help="folder to store the training data in", metavar="TRAINING PATH")

    (options, args) = parser.parse_args()

    if not (options.sampledir and options.podfile and options.output):
        parser.error('You need to specify all inputs/outputs. Use -h for help.')
        sys.stdout.write('Baeh')
    
    createTrainingData(options.podfile, options.sampledir, options.output)

if __name__ == '__main__':
    main()