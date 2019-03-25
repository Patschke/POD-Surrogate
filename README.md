# Implementation of a POD based surrogate model for aerodynamic shape optimization

Please refer to the thesis for a in depth description of the methods involved. 

Make sure the imports of SU2, Armadillo, PyDOE and RoDeO within the source files work for your setup. 


# TL;DR
0. Requirements: SU2, RoDeO, Armadillo, PyDOE a SU2 compatible config file and mesh for your model. 
1. Compile model_order_reduction.hpp, model_training.cpp, model_evaluation.cpp and model_order_decompression
2. Surrogate model creation:
```
python acquire_samples.py -f [CONFIG.cfg] -p [SAMPLE_DIR] -s [SAMPLENR]
mkdir [REDUCED_MODEL_DIR]
model_order_reduction -i[SAMPLE_DIR] -o[REDUCED_MODEL_DIR] -n[MODECOUNT]
python prepare_training_data.py -f [REDUCED_MODEL_DIR]/PODCoefficients.txt -p [SAMPLE_DIR] -o [TRAINING_DIR]
mkdir [KRIGING_DIR]
model_training -i[TRAINING_DIR] -o[KRIGING_DIR] -n[MODECOUNT] -i[ITERATIONS] -r[REGPARAM]
```
3. Surrogate model evaluation (INPUT.txt containing line separated DV_VALUES)
```
model_evaluation -i[INPUT.txt] -m[KRIGIN_DIR] -p[REDUCED_MODEL_DIR] -o[FLOW_VECTOR.txt]
python  -i [FLOW_VECTOR.txt] -t [SAMPLE_DIR]/flow_file_template.vtk -o [OUTPUT_FLOWFILE.vtk] 
```
4. Enjoy your surrogate model