#pragma once
#include <armadillo>

using namespace arma;
int main(int argc, char** argv);
double evaluateKrigingModel(mat data, vec kring, double reg_param, vec input);