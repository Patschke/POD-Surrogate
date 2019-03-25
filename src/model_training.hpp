#pragma once
#include <armadillo>

using namespace arma;
int main(int argc, char** argv);
void trainMode(std::string input_file, std::string output_file, std::string hyperparamfile, double regparam, int max_iterations);
