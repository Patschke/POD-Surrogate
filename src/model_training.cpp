#pragma once

//#include <iostream>
#include <vector>
#include "model_training.hpp"
#include "kriging_training.hpp"
#include "Rodeo_macros.hpp"
#include<omp.h>

using namespace std;
using namespace arma;


int main(int argc, char** argv) {
	// Read input args
	std::string inputpath = "";// "../POD_Surrogate/FinalDataset/TrainingData/";
	std::string outputpath = "";// "../POD_Surrogate/FinalDataset/KrigingModel/";

	int numberofpodmodes = 20;
	int numberofiterations = 5000;
	double regparam = 2;
	for (int i = 1; i < argc; i++) { // Skip first arg, as it is the program name
		std::string arg = argv[i];
		if (arg.substr(0, 2) == "-i") {
			inputpath = arg.substr(2);
			if (inputpath.back() != '/')
				inputpath += "/";
		}
		else if (arg.substr(0, 2) == "-o") {
			outputpath = arg.substr(2);
			if (outputpath.back() != '/')
				outputpath += "/";
		}
		else if (arg.substr(0, 2) == "-n") {
			numberofpodmodes = std::stoi(arg.substr(2));
		}
		else if (arg.substr(0, 2) == "-i") {
			numberofiterations = std::stoi(arg.substr(2));
		}
		else if (arg.substr(0, 2) == "-r") {
			regparam = std::stod(arg.substr(2));
		}
		else {
			cout << "Unknown option " << arg << "\n";
			cout << "Aborting.";
			return -1;
		}
	}

	// Check if input args have been set
	if (inputpath == "" || outputpath == "") {
		cout << "Usage:\n -isamplepath\n -ooutputpath\n [-n(Number of Modes, default 20)]\n [-i(Number of Training Iterations)]\n [-r(Regression Parameter)]\n No spaces between option and content permitted, e.g. -n20\n Aborting.";
		return -1;
	}

	// Train modes
#pragma omp parallel for
	for (int i = 0; i < numberofpodmodes; i++) {
		if (omp_get_thread_num() == 0) 
			printf("Kriging Model training using %d thread(s)\n", omp_get_num_threads());
			
		trainMode(inputpath + "TrainingDataMode" + std::to_string(i) + ".csv",
			outputpath + "KriginParams" + std::to_string(i) + ".vec",
			outputpath + "Hyperparams" + std::to_string(i) + ".csv", 
			regparam, numberofiterations);
	}

	return 0;
}

void trainMode(std::string input_file, std::string output_file, std::string hyperparamfile, double regparam, int max_iterations) {
	int linear_regression = LINEAR_REGRESSION_OFF;
	vec regression_weights, krigin_params;
	mat data;
	if (!data.load(input_file, csv_ascii)) {
		cout << "Invalid training data " << input_file << "\n";
		cout << "Skipping mode";
		return;
	}
	int dim = data.n_cols - 1;
	krigin_params = zeros(dim * 2);
	int data_file_format = CSV_ASCII;

	// Training Method from RoDeO
	train_kriging_response_surface(input_file, hyperparamfile, linear_regression, regression_weights, krigin_params, regparam, max_iterations, data_file_format);

	//regression_weights.save("RegressionWeights" + std::to_string(podMode) + ".vec"); //Empty, as lin reg is off
	krigin_params.save(output_file);
}


