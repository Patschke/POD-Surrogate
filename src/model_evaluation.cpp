#include <vector>
#include "model_evaluation.hpp"
#include "kriging_training.hpp"
#include "Rodeo_macros.hpp"

using namespace std;
using namespace arma;


/*int main(int bla, char** blabla) {
    evaluate a set of inputs
	char* args[] = { NULL, NULL, NULL};

	for (int sample = 49; sample < 50; sample++) {
		cout << "Sample " << sample << "\n";
		args[1] = _strdup(("-i../POD_Surrogate/FinalDataset/ValidationSamples/FullSimulation/Sample" + std::to_string(sample) + "/snapshot_settings.txt").c_str());
		args[2] = _strdup(("-o../POD_Surrogate/FinalDataset/ValidationSamples/SurrogateOutput/sample" + std::to_string(sample) + ".txt").c_str());
		main7(3, args);
	}

	return 0;
}*/

int main(int argc, char** argv) {
	// Read input args
	std::string inputtraining = "";// "../POD_Surrogate/FinalDataset/TrainingData/";
	std::string inputmodel = "";// "../POD_Surrogate/FinalDataset/KrigingModel/";
	std::string inputpod = "";// "../POD_Surrogate/FinalDataset/POD/";
	std::string inputparameter = "";// "../POD_Surrogate/FinalDataset/ReconstructedFlowfields/Kriging/snapshot_settings.txt";
	std::string outputfile = "";// "../POD_Surrogate/FinalDataset/ReconstructedFlowfields/Kriging/flow20.txt";

	int numberofpodmodes = 20;
	double regparam = 2;
	for (int i = 1; i < argc; i++) { // Skip first arg, as it is the program name
		std::string arg = argv[i];
		if (arg.substr(0, 2) == "-i") {
			inputparameter = arg.substr(2);
		}
		else if (arg.substr(0, 2) == "-t") {
			inputtraining = arg.substr(2);
			if (inputtraining.back() != '/')
				inputtraining += "/";
		}
		else if (arg.substr(0, 2) == "-m") {
			inputmodel = arg.substr(2);
			if (inputmodel.back() != '/')
				inputmodel += "/";
		}
		else if (arg.substr(0, 2) == "-o") {
			outputfile = arg.substr(2);
		}
		else if (arg.substr(0, 2) == "-p") {
			inputpod = arg.substr(2);
			if (inputpod.back() != '/')
				inputpod += "/";
		}
		else if (arg.substr(0, 2) == "-n") {
			numberofpodmodes = std::stoi(arg.substr(2));
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
	if (inputparameter == "" || outputfile == "" || inputtraining == "" || inputmodel == "" || inputpod == "") {
		cout << "Usage:\n -idesignparameterfile\n -ttrainingdatapath\n -mmodeldatapath\n -ppodpath\n -ooutputpath\n [-n(Number of Modes, default 20)]\n [-r(Regression Parameter)]\n No spaces between option and content permitted, e.g. -n20\n Aborting.";
		return -1;
	}

	// Evaluate kriging model for modes and construct outputvector using pod modes
	cout << "Loading POD data\n";
	mat podBasis;
	podBasis.load(inputpod + "podModes.mat");

	// We subtract the snapshot mean before doing pod, so we need to add it here again
	rowvec meanvalue;
    meanvalue.load(inputpod + "snapshotMeanValue.vec");
	vec out = conv_to<vec>::from(meanvalue);

	vec input;
	input.load(inputparameter, raw_ascii);
	//input = zeros(38);

	mat trainingsdata;
	vec krigingparams;

	// process kriging models to obtain mode coefficients, add pod modevector times coefficient to output 
	for (int i = 0; i < numberofpodmodes; i++) {
		cout << "Loading Kriging Model for Mode: " << i << "\n";
		trainingsdata.load(inputtraining + "TrainingDataMode" + std::to_string(i) + ".csv", raw_ascii);
		krigingparams.load(inputmodel + "KriginParams" + std::to_string(i) + ".vec");


		cout << "Evaluating Kriging Model for Mode: " << i << "\n";
		double modecoefficient = evaluateKrigingModel(trainingsdata, krigingparams, regparam, input);

		cout << "Updating output\n";
		out += podBasis.col(i) * modecoefficient;
	}
	cout << "Writing outputfile\n";

	out.save(outputfile, raw_ascii);

	return 0;
}

double evaluateKrigingModel(mat trainingsdata, vec krigingweights, double reg_param, vec input)
{
	// The prediction method in RoDeO (calculate_f_tilde) needs some input we need to calculate first
	// Method following implementation from RoDeO

	int nrows = trainingsdata.n_rows;
	int ncols = trainingsdata.n_cols;
	int dim = ncols - 1;

	int dimension_of_R = nrows;

	double beta0 = 0.0;

	mat X = trainingsdata.submat(0, 0, nrows - 1, ncols - 2);

	vec ys = trainingsdata.col(dim);

	vec x_max(dim);
	x_max.fill(0.0);

	vec x_min(dim);
	x_min.fill(0.0);

	for (int i = 0; i < dim; i++) {
		x_max(i) = trainingsdata.col(i).max();
		x_min(i) = trainingsdata.col(i).min();

	}

	/* normalize data */
	for (int i = 0; i < nrows; i++) {
		for (int j = 0; j < dim; j++) {
			X(i, j) = (X(i, j) - x_min(j)) / (x_max(j) - x_min(j));
		}
	}

	vec normalizedInput;
	normalizedInput = zeros(dim);

	for (int i = 0; i < dim; i++) {
		normalizedInput(i) = (input(i) - x_min(i)) / (x_max(i) - x_min(i));
	}

	// Removed some code for linear regression here, as we don't use it on the training either

	/* allocate the correlation matrix */
	mat R = zeros(dimension_of_R, dimension_of_R);

	// Set kriging modell data
	vec regression_weights;
	// calculate_f_tile should accept an linear_regression on/off switch, but checks the length of regression_weights instead. 
	// So we need to set it to a vec with length 0, as we disabled lin reg on the training. Bäh. Should be fixed in RoDeO
	regression_weights = zeros(0);

	vec theta = krigingweights.head(dim);
	vec gamma = krigingweights.tail(dim);

	reg_param = pow(10.0, -1.0*reg_param);

	/* evaluate the correlation matrix for the given theta and gamma */
	compute_R_matrix(theta, gamma, reg_param, R, X);

	vec I = ones(dimension_of_R);
	vec R_inv_ys(dimension_of_R); /* R^-1 * ys */
	vec R_inv_I(dimension_of_R); /* R^-1 * I */
	vec R_inv_ys_min_beta(dimension_of_R); /* R^-1 * (ys-beta*I) */


	mat Rinv = inv(R);
	beta0 = (1.0 / dot(I, Rinv*I)) * (dot(I, Rinv*ys));
	R_inv_ys_min_beta = Rinv * (ys - beta0 * I);

	// finally we have everything to call calculate_f_tilde

	return  calculate_f_tilde(conv_to<rowvec>::from(normalizedInput),
		X,
		beta0,
		regression_weights,
		R_inv_ys_min_beta,
		krigingweights);
}