#include <iostream>

#include <TROOT.h>
#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>
#include <TMath.h>

#include "../inc/cpfit.hpp"

int cp_fit_func(KLOE::interference &event, std::vector<std::vector<Double_t>> &relativeErr, std::vector<std::vector<Double_t>> &real, std::vector<std::vector<Double_t>> &imaginary, ErrorHandling::ErrorLogs &logger, KLOE::pm00 &Obj)
{
	const Double_t
			x_min = Utils::properties["variables"]["CPFit"]["histoResults"]["rangeX"][0],
			x_max = Utils::properties["variables"]["CPFit"]["histoResults"]["rangeX"][1],
			res_deltaT = Utils::properties["variables"]["Resolutions"]["deltaT"];
	const UInt_t
			nbins = 1 + ((x_max - x_min) / res_deltaT);

	Double_t split[3] = {-30.0, 0.0, 30.0};

	Bool_t scanFlag = (Bool_t)Utils::properties["variables"]["CPFit"]["cuts"]["cutScanMode"]["flag"];
	Int_t numberOfPoints;
	Double_t
			cutLimits[2] = {0.0},
			cutStep = 0.0;

	std::vector<Double_t>
			cutLimit(numberOfPoints);
	std::vector<std::vector<Double_t>>
			errValue(2),
			realValue(2),
			imaginaryValue(2);
	// ---------------------------------------------------------

	ROOT::Math::Minimizer *minimum =
			ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");

	// set tolerance , etc...
	minimum->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
	minimum->SetTolerance(0.1);
	minimum->SetPrintLevel(0);
	minimum->SetStrategy(2);

	const UInt_t num_of_vars = 11;

	ROOT::Math::Functor minimized_function(&event, &KLOE::interference::interf_chi2, num_of_vars);

	minimum->SetFunction(minimized_function);

	const Double_t init_vars[num_of_vars] = {
			PhysicsConstants::Re,
			PhysicsConstants::Im_nonCPT,
			Utils::properties["variables"]["CPFit"]["initParams"]["Norm"]["Signal"],
			Utils::properties["variables"]["CPFit"]["initParams"]["Norm"]["Regeneration"]["FarLeft"],
			Utils::properties["variables"]["CPFit"]["initParams"]["Norm"]["Regeneration"]["CloseLeft"],
			Utils::properties["variables"]["CPFit"]["initParams"]["Norm"]["Regeneration"]["CloseRight"],
			Utils::properties["variables"]["CPFit"]["initParams"]["Norm"]["Regeneration"]["FarRight"],
			Utils::properties["variables"]["CPFit"]["initParams"]["Norm"]["Omegapi0"],
			Utils::properties["variables"]["CPFit"]["initParams"]["Norm"]["Threepi0"],
			Utils::properties["variables"]["CPFit"]["initParams"]["Norm"]["Semileptonic"],
			Utils::properties["variables"]["CPFit"]["initParams"]["Norm"]["Other"]},
								 step[num_of_vars] = {1E-5, 1E-5, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05};

	Double_t
			limit_upper = Utils::properties["variables"]["CPFit"]["limitPer"]["upper"],
			limit_lower = Utils::properties["variables"]["CPFit"]["limitPer"]["lower"];

	minimum->SetVariable(0, "Real part", init_vars[0], step[0]);
	minimum->SetVariable(1, "Imaginary part", init_vars[1], step[1]);
	minimum->SetLimitedVariable(2, "Norm signal", init_vars[2], step[2], init_vars[2] - 1.0 * init_vars[2], init_vars[2] + 5.0 * init_vars[2]);

	if (x_min > -30.0 && x_max < 30.0)
	{
		minimum->SetFixedVariable(3, "Norm left DC wall", init_vars[3]);
		minimum->SetLowerLimitedVariable(4, "Norm left beam pipe", init_vars[4], step[4], 0.0);
		minimum->SetLowerLimitedVariable(5, "Norm right beam pipe", init_vars[5], step[5], 0.0);
		minimum->SetFixedVariable(6, "Norm right DC wall", init_vars[6]);
	}
	else
	{
		minimum->SetLowerLimitedVariable(3, "Norm left DC wall", init_vars[3], step[3], 0.0);
		minimum->SetLowerLimitedVariable(4, "Norm left beam pipe", init_vars[4], step[4], 0.0);
		minimum->SetLowerLimitedVariable(5, "Norm right beam pipe", init_vars[5], step[5], 0.0);
		minimum->SetLowerLimitedVariable(6, "Norm right DC wall", init_vars[6], step[6], 0.0);
	}

	minimum->SetLimitedVariable(7, "Norm omega", init_vars[7], step[7], init_vars[7] - limit_lower * init_vars[7], init_vars[7] + limit_upper * init_vars[7]);
	minimum->SetLimitedVariable(8, "Norm three", init_vars[8], step[8], init_vars[8] - limit_lower * init_vars[8], init_vars[8] + limit_upper * init_vars[8]);
	minimum->SetLimitedVariable(9, "Norm semi", init_vars[9], step[9], init_vars[9] - limit_lower * init_vars[9], init_vars[9] + limit_upper * init_vars[9]);
	minimum->SetLimitedVariable(10, "Norm other bcg", init_vars[10], step[10], init_vars[10] - limit_lower * init_vars[10], init_vars[10] + limit_upper * init_vars[10]);

	minimum->Minimize();

	relativeErr[0].push_back(abs(minimum->Errors()[0] / minimum->X()[0]));
	relativeErr[1].push_back(abs(minimum->Errors()[1] / minimum->X()[1]));

	real[0].push_back(abs(minimum->X()[0]));
	real[1].push_back(abs(minimum->Errors()[0]));

	imaginary[0].push_back(abs(minimum->X()[1]));
	imaginary[1].push_back(abs(minimum->Errors()[1]));

	return 0;
}
