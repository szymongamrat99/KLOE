#include <TCanvas.h>
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TFitResult.h>
#include <TPaveText.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <ROOT/RVec.hxx>

#include <neutral_mom.h>
#include <pi0_photon_pair.h>
#include <triple_gaus.h>
#include <interference.h>
#include <kloe_class.h>
#include <lorentz_transf.h>

#include "../inc/omegarec.hpp"

int plots(int first_file, int last_file, int loopcount, int M, int range, Controls::DataType data_type)
{
	std::string path = workdirPath + "/" + "scripts/ROOT_files/DBV-26/MK0/";

	std::string fullname = "",
							dirnamemc = "MONTE_CARLO", dirnamedata = "DATA",
							filenamemc = "mc_stream62_mccard2_", filenamedata = "data_stream42_",
							extension = ".root";
	// std::string TTreeName = path + "/" + dirnamedata + "/" + filenamedata + "*" + extension;

	std::string TTreeName = path + "prod2root_mk0_all_phys2_30401.root";

	// ROOT::EnableImplicitMT();
	ROOT::RDataFrame dfData("h1", TTreeName);

	auto colType = dfData.GetColumnType("XCl");
	// Print column type
	std::cout << "Column " << colType << " has type " << colType << std::endl;

	TCanvas *c1 = new TCanvas("c1", "c1");

	auto dfDataAlias = dfData.Define("NClu", [](int &v) {return v;}, {"nClu"});

	dfDataAlias.GetDefinedColumnNames();

	c1->Print("test.png");

	return 0;
}
