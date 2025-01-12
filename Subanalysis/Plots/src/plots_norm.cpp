#include <TCanvas.h>
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TFitResult.h>
#include <TPaveText.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <ROOT/RDF/HistoModels.hxx>

#include <neutral_mom.h>
#include <pi0_photon_pair.h>
#include <triple_gaus.h>
#include <interference.h>
#include <kloe_class.h>
#include <lorentz_transf.h>

#include "../inc/plots.hpp"

int plotsNorm(int first_file, int last_file, int loopcount, int M, int range, Controls::DataType data_type)
{
	TTree
			*treeMctruth,
			*treeOmega,
			*treeTriangle;

	TFile
			*fileMctruth,
			*fileOmega,
			*fileTriangle;

	TChain *chain = new TChain("INTERF/h1");
	chain_init(chain, first_file, last_file);

	fileMctruth = new TFile(mctruthPath);
	treeMctruth = (TTree *)fileMctruth->Get(gen_vars_tree);

	fileOmega = new TFile(omegaRecPath);
	treeOmega = (TTree *)fileOmega->Get(omegarec_tree);

	fileTriangle = new TFile(trianglePath);
	treeTriangle = (TTree *)fileTriangle->Get(neutrec_triangle_tree);

	treeTriangle->Print();

	chain->AddFriend(treeMctruth, "mc");
	chain->AddFriend(treeOmega, "omega");
	chain->AddFriend(treeTriangle, "triangle");

	ROOT::EnableImplicitMT(5);
	ROOT::RDataFrame dfData(*chain);

	std::vector<std::string>
			alias = {
					"PxOmega",
					"PyOmega",
					"PzOmega",
					"EOmega",
					"TotalMomOmega",
					"MOmega",
					"AngleOmegaPhiLAB",
					"AngleOmegaZAxisLAB",
					"AngleKneZAxisLAB",
					"DistIPAvgIPBoost"};

	std::map<const std::string, Int_t>
			bin_num = {
					{alias[0], 100},
					{alias[1], 100},
					{alias[2], 100},
					{alias[3], 100},
					{alias[4], 100},
					{alias[5], 100},
					{alias[6], 100},
					{alias[7], 100},
					{alias[8], 100},
					{alias[9], 100}};

	std::map<const std::string, Double_t>
			range_min = {
					{alias[0], -500.0},
					{alias[1], -500.0},
					{alias[2], -500.0},
					{alias[3], 600.0},
					{alias[4], 0.0},
					{alias[5], 500.0},
					{alias[6], -10.0},
					{alias[7], -10.0},
					{alias[8], -10.0},
					{alias[9], -0.5}},
			range_max = {
				{alias[0], 500.0}, 
				{alias[1], 500.0}, 
				{alias[2], 500.0}, 
				{alias[3], 1000.0}, 
				{alias[4], 400.0}, 
				{alias[5], 1000.0},
				{alias[6], 190.0},
				{alias[7], 190.0},
				{alias[8], 190.0},
				{alias[9], 2.5}}, 
			bin_width;

	std::map<const std::string, std::string>
			xTitle = {
					{alias[0], "p_{x,#omega} [MeV/c]"},
					{alias[1], "p_{y,#omega} [MeV/c]"},
					{alias[2], "p_{z,#omega} [MeV/c]"},
					{alias[3], "E_{#omega} [MeV]"},
					{alias[4], "p_{tot,#omega} [MeV/c]"},
					{alias[5], "m^{inv}_{#pi^{+}#pi^{-}#pi^{0}} [MeV/c^{2}]"},
					{alias[6], "#angle(#omega,#phi) [#circ]"},
					{alias[7], "#angle(#vec{p}_{#omega},#hat{z}) [#circ]"},
					{alias[8], "#angle(#vec{p}^{triangle}_{K#rightarrow#pi^{0}#pi^{0}},#hat{z}) [#circ]"},
					{alias[9], "|#vec{IP}_{boost} - #vec{IP}_{bhabha}| [cm]"}},
			yTitle,
			flags = {
				{"data", "mcflag == 0"},
				{"mc", "mcflag == 1"},
				{"omegaRec", "doneomega == 1"},
				{"triangleRec", "done_triangle == 1"}
			};

	std::map<const std::string, ROOT::RDF::TH1DModel>
			model;

	for (Int_t i = 0; i < alias.size(); i++)
	{
		bin_width[alias[i]] = abs(range_max[alias[i]] - range_min[alias[i]]) / bin_num[alias[i]];
		yTitle[alias[i]] = (std::string)Form("Counts/%.2f", bin_width[alias[i]]);

		model[alias[i]] = {alias[i].c_str(), "", bin_num[alias[i]], range_min[alias[i]], range_max[alias[i]]};
	}

	std::vector<Int_t>
			mctruth = {1, 3, 4, 5, 6, 7};

	std::vector<std::string>
			filterMctruth(channNum);

	for(Int_t i = 0; i < filterMctruth.size(); i++)
	{
		filterMctruth[i] = "mc.mctruth == " + std::to_string(mctruth[i]) + " && " + 
																									 flags["omegaRec"] + " && " + 
																									 flags["triangleRec"];
	}

	std::vector<TCanvas *> canvas;
	TString canvas_name = "";

	Int_t index;

	auto getArrElement = [&index](ROOT::RVec<Float_t> &v)
	{
		return v[index];
	};

	auto getAngleZAxis = [](ROOT::RVec<Float_t> &v)
	{
		return acos(v[2] / sqrt(pow(v[0], 2) + pow(v[1], 2) + pow(v[2], 2))) * 180. / M_PI;
	};

	auto getAngleVector = [](ROOT::RVec<Float_t> &v, ROOT::RVec<Float_t> &w)
	{
		return acos((v[0] * w[0] + v[1] * w[1] + v[2] * w[2]) / (sqrt(pow(v[0], 2) + pow(v[1], 2) + pow(v[2], 2)) * sqrt(pow(w[0], 2) + pow(w[1], 2) + pow(w[2], 2)))) * 180. / M_PI;
	};

	auto getAngleComponents = [](ROOT::RVec<Float_t> &v, Float_t w1, Float_t w2, Float_t w3)
	{
		return acos((v[0] * w1 + v[1] * w2 + v[2] * w3) / (sqrt(pow(v[0], 2) + pow(v[1], 2) + pow(v[2], 2)) * sqrt(pow(w1, 2) + pow(w2, 2) + pow(w3, 2)))) * 180. / M_PI;
	};

	auto getDistanceComponents = [](ROOT::RVec<Float_t> &v, Float_t w1, Float_t w2, Float_t w3)
	{
		return sqrt(pow(v[0] - w1, 2) + pow(v[1] - w2, 2) + pow(v[2] - w3, 2));
	};

	for (index = 0; index < alias.size(); index++)
	{
		canvas_name = "canvas_" + std::to_string(index);
		canvas.push_back(new TCanvas(canvas_name, canvas_name));

		std::vector<TH1 *> hist;
		TString hist_name = "";

		std::vector<Double_t>
				maxHeight;

		TLegend *legend = new TLegend(0.7,0.7,0.95,0.95);

		for (Int_t j = 0; j < channNum; j++)
		{
			hist_name = "hist_" + std::to_string(j);
			hist.push_back(new TH1D(hist_name, channName[j], bin_num[alias[index]], range_min[alias[index]], range_max[alias[index]]));

			if(index == 6)
			{
				dfData.Define(alias[index].c_str(), getAngleComponents, {"omega.omega", "Bpx", "Bpy", "Bpz"})
					.Filter(filterMctruth[j].c_str())
					.Histo1D(model[alias[index]], alias[index])
					.GetPtr()
					->Copy(*hist[j]);
			}
			else if(index == 7)
			{
				dfData.Define(alias[index].c_str(), getAngleZAxis, {"omega.omega"})
					.Filter(filterMctruth[j].c_str())
					.Histo1D(model[alias[index]], alias[index])
					.GetPtr()
					->Copy(*hist[j]);
			}
			else if(index == 8)
			{
				dfData.Define(alias[index].c_str(), getAngleZAxis, {"triangle.fourKnetriangle"})
					.Filter(filterMctruth[j].c_str())
					.Histo1D(model[alias[index]], alias[index])
					.GetPtr()
					->Copy(*hist[j]);
			}
			else if(index == 9)
			{
				dfData.Define(alias[index].c_str(), getDistanceComponents, {"ip", "Bx", "By", "Bz"})
					.Filter(filterMctruth[j].c_str())
					.Histo1D(model[alias[index]], alias[index])
					.GetPtr()
					->Copy(*hist[j]);
			}
			else
			{
				dfData.Define(alias[index].c_str(), getArrElement, {"omega.omega"})
					.Filter(filterMctruth[j].c_str())
					.Histo1D(model[alias[index]], alias[index])
					.GetPtr()
					->Copy(*hist[j]);
			}

			maxHeight.push_back(hist[j]->GetMaximum());
		};

		std::vector<Int_t>
				indices(maxHeight.size());

		TMath::Sort((Int_t)maxHeight.size(), maxHeight.data(), indices.data(), kTRUE);

		std::map<const std::string, std::string>
				imgNames;

		canvas[index]->cd();

		for (Int_t j = 0; j < channNum; j++)
		{
			hist[indices[j]]->SetLineColor(channColor[indices[j]]);
			legend->AddEntry(hist[j], channName[j], "l");

			if (j == 0)
			{
				if (maxHeight[indices[0]] / maxHeight[indices[channNum - 1]] > 10.)
				{
					canvas[index]->SetLogy(1);
					hist[indices[j]]->GetYaxis()->SetRangeUser(1E-1, maxHeight[indices[j]] * 10.0);
				}

				hist[indices[j]]->SetStats(0);
				hist[indices[j]]->GetYaxis()->SetMaxDigits(3);
				hist[indices[j]]->GetXaxis()->SetTitle(xTitle[alias[index]].c_str());
				hist[indices[j]]->GetYaxis()->SetTitle(yTitle[alias[index]].c_str());
				hist[indices[j]]->Draw();
			}
			else
			{
				hist[indices[j]]->Draw("SAME");
			}
		}

		legend->Draw();

		imgNames[alias[index]] = plots_dir + img_dir + alias[index] + ext_img;

		canvas[index]->Print(imgNames[alias[index]].c_str());
	}

	return 0;
}
