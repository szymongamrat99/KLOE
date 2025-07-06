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

	KLOE::pm00 event;

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
					"AngleKnePhiLAB",
					"AngleKneZAxisLAB",
					"DistIPAvgIPBoost",
					"TimeDifference"
					};

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
					{alias[9], 100},
					{alias[10], 100},
					{alias[11], 100}};

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
					{alias[9], -10.0},
					{alias[10], -0.5},
					{alias[11], -90.0}},
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
				{alias[9], 190.0},
				{alias[10], 2.5},
				{alias[11], 90.0}}, 
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
					{alias[8], "#angle(#vec{p}^{triangle}_{K#rightarrow#pi^{0}#pi^{0}},#phi) [#circ]"},
					{alias[9], "#angle(#vec{p}^{triangle}_{K#rightarrow#pi^{0}#pi^{0}},#hat{z}) [#circ]"},
					{alias[10], "|#vec{IP}_{boost} - #vec{IP}_{bhabha}| [cm]"},
					{alias[11], "#Deltat [#tau_{S}]"}},
			yTitle,
			flags = {
				{"data", "mcflag == 0"},
				{"mc", "mcflag == 1"},
				{"omegaRec", "doneomega == 1"},
				{"triangleRec", "done_triangle == 1"}
			};

	std::map<const std::string, ROOT::RDF::TH1DModel>
			model;

	std::map<const std::string, ROOT::RDF::TH2DModel>
			model2D;

	model2D["CombinedAnglesZ"] = {"CombinedAnglesZ", "", bin_num[alias[7]], range_min[alias[7]], range_max[alias[7]], bin_num[alias[9]], range_min[alias[9]], range_max[alias[9]]};

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

	auto setLorentzVectorMom = [](ROOT::RVec<Float_t> &v)
	{
		TLorentzVector vLorentz(v[0], v[1], v[2], v[3]);

		return vLorentz;
	};

	auto setLorentzVectorPos = [](ROOT::RVec<Float_t> &v, ROOT::RVec<Float_t> &w)
	{
		TLorentzVector vLorentz(v[6] - w[0], v[7] - w[1], v[8] - w[2], cVel * v[9]);

		return vLorentz;
	};

	auto setLorentzVectorPosAlt = [](ROOT::RVec<Float_t> &v, ROOT::RVec<Float_t> &w, Float_t TimeCh)
	{
		TLorentzVector vLorentz(v[6] - w[0], v[7] - w[1], v[8] - w[2], cVel * TimeCh);

		return vLorentz;
	};

	auto getTimeCh = [](ROOT::RVec<Float_t> &v, ROOT::RVec<Float_t> &w)
	{
		Float_t 
			boost = cVel * sqrt(pow(v[0], 2) + pow(v[1], 2) + pow(v[2], 2)) / v[3],
			time = sqrt(pow(v[6] - w[0], 2) + pow(v[7] - w[1], 2) + pow(v[8] - w[2], 2)) / boost;

		return time;
	};

	auto getDeltaT = [&event](TLorentzVector momKch, TLorentzVector posKch, TLorentzVector momKne, TLorentzVector posKne)
	{
		return event.DeltaT(&momKch,&posKch,&momKne,&posKne);
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
				dfData.Define(alias[index].c_str(), getAngleComponents, {"triangle.fourKnetriangle", "Bpx", "Bpy", "Bpz"})
					.Filter(filterMctruth[j].c_str())
					.Histo1D(model[alias[index]], alias[index])
					.GetPtr()
					->Copy(*hist[j]);
			}
			else if(index == 9)
			{
				dfData.Define(alias[index].c_str(), getAngleZAxis, {"triangle.fourKnetriangle"})
					.Filter(filterMctruth[j].c_str())
					.Histo1D(model[alias[index]], alias[index])
					.GetPtr()
					->Copy(*hist[j]);
			}
			else if(index == 10)
			{
				dfData.Define(alias[index].c_str(), getDistanceComponents, {"ip", "Bx", "By", "Bz"})
					.Filter(filterMctruth[j].c_str())
					.Histo1D(model[alias[index]], alias[index])
					.GetPtr()
					->Copy(*hist[j]);
			}
			else if(index == 11)
			{
				dfData.Define("Kne4mom", setLorentzVectorMom, {"fourKnetriangle"})
					.Define("Kne4pos", setLorentzVectorPos, {"fourKnetriangle", "ip"})
					.Define("Kch4mom", setLorentzVectorMom, {"Kchboost"})
					.Define("TimeCh", getTimeCh, {"Kchboost", "ip"})
					.Define("Kch4pos", setLorentzVectorPosAlt, {"Kchboost", "ip", "TimeCh"})
					.Define(alias[index], getDeltaT, {"Kch4mom", "Kch4pos", "Kne4mom", "Kne4pos"})
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

		imgNames[alias[index]] = SystemPath::plots_dir + SystemPath::img_dir + alias[index] + ext_img;

		canvas[index]->Print(imgNames[alias[index]].c_str());
	}

	//

	UInt_t hist2DNum = 1;

	std::vector<TCanvas*> canva2D;
	TString canva2D_name = "";
	
	for (Int_t i = 0; i < channNum * hist2DNum; i++)
	{
		canva2D_name = "canva2D_" + std::to_string(i);
		canva2D.push_back(new TCanvas(canva2D_name, canva2D_name));
	};

	std::vector<TH2*> hist2D;
	TString hist2D_name = "";
	
	for (Int_t i = 0; i < channNum; i++)
	{
		hist2D_name = "hist2D_" + std::to_string(i);
		hist2D.push_back(new TH2D());

		dfData.Define(alias[7].c_str(), getAngleZAxis, {"omega.omega"})
					.Define(alias[9].c_str(), getAngleZAxis, {"triangle.fourKnetriangle"})
					.Filter(filterMctruth[i].c_str())
					.Histo2D(model2D["CombinedAnglesZ"], alias[7], alias[9])
					.GetPtr()
					->Copy(*hist2D[i]);
		
		canva2D[i]->cd();

		hist2D[i]->SetStats(0);

		hist2D[i]->GetXaxis()->SetTitle(xTitle[alias[7]].c_str());
		hist2D[i]->GetYaxis()->SetTitle(xTitle[alias[9]].c_str());

		hist2D[i]->Draw("COLZ");

		std::string hist2DName = (std::string)SystemPath::plots_dir + (std::string)SystemPath::img_dir + "CombinedAnglesZ" + (std::string)hist2D_name + (std::string)ext_img;

		canva2D[i]->Print(hist2DName.c_str());

	};



	return 0;
}
