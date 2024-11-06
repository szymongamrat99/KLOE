#include <TCut.h>

#include "../inc/efficiency.hpp"

#include <const.h>

int efficiencyScan(UInt_t first_file, UInt_t last_file)
{
	ErrorHandling::ErrorLogs logger;

	TFile
			*file_mctruth,
			*file_cutvars,
			*cutvars_csv;

	TTree
			*tree_mctruth,
			*tree_cutvars;

	TString
			mctruth_name = gen_vars_dir + root_files_dir + mctruth_filename + first_file + "_" + last_file + ext_root,
			cutvars_csv_name = "CutVars";

	// Wrapper to read the CSV file properly
	std::ifstream CSVfile(efficiency_dir + input_dir + cutvars_csv_name + ext_csv);

	std::string
			line;

	Int_t
			lineNum = 0,
			headerSize = 0;

	TObjArray
			*tokenize;

	std::vector<TString>
			lineStr,
			header;

	while (getline(CSVfile, line))
	{
		lineNum++;
		lineStr.push_back(line);
	}

	CSVfile.close();

	std::vector<std::map<TString, TString>> cutDetails(lineNum - 1);

	for (Int_t i = 0; i < lineNum; i++)
	{
		tokenize = lineStr[i].Tokenize(";");

		if (i == 0)
		{
			// To get the header
			headerSize = tokenize->GetEntries();
			for (Int_t j = 0; j < headerSize; j++)
			{
				header.push_back(((TObjString *)tokenize->At(j))->String());
			}
		}
		else
		{
			// For other lines
			for (Int_t j = 0; j < headerSize; j++)
			{
				cutDetails[i - 1][header[j]] = ((TObjString *)tokenize->At(j))->String();
			}
		}
	}
	//

	TChain *chain = new TChain("INTERF/h1");
	chain_init(chain, first_file, last_file);

	file_mctruth = new TFile(mctruth_name);
	tree_mctruth = (TTree *)file_mctruth->Get(gen_vars_tree);

	chain->AddFriend(tree_mctruth);

	UInt_t
			sel_ev[channNum],
			total_ev[channNum];

	TCut
			cut,
			channelChoice[channNum];

	TString
			cutStr;

	std::vector<UInt_t>
			points_num;

	std::vector<Float_t>
			max_lim,
			min_lim,
			step;

	std::vector<std::vector<Float_t>>
			x_val[channNum],
			eff[channNum],
			purity(lineNum - 1);

	for (Int_t i = 0; i < channNum; i++)
	{
		x_val[i].resize(lineNum - 1);
		eff[i].resize(lineNum - 1);
	}

	std::vector<TString>
			varName,
			benchmark,
			sign;

	std::vector<Bool_t>
			absVal;

	std::vector<TGraph *>
			eff_graphs[channNum],
			purity_graph;

	for (Int_t i = 0; i < lineNum - 1; i++)
	{
		for (Int_t j = 0; j < headerSize; j++)
		{
			if (header[j] == "VariableName")
				varName.push_back(cutDetails[i][header[j]]);
			if (header[j] == "Benchmark")
				benchmark.push_back(cutDetails[i][header[j]]);
			if (header[j] == "AbsoluteValue")
				absVal.push_back((Bool_t)cutDetails[i][header[j]]);
			if (header[j] == "Sign")
				sign.push_back(cutDetails[i][header[j]]);
			if (header[j] == "LowerLimit")
				min_lim.push_back(std::atof(cutDetails[i][header[j]]));
			if (header[j] == "UpperLimit")
				max_lim.push_back(std::atof(cutDetails[i][header[j]]));
			if (header[j] == "NumberOfPoints")
				points_num.push_back(std::atoi(cutDetails[i][header[j]]));
		}

		step.push_back((max_lim[i] - min_lim[i]) / (Float_t)points_num[i]);
	}

	// Total num of events
	for (Int_t i = 0; i < channNum; i++)
	{
		channelChoice[i] = gen_vars_tree + ".mctruth == " + channelInt[i];
		total_ev[i] = chain->GetEntries(channelChoice[i]);
	}

	for (Int_t k = 0; k < lineNum - 1; k++)
	{
		for (Int_t j = 0; j <= points_num[k]; j++)
		{
			for (Int_t i = 0; i < channNum; i++)
			{
				x_val[i][k].push_back(min_lim[k] + j * step[k]);

				if (absVal[k])
					cutStr = "abs(" + varName[k] + " - " + benchmark[k] + ")" + sign[k];
				else
					cutStr = "(" + varName[k] + " - " + benchmark[k] + ")" + sign[k];

				cut = (cutStr + std::to_string(x_val[i][k][j]) )&& channelChoice[i];

				sel_ev[i] = chain->GetEntries(cut);

				eff[i][k].push_back(sel_ev[i] / (Float_t)total_ev[i]);
			}

			purity[k].push_back(sel_ev[0] / (Float_t)(sel_ev[0] + sel_ev[1] + sel_ev[2] + sel_ev[3] + sel_ev[4] + sel_ev[5]));

			if (purity[k][j] > 1 || std::isnan(purity[k][j]) || std::isinf(purity[k][j]))
				purity[k].push_back(0.);
		}

		auto legend = new TLegend(0.48, 0.1, 0.85, 0.3);

		for (Int_t i = 0; i < channNum; i++)
		{
			eff_graphs[i].push_back(new TGraph(points_num[k] + 1, x_val[i][k].data(), eff[i][k].data()));
			eff_graphs[i][k]->SetTitle(channName[i] + " efficiency");
			eff_graphs[i][k]->SetLineColor(channColor[i]);
			eff_graphs[i][k]->SetLineWidth(3);

			legend->AddEntry(eff_graphs[i][k], channName[i] + " eff");
		}

		purity_graph.push_back(new TGraph(points_num[k] + 1, x_val[0][k].data(), purity[k].data()));
		purity_graph[k]->SetLineColor(kBlack);
		purity_graph[k]->SetTitle(channName[0] + " purity");
		purity_graph[k]->SetLineWidth(3);
		purity_graph[k]->SetLineStyle(9);

		legend->AddEntry(purity_graph[k], channName[0] + " purity");

		TCanvas *canva = new TCanvas();
		canva->SetRightMargin(0.15);
		canva->Range(min_lim[k], 0, max_lim[k], 1);

		TMultiGraph *mg = new TMultiGraph();

		for (Int_t i = 0; i < channNum; i++)
			mg->Add(eff_graphs[i][k]);

		mg->Add(purity_graph[k]);

		TString x_title;

		x_title = cutStr;

		mg->GetXaxis()->SetMaxDigits(3);
		mg->GetXaxis()->SetTitle(x_title);
		mg->GetYaxis()->SetTitle("Efficiency");
		mg->GetXaxis()->CenterTitle(1);
		mg->GetYaxis()->CenterTitle(1);
		mg->GetXaxis()->SetLimits(min_lim[k], max_lim[k]);
		mg->GetYaxis()->SetRangeUser(0, 1);
		mg->Draw("AL");

		TGaxis *axis = new TGaxis(gPad->GetUxmax(), gPad->GetUymin(), gPad->GetUxmax(), gPad->GetUymax(), 0, 1, 110, "+L");
		axis->SetTitle("Purity");
		axis->CenterTitle(1);
		axis->SetVertical(1);
		axis->SetTitleOffset(1.5);
		axis->SetTitleFont(62);
		axis->SetTitleSize(0.05);
		axis->Draw();

		legend->Draw();

		TString pngname = img_dir + "eff_" + cutStr + ".png";
		canva->Print(pngname);
	}

	return 0;
}
