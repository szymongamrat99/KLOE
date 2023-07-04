#include "../../Include/const.h"

void efficiency_plots( const UInt_t points_num = 1, TString var = "", TString sign = "", UInt_t unit_ind = 0,
                       Float_t min_lim = 0, Float_t max_lim = 0, TCut additional = "" )
{
    TString fullname = "",
            dirnamemc = "230623_mc", dirnamedata = "230623_data",
            filenamemc = "mc_stream62_mccard2_", filenamedata = "data_stream42_",
            extension = ".root";

    TChain chain("INTERF/h1");

    for(Int_t i = 1; i <= 50; i++)
    {
	if( i != 100 && i != 64 && i != 529 )	
	{	
        	fullname = "../../ROOT_files/" + dirnamemc + "/backup/" + filenamemc + i + extension;
        	chain.Add(fullname);
	}
    }

    UInt_t sel_ev[chann_num - 2], total_ev[chann_num - 2];
    TCut cut;

    //Total num of events
    total_ev[0] = chain.GetEntries("mctruth == 1");
    total_ev[1] = chain.GetEntries("mctruth == 3");
    total_ev[2] = chain.GetEntries("mctruth == 4");
    total_ev[3] = chain.GetEntries("mctruth == 5");
    total_ev[4] = chain.GetEntries("mctruth == 6");
    total_ev[5] = chain.GetEntries("mctruth_pipi == 1");
    total_ev[6] = chain.GetEntries("mctruth == 7");

    Float_t step, x_val[chann_num - 2][points_num + 1], eff[chann_num - 2][points_num + 1], purity[points_num + 1];

    step = (max_lim - min_lim)/(Float_t)points_num;

    for(Int_t j = 0; j <= points_num; j++)
    {
        for(Int_t i = 0; i < chann_num - 2; i++)
        {
            x_val[i][j] = min_lim + j*step;

            cut = var + sign + x_val[i][j];

            if(i == 0) sel_ev[i] = chain.GetEntries(cut && "mctruth == 1" && additional);
            else if(i == 1) sel_ev[i] = chain.GetEntries(cut && "mctruth == 3" && additional);
            else if(i == 2) sel_ev[i] = chain.GetEntries(cut && "mctruth == 4" && additional);
            else if(i == 3) sel_ev[i] = chain.GetEntries(cut && "mctruth == 5" && additional);
            else if(i == 4) sel_ev[i] = chain.GetEntries(cut && "mctruth == 6" && additional);
	        else if(i == 5) sel_ev[i] = chain.GetEntries(cut && "mctruth_pipi == 1" && additional);
            else if(i == 6) sel_ev[i] = chain.GetEntries(cut && "mctruth == 7" && additional);

            eff[i][j] = sel_ev[i]/(Float_t)total_ev[i];
        }

        purity[j] = sel_ev[5]/(Float_t)(sel_ev[0] + sel_ev[1] + sel_ev[2] + sel_ev[3] + sel_ev[4] + sel_ev[5] + sel_ev[6]);
	if(purity[j] > 1 || std::isnan(purity[j]) || std::isinf(purity[j]) ) purity[j] = 1.;
    }

    auto legend = new TLegend(0.48,0.7,0.85,0.9);

    TGraph *eff_graphs[chann_num - 2];
    TGraph *purity_graph;

    for(Int_t i = 0; i < chann_num - 2; i++)
    {
        eff_graphs[i] = new TGraph(points_num + 1, x_val[i], eff[i]);
        eff_graphs[i]->SetTitle(chann_name[i] + " efficiency");
        eff_graphs[i]->SetLineColor(chann_color[i]);
        eff_graphs[i]->SetLineWidth(3);

        legend->AddEntry(eff_graphs[i], chann_name[i] + " eff");
    }

    purity_graph = new TGraph(points_num + 1, x_val[0], purity);
    purity_graph->SetLineColor(kBlack);
    purity_graph->SetTitle(chann_name[4] + " purity");
    purity_graph->SetLineWidth(3);
    purity_graph->SetLineStyle(9);

    legend->AddEntry(purity_graph, chann_name[4] + " purity");

    TCanvas *canva = new TCanvas();
    canva->SetRightMargin(0.15);
    canva->Range(min_lim,0,max_lim,1);

    TMultiGraph *mg = new TMultiGraph();

    for(Int_t i = 0; i < chann_num - 2; i++) mg->Add(eff_graphs[i]);

    mg->Add(purity_graph);

    TString x_title;

    if(sign == ">") x_title = "|#alpha^{CM}_{#pi^{+},#pi^{-}} - 141#circ| more than [#circ]";
    if(sign == "<") x_title = "m^{inv}_{#pi^{+}#pi^{-},1} cut window half-width [MeV/c^{2}]";

    mg->GetXaxis()->SetMaxDigits(3);
    mg->GetXaxis()->SetTitle(x_title);
    mg->GetYaxis()->SetTitle("Efficiency");
    mg->GetXaxis()->CenterTitle(1);
    mg->GetYaxis()->CenterTitle(1);
    mg->GetXaxis()->SetLimits(min_lim,max_lim);
    mg->GetYaxis()->SetRangeUser(0,1);
    mg->Draw("AL");

    TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(), gPad->GetUxmax(), gPad->GetUymax(),0,1,110,"+L");
    axis->SetTitle("Purity");
    axis->CenterTitle(1);
    axis->SetVertical(1);
    axis->SetTitleOffset(1.5);
    axis->SetTitleFont(62);
    axis->SetTitleSize(0.05);
    axis->Draw();

    legend->Draw();

    TString pngname = "eff_" + var +".png";
    canva->Print(pngname);
    

}
