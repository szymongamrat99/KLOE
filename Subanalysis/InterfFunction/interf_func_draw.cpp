#include <TCanvas.h>
#include <TF1.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TMath.h>
#include <TLegend.h>
#include <interf_function.h>

#include <const.h>

int main()
{
    KLOE::setGlobalStyle();

    TCanvas *canva = new TCanvas("", "", 750, 750);
    canva->SetMargin(0.15, 0.15, 0.15, 0.1);
    TF1 *func_with_im = new TF1("Interference function2", &interf_function, -20.0, 20.0, 2);

    canva->SetGrid(1, 1);

    gStyle->SetGridStyle(3);
    gStyle->SetGridWidth(2);

    func_with_im->SetNpx(1E4);
    func_with_im->SetLineColor(kBlack);
    func_with_im->SetLineWidth(3);

    func_with_im->SetTitle("");

    func_with_im->GetYaxis()->SetRangeUser(0, 0.3);
    func_with_im->GetYaxis()->SetTitle("I [-]");

    func_with_im->GetXaxis()->SetTitle("#Deltat [#tau_{S}]");

    func_with_im->Draw();

    // TLegend *legend = new TLegend(0.6, 0.73, 0.9, 0.9,"");

    // TString name_0, name_1;

    // const Int_t expRe = floor(log10(abs(0.005))), expIm = floor(log10(abs(0.05)));
    // const Double_t frontRe = 0.005 / pow(10, expRe), frontIm = 0.05 / pow(10, expIm);  

    // name_0 = "#varepsilon'/#varepsilon = 0";
    // name_1 = Form("#splitline{Re(#varepsilon'/#varepsilon) = %.2f#times10^{%d}}{Im(#varepsilon'/#varepsilon) = %.2f#times10^{%d}}",frontRe,expRe,frontIm,expIm);

    // legend->AddEntry(func_wo_im, name_0, "l");
    // legend->AddEntry(func_with_im, name_1, "l");

    // legend->SetTextSize(0.025);
    // legend->SetTextFont(42);

    // legend->Draw();

    canva->Print("interf_func" + Paths::ext_img);
}