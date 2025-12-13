#include <TCanvas.h>
#include <TF1.h>
#include <TF2.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TMath.h>
#include <TLegend.h>

#include <interf_function.h>
#include <const.h>

int main(int argc, char* argv[])
{
    Float_t Tpmpm1 = 300,
        Tpmpm2 = 300,
        Tpm00 = 300,
        T00pm = 300;

  if (argc >= 5)
  {
    Tpmpm1 = atof(argv[1]);
    Tpmpm2 = atof(argv[2]);
    Tpm00 = atof(argv[3]);
    T00pm = atof(argv[4]);
  }

  KLOE::setGlobalStyle();

  TF2 *func_00pm = new TF2("I(#pi^{0}#pi^{0},t_{1},#pi^{+}#pi^{-},t_{2});t_{1} [#tau_{S}]; t_{2} [#tau_{S}]", &interf_function_00pm, 0.0, 20, 0.0, 20, 0);

  TF2 *func_pm00 = new TF2("I(#pi^{+}#pi^{-},t_{1},#pi^{0}#pi^{0},t_{2});t_{1} [#tau_{S}]; t_{2} [#tau_{S}]", &interf_function_pm00, 0.0, 20, 0.0, 20, 0);

  TF2 *func_pmpm = new TF2("I(#pi^{+}#pi^{-},t_{1},#pi^{+}#pi^{-},t_{2});t_{1} [#tau_{S}]; t_{2} [#tau_{S}]", &interf_function_pmpm, 0.0, 20, 0.0, 20, 0);

  Double_t norm00pm = 1; // func_00pm->Integral(0.0, 20, 0.0, 20);
  Double_t normpm00 = 1; // func_pm00->Integral(0.0, 20, 0.0, 20);
  Double_t normpmpm = 1; // func_pmpm->Integral(0.0, 20, 0.0, 20);

  TCanvas *c_func_00pm = new TCanvas("c_func_00pm", "Interference function 00pm", 800, 600);
  c_func_00pm->SetLogz(1);
  func_00pm->SetNpx(100);
  func_00pm->SetNpy(100);
  func_00pm->Draw("COLZ");
  c_func_00pm->SaveAs(Paths::img_dir + "interf_func_00pm.png");

  TCanvas *c_func_pm00 = new TCanvas("c_func_pm00", "Interference function pm00", 800, 600);
  c_func_pm00->SetLogz(1);
  func_pm00->SetNpx(100);
  func_pm00->SetNpy(100);
  func_pm00->Draw("COLZ");
  c_func_pm00->SaveAs(Paths::img_dir + "interf_func_pm00.png");

  TCanvas *c_func_pmpm = new TCanvas("c_func_pmpm", "Interference function pmpm", 800, 600);
  c_func_pmpm->SetLogz(1);
  func_pmpm->SetNpx(100);
  func_pmpm->SetNpy(100);
  func_pmpm->Draw("COLZ");
  c_func_pmpm->SaveAs(Paths::img_dir + "interf_func_pmpm.png");

  auto func_00pm_1D = [&](Double_t *x, Double_t *par)
  {
    Double_t current_x = x[0];

    auto temp_f1_lambda = [&](Double_t *t, Double_t *par)
    {
      return func_00pm->Eval(current_x, t[0]) / norm00pm;
    };

    TF1 temp_f1("temp_f1_integral", temp_f1_lambda, 0.0, par[1], 0);

    Double_t integral_result = temp_f1.Integral(0.0, par[1]);

    return integral_result;
  };

  auto func_pm00_1D = [&](Double_t *x, Double_t *par)
  {
    Double_t current_x = x[0];

    auto temp_f1_lambda = [&](Double_t *t, Double_t *par)
    {
      return func_pm00->Eval(current_x, t[0]) / normpm00;
    };

    // Utwórz tymczasowy TF1 z funkcji temp_f1_lambda
    // Uwaga: Nowe obiekty TF1 powinny mieć unikalne nazwy!
    TF1 temp_f1("temp_f1_integral", temp_f1_lambda, 0.0, par[1], 0);

    // Użyj TF1::Integral do obliczenia całki po y (od y_min do y_max)
    Double_t integral_result = temp_f1.Integral(0.0, par[1]);

    return integral_result;
  };

  auto func_pmpm_1D = [&](Double_t *x, Double_t *par)
  {
    Double_t current_x = x[0];

    auto temp_f1_lambda = [&](Double_t *t, Double_t *par)
    {
      // F_xy->GetParameter(i) to wartość p[i]
      // par[0] to wartość x - używamy, gdyby x było w p
      return func_pmpm->Eval(current_x, t[0]) / normpmpm;
    };

    // Utwórz tymczasowy TF1 z funkcji temp_f1_lambda
    // Uwaga: Nowe obiekty TF1 powinny mieć unikalne nazwy!
    TF1 temp_f1("temp_f1_integral", temp_f1_lambda, 0.0, par[0], 0);

    // Użyj TF1::Integral do obliczenia całki po y (od y_min do y_max)
    Double_t integral_result = temp_f1.Integral(0.0, par[0]);

    return integral_result;
  };

  auto func_pmpm_00pm_1D = [&](Double_t *x, Double_t *par)
  {
    return func_pmpm_1D(x, par) / func_00pm_1D(x, par);
  };

  auto func_pmpm_pm00_1D = [&](Double_t *x, Double_t *par)
  {
    return func_pmpm_1D(x, par) / func_pm00_1D(x, par);
  };

  auto func_RA_RB_1D = [&](Double_t *x, Double_t *par)
  {
    return func_pmpm_00pm_1D(x, par) / func_pmpm_pm00_1D(x, par);
  };

  ////////////////////////////////////////////////////////

  TF1 *pm00 = new TF1("pm00", func_pm00_1D, 0, 20, 0);

  pm00->SetParameter(1, Tpm00);

  // 4. Rysowanie wyników
  TCanvas *c_pm00 = new TCanvas("c_pm00", "Integral of TF2 over y", 800, 600);

  pm00->SetNpx(1000);
  pm00->Draw();

  c_pm00->SaveAs(Paths::img_dir + "interf_func_pm001D_draw.png");

  ///////////////////////////////////////////////////////////////////

  TF1 *oopm = new TF1("00pm", func_00pm_1D, 0, 20, 0);
  oopm->SetParameter(1, T00pm);

  // 4. Rysowanie wyników
  TCanvas *c_00pm = new TCanvas("c_00pm", "Integral of TF2 over y", 800, 600);

  oopm->SetNpx(1000);
  oopm->Draw();

  c_00pm->SaveAs(Paths::img_dir + "interf_func_00pm1D_draw.png");

  ////////////////////////////////////////////////////////////////////

  TF1 *pmpm = new TF1("pmpm", func_pmpm_1D, 0, 20, 0);
  pmpm->SetParameter(0, Tpmpm1);

  // 4. Rysowanie wyników
  TCanvas *c_pmpm = new TCanvas("c_pmpm", "Integral of TF2 over y", 800, 600);

  pmpm->SetNpx(1000);
  pmpm->Draw();

  c_pmpm->SaveAs(Paths::img_dir + "interf_func_pmpm1D_draw.png");

  ////////////////////////////////////////////////////////////////////

  Double_t *x = new Double_t[1];
  Double_t *par = new Double_t[2];

  x[0] = 0;

  TF1 *const_line = new TF1("const_line", "1 + 6 * [0]", 0, 20);
  const_line->SetParameter(0, PhysicsConstants::Re);
  const_line->SetLineColor(kBlack);
  const_line->SetLineStyle(2);
  const_line->SetLineWidth(2);

  TF1 *G_x = new TF1("R_{A}(t_{1});t_{1} [#tau_{S}];R_{A}(t_{1}) [-]", func_pmpm_00pm_1D, 0, 20, 0);

  G_x->SetParameter(0, Tpmpm1);
  G_x->SetParameter(1, T00pm);

  // 4. Rysowanie wyników
  TCanvas *c = new TCanvas("c", "Integral of TF2 over y", 800, 600);

  G_x->SetNpx(1000);
  G_x->Draw();
  const_line->Draw("SAME");

  c->SaveAs(Paths::img_dir + "interf_func_RA_draw.png");

  TCanvas *cg = new TCanvas("cg", "Integral of TF2 over y", 800, 600);

  TF1 *F_x = new TF1("R_{B}(t_{1});t_{1} [#tau_{S}];R_{B}(t_{1}) [-]", func_pmpm_pm00_1D, 0, 20, 0);

  F_x->SetParameter(0, Tpmpm2);
  F_x->SetParameter(1, Tpm00);

  F_x->SetNpx(1000);
  F_x->Draw();
  const_line->Draw("SAME");

  cg->SaveAs(Paths::img_dir + "interf_func_RB_draw.png");

  // 4. Rysowanie wyników
  TCanvas *c1 = new TCanvas("c1", "Integral of TF2 over y", 800, 600);
  TF1 *E_x = new TF1("R_{C}(t_{1});t_{1} [#tau_{S}];R_{C}(t_{1}) [-]", func_RA_RB_1D, 0, 20, 0);

  E_x->SetNpx(1000);
  E_x->Draw();

  TF1 *const_line2 = new TF1("const_line2", "1 + 6 * [0]", 0, 20);
  const_line2->SetParameter(0, PhysicsConstants::Re);
  const_line2->SetLineColor(kBlack);
  const_line2->SetLineStyle(2);
  const_line2->SetLineWidth(2);

  TF1 *const_line3 = new TF1("const_line3", "1 - 6 * [0]", 0, 20);
  const_line3->SetParameter(0, PhysicsConstants::Re);
  const_line3->SetLineColor(kBlack);
  const_line3->SetLineStyle(2);
  const_line3->SetLineWidth(2);

  // 4. Rysowanie wyników
  // F_x->SetTitle("F(x) = #int_{-1}^{1} F(x, y) dy");
  const_line2->Draw("SAME");
  const_line3->Draw("SAME");

  c1->SaveAs(Paths::img_dir + "interf_func_RC_double_draw.png");
}