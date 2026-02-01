#pragma once
#include <const.h>
#include <TF2.h>
#include <TF1.h>

namespace KLOE
{
  class BttfInterf
  {
  private:
    TF2
        *func_00pm,
        *func_pm00,
        *func_pmpm;
    TF1
        *func_00pm_1D,
        *func_pm00_1D,
        *func_pmpm_1D_RA,
        *func_pmpm_1D_RB;

    Double_t
        TpmpmRA,
        TpmpmRB,
        T00pm,
        Tpm00;

  public:
    BttfInterf(Double_t func_00pm_int_limit, Double_t func_pm00_int_limit, Double_t func_pmpm_int_limit_RA, Double_t func_pmpm_int_limit_RB) : TpmpmRA(func_pmpm_int_limit_RA), TpmpmRB(func_pmpm_int_limit_RB), T00pm(func_00pm_int_limit), Tpm00(func_pm00_int_limit)
    {
      func_00pm = new TF2("I(#pi^{0}#pi^{0},t_{1},#pi^{+}#pi^{-},t_{2});t_{1} [#tau_{S}]; t_{2} [#tau_{S}]", &interf_function_00pm, 0, 300.0, 0, 300.0, 2);
      func_pm00 = new TF2("I(#pi^{+}#pi^{-},t_{1},#pi^{0}#pi^{0},t_{2});t_{1} [#tau_{S}]; t_{2} [#tau_{S}]", &interf_function_pm00, 0, 300.0, 0, 300.0, 2);
      func_pmpm = new TF2("I(#pi^{+}#pi^{-},t_{1},#pi^{+}#pi^{-},t_{2});t_{1} [#tau_{S}]; t_{2} [#tau_{S}]", &interf_function_pmpm, 0, 300.0, 0, 300.0, 2);
    }
  };
}