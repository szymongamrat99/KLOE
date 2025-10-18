#ifndef INTERFERENCE_H
#define INTERFERENCE_H

#include "kloe_class.h"
#include <const.h>
#include "TMath.h"
#include "TString.h"
#include "TH1.h"
#include "TGraph.h"
#include <string>
#include <iostream>
#include <TGraphErrors.h>
#include <TFile.h>
#include <vector>

namespace KLOE
{
  class interference : public pm00
  {
    public:

      std::vector<Double_t> time_diff[KLOE::channNum], time_diff_gen, time_diff_data;
      std::vector<Double_t> time_diff_rand_mc[KLOE::channNum], time_diff_rand_data[KLOE::channNum], time_diff_gen_rand_mc, time_diff_gen_rand_data;

      Double_t *exclusions, left_x_split, center_x_split, right_x_split;

      std::vector<Double_t> init_vars, step;

      Double_t *corr_vals, *eff_vals, *resi_vals, tmp_norm[8]; 

      //////////////////////////////////////////////////////////////////////////////////////////////////////

      interference() {};

      interference(TString mode_init, Bool_t corr_check, UInt_t bin_num_init, Double_t x_min_init, Double_t x_max_init, Double_t *split) : mode(mode_init), pm00()
      {

        if (mode == "split") num_of_vars = 11;
        else if (mode == "exclude") num_of_vars = 8;
        else if (mode == "mc") num_of_vars = 3;
        else if (mode == "bcg") num_of_vars = 8;
        else if (mode == "final") num_of_vars = 3;

        /*for(Int_t i = 0; i < num_of_vars; i++)
        {
          if(i == 0) init_vars.push_back(PhysicsConstants::Re);
          if(i == 1) init_vars.push_back(M_PI*PhysicsConstants::Im_nonCPT/180.);
          else init_vars.push_back(1.0);
          
          step.push_back(abs(init_vars[i]/10.));
        }*/

        pm00();

        bin_number = bin_num_init;
        x_min = x_min_init;
        x_max = x_max_init;

        //////////////////////////////////////////////////////////////////////////////////

        left_x_split = split[0];
        center_x_split = split[1];
        right_x_split = split[2];

        //////////////////////////////////////////////////////////////////////////////////

        resi_vals = new Double_t[bin_number];

        // Channels of MC: pm00, regen, omega, three, semi, other bcg (6)

        for (Int_t i = 0; i < KLOE::channNum; i++)
          frac.push_back(new TH1D(("Fitted histo " + std::to_string(i)).c_str(), "", bin_number, x_min, x_max));

        // MC sum histogram

        mc_sum = new TH1D("MC sum", "", bin_number, x_min, x_max);

        // Fractions of MC for 1/2 MC - 1/2 fake DATA fit

        if(mode == "mc" || mode == "bcg")
        {
          for (Int_t i = 0; i < KLOE::channNum; i++)
            frac_data.push_back(new TH1D(("MC 'data' fracs " + std::to_string(i)).c_str(), "", bin_number, x_min, x_max));
        }

        // Total histogram for Data

        data = new TH1D("DATA histogram", "", bin_number, x_min, x_max);

        data_sub = new TH1D("DATA subtraction histogram", "", bin_number, x_min, x_max);

        mc_sub = new TH1D("MC subtraction histogram", "", bin_number, x_min, x_max);

        for (Int_t i = 0; i < bin_number; i++)
		    {
			    b_mcsum.push_back(0.);
			    e_mcsum.push_back(0.);

          b_data.push_back(0.);
			    e_data.push_back(0.);
		    }

        for (Int_t i = 0; i < KLOE::channNum; i++)
          for (Int_t j = 0; j < bin_number; j++)
          {
            b[i].push_back(0.);
            e[i].push_back(0.);
          }


        //////////////////////////////////////////////////////////////////////////////////////////////////

        corr_vals = new Double_t[bin_number];

        if(corr_check)
        {
          TFile file("../Efficiency_analysis/correction_factor.root");

          corr_factor = (TGraphErrors *)file.Get("correction_factor");
          corr_vals = corr_factor->GetY();

          file.Close();
        }
        else
        {
          for(Int_t i = 0; i < bin_number; i++) corr_vals[i] = 1.;
        }

        ///////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////

        eff_vals = new Double_t[bin_number];

        ///////////////////////////////////////////////////////////////////////////////////
      };

      //! Overloading of constructor for exclusion method

      //////////////////////////////////////////////////////////////////////////////////////////////////////

      void bin_extraction(UInt_t channel, TH1 *histogram);

      void mc_randomize();

      //////////////////////////////////////////////////////////////////////////////////////////////////////

      Double_t interf_function(const Float_t x, Int_t check = 1, const Double_t *par = 0);

      Double_t interf_chi2_split(const Double_t *xx);
      Double_t interf_chi2_excluded(const Double_t *xx);

      Double_t interf_chi2_mc(const Double_t *xx);
      Double_t interf_chi2_mc_data(const Double_t *xx);
      Double_t interf_chi2_bcg(const Double_t *xx);

      Double_t interf_chi2(const Double_t *xx);

      void ClearVectors();

      ////////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////////

    private:

      TString mode; //! "split", "window", "excluded, "mc", "bcg", "all"
      TGraphErrors *corr_factor, *eff_factor;

      /////////////////////////////////////////////////////////////////////////////////////////////////////

      std::vector<Double_t> b[6], e[6];
      std::vector<Double_t> b_mcsum, e_mcsum;
      std::vector<Double_t> b_data, e_data;

      std::vector<Double_t> fit_pars;

      UInt_t num_of_vars;

      /////////////////////////////////////////////////////////////////////////////////////////////////////

  };
} // namespace KLOE

#endif