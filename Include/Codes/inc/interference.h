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
#include <array>

namespace KLOE
{
  class interference
  {
  public:
    std::map<TString, std::vector<Double_t>> time_diff;
    std::vector<Double_t> time_diff_gen;
    std::map<TString, std::vector<Double_t>> time_diff_rand_mc, time_diff_rand_data;
    std::vector<Double_t> time_diff_gen_rand_mc, time_diff_gen_rand_data;

    Double_t *exclusions, left_x_split, center_x_split, right_x_split;

    std::vector<Double_t> init_vars, step;

    Double_t *corr_vals, *eff_vals, *resi_vals; 
    
    std::map<TString, std::vector<Double_t>> tmp_norm;

    //////////////////////////////////////////////////////////////////////////////////////////////////////

    interference() {};

    interference(TString mode, Bool_t corr_check, UInt_t bin_number, Double_t x_min, Double_t x_max, Double_t *split) : _mode(mode), _corr_check(corr_check), _bin_number(bin_number), _x_min(x_min), _x_max(x_max)
    {

      if (_mode == "split")
        num_of_vars = 11;
      else if (_mode == "exclude")
        num_of_vars = 8;
      else if (_mode == "mc")
        num_of_vars = 3;
      else if (_mode == "bcg")
        num_of_vars = 8;
      else if (_mode == "final")
        num_of_vars = 3;

      //////////////////////////////////////////////////////////////////////////////////

      left_x_split = split[0];
      center_x_split = split[1];
      right_x_split = split[2];

      for (const auto &name : KLOE::channName)
      {
        e[name.second] = std::vector<Double_t>(bin_number, 0.);
        b[name.second] = std::vector<Double_t>(bin_number, 0.);
      }

      //////////////////////////////////////////////////////////////////////////////////

      resi_vals = new Double_t[bin_number];

      // Channels of MC: pm00, regen, omega, three, semi, other bcg (6)

      for (const auto &name : KLOE::channName)
      {
        _frac[name.second] = new TH1D("Fitted histo " + name.second,
                                      "",
                                      bin_number,
                                      x_min,
                                      x_max);

        _frac[name.second]->SetLineColor(KLOE::channColor.at(name.second));
      }

      // Fractions of MC for 1/2 MC - 1/2 fake DATA fit
      if (_mode == "mc" || _mode == "bcg")
      {
        for (const auto &name : KLOE::channName)
        {
          _frac_data[name.second] = new TH1D("MC 'data' fracs " + name.second, 
                                              "", 
                                              bin_number, 
                                              x_min, 
                                              x_max);

          _frac_data[name.second]->SetLineColor(KLOE::channColor.at(name.second));
        }
      }

     _data_sub = new TH1D("DATA subtraction histogram", "", _bin_number, _x_min, _x_max);

      _mc_sub = new TH1D("MC subtraction histogram", "", _bin_number, _x_min, _x_max);


      for (const auto &name : KLOE::channName)
      {
        b[name.second] = std::vector<Double_t>(_bin_number, 0.);
        e[name.second] = std::vector<Double_t>(_bin_number, 0.);
      }

      //////////////////////////////////////////////////////////////////////////////////////////////////

      corr_vals = new Double_t[_bin_number];

      if (_corr_check == true)
      {
        TFile file("../Efficiency_analysis/correction_factor.root");

        corr_factor = (TGraphErrors *)file.Get("correction_factor");
        corr_vals = corr_factor->GetY();

        file.Close();
      }
      else
      {
        for (Int_t i = 0; i < bin_number; i++)
          corr_vals[i] = 1.;
      }

      ///////////////////////////////////////////////////////////////////////////////////
      ///////////////////////////////////////////////////////////////////////////////////

      eff_vals = new Double_t[bin_number];

      ///////////////////////////////////////////////////////////////////////////////////
    };

    //! Overloading of constructor for exclusion method

    //////////////////////////////////////////////////////////////////////////////////////////////////////

    void bin_extraction(TString channel, TH1 *histogram);

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

    TH1 *getFracHistogram(TString channel)
    {
      return _frac[channel];
    };

    Int_t getBinNumber()
    {
      return _bin_number;
    };

    Bool_t channOmit(TString channel)
    {
      if (ToLower(channel) == "data" || ToLower(channel) == "mc sum" || ToLower(channel) == "pi+pi-pi+pi-")
        return true;
      else
        return false;
    };

  private:
    TString _mode; //! "split", "window", "excluded, "mc", "bcg", "all"
    TGraphErrors *corr_factor, *eff_factor;

    Int_t _bin_number;
    Double_t _x_min, _x_max;

    /////////////////////////////////////////////////////////////////////////////////////////////////////

    std::map<TString, std::vector<Double_t>> b, e;

    std::map<TString, TH1 *> _frac, _frac_data;

    TH1 *_data_sub, *_mc_sub;

    std::vector<Double_t> fit_pars;

    UInt_t num_of_vars;

    Bool_t _corr_check;

    /////////////////////////////////////////////////////////////////////////////////////////////////////
  };
} // namespace KLOE

#endif