#include <algorithm>
#include <random>

#include <TH1.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>
#include <TRatioPlot.h>

#include <interference.h>

namespace KLOE
{
  Double_t interference::interf_function(const Double_t x, Int_t check, const Double_t *par)
  {
    Double_t Value = 0;
    Double_t Epsilon = 0, RePart = 0, Dphi = 0, TauKs = 0, TauKl = 0, MassDiff = 0;
    Double_t ImPart = 0, GammaKl = 0, GammaKs = 0, Gamma = 0, DMass = 0;

    Double_t dt = x;

    // Parameters from PDG2023

    Epsilon = PhysicsConstants::mod_epsilon;
    Dphi = PhysicsConstants::phi_pm_nonCPT - PhysicsConstants::phi_00_nonCPT; // phi(+-)-phi(00) (degrees)
    TauKs = PhysicsConstants::tau_S_nonCPT;                                   // PDG fit not assuming CPT (s)
    TauKl = PhysicsConstants::tau_L;                                          // Kl mean life (s)
    MassDiff = PhysicsConstants::delta_mass_nonCPT;                           // M(Kl)-M(Ks) ( (h/2pi)s-1 ):

    if (check == 0)
    {
      // Pobieramy indeksy
      Int_t idxRe = fParamIndices["Re"];
      Int_t idxIm = fParamIndices["Im"];

      // Jeśli indeks jest poprawny (nie -1), bierzemy wartość z minimizera.
      // W przeciwnym razie bierzemy 0.0 (lub inną wartość domyślną).
      RePart = (idxRe != -1) ? par[idxRe] : 0.0;
      ImPart = (idxIm != -1) ? par[idxIm] : 0.0;
    }
    else
    {
      RePart = PhysicsConstants::Re;
      ImPart = PhysicsConstants::Im_nonCPT; // M_PI * (Dphi / 3.) / 180.; // Im(epsilon'/epsilon) = Dphi/3;
    }

    // All parameters are calculated taking into account that DT is in TauKs units
    GammaKs = 1.;
    GammaKl = TauKs / TauKl;
    Gamma = GammaKs + GammaKl;
    DMass = MassDiff * TauKs;

    if (dt >= 0.)
    {
      Value = (1. + 2. * RePart) * exp(-GammaKl * dt) +
              (1. - 4. * RePart) * exp(-GammaKs * dt) -
              2. * exp(-0.5 * Gamma * dt) *
                  ((1. - RePart) * cos(DMass * dt) +
                   3. * ImPart * sin(DMass * dt));
    }
    else
    {
      Value = (1. + 2. * RePart) * exp(-GammaKs * std::abs(dt)) +
              (1. - 4. * RePart) * exp(-GammaKl * std::abs(dt)) -
              2. * exp(-0.5 * Gamma * std::abs(dt)) *
                  ((1. - RePart) * cos(DMass * std::abs(dt)) -
                   3. * ImPart * sin(DMass * std::abs(dt)));
    }

    return (std::pow(Epsilon, 2) / (2. * Gamma)) * Value * 100000;
  }

  Double_t interference::double_exponential(const Double_t x, Int_t check, const Double_t *par)
  {
    Double_t tau_S_org = 0.0587;
    Double_t tau_L_org = 49.75;

    Double_t gammaS = 1.0;                   // Wartość gamma_S do ustawienia zakresów
    Double_t gammaL = tau_S_org / tau_L_org; // Wartość gamma_L do ustawienia zakresów

    Double_t dt = x * PhysicsConstants::tau_S_nonCPT / tau_S_org;

    Double_t RePart = 0;

    if (check == 0)
    {
      RePart = par[fParamIndices["Re"]];
    }
    else
    {
      RePart = 0.0; // PhysicsConstants::Re;
    }

    Double_t value = 0.;

    if (dt >= 0)
    {
      value = (1. + 2. * RePart) * exp(-gammaL * dt) +
              (1. - 4. * RePart) * exp(-gammaS * dt);
    }
    else
    {
      value = (1. + 2. * RePart) * exp(-gammaS * std::abs(dt)) +
              (1. - 4. * RePart) * exp(-gammaL * std::abs(dt));
    }

    return value;
  }

  Double_t interference::fit_function(const Double_t x, Int_t check, const Double_t *par)
  {
    return interf_function(x, check, par) / double_exponential(x, 1, par);
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////

  void interference::bin_extraction(TString channel, TH1 *histogram)
  {
    for (Int_t i = 0; i < _bin_number; i++)
      b[channel][i] = (histogram->GetBinContent(i + 1));

    for (Int_t i = 0; i < _bin_number; i++)
      e[channel][i] = (histogram->GetBinError(i + 1));
  };

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////

  void interference::mc_randomize()
  {
    UInt_t rnd_ind;
    srand(time(NULL));

    for (Int_t i = 0; i < KLOE::channNum; i++)
    {
      for (Int_t j = 0; j < time_diff[i].size(); j++)
      {
        rnd_ind = rand() % 2;

        if (rnd_ind == 0)
        {
          if (i == 0)
          {
            time_diff_rand_mc[i].push_back(time_diff[i][j]);
            time_diff_gen_rand_mc.push_back(time_diff_gen[j]);
          }
          else
            time_diff_rand_mc[i].push_back(time_diff[i][j]);
        }
        else if (rnd_ind == 1)
        {
          if (i == 0)
          {
            time_diff_rand_data[i].push_back(time_diff[i][j]);
            time_diff_gen_rand_data.push_back(time_diff_gen[j]);
          }
          else
            time_diff_rand_data[i].push_back(time_diff[i][j]);
        }
      }
    }
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////

  //! Fitting with splitted regeneration into 4 parts
  Double_t interference::interf_chi2_split(const Double_t *xx)
  {
    /////////////////////////////////////////////////////////////////////////////////////////////
    for (const auto &name : KLOE::channName)
    {
      if (ToLower(name.second) == "data")
        goto skip_data;

      if (channOmit(name.second))
        continue;

      for (UInt_t j = 0; j < time_diff[name.second].size(); j++)
      {
        if (ToLower(name.second) == "signal")
          _frac[name.second]->Fill(time_diff[name.second][j], fit_function(time_diff_gen[j], 0, xx)); //! Filling Signal
        else
          _frac[name.second]->Fill(time_diff[name.second][j]); //! Filling background
      }

      //! Using correction factor and efficiency
      if (ToLower(name.second) == "signal")
      {
        _frac[name.second]->Scale(_frac[name.second]->GetEntries() / _frac[name.second]->Integral(0, _bin_number + 1));

        for (UInt_t j = 0; j < _bin_number; j++)
        {
          _frac[name.second]->SetBinContent(j + 1, _frac[name.second]->GetBinContent(j + 1) * corr_vals[j]);
        }

        interference::bin_extraction(name.second, _frac[name.second]);
      }
      else
        interference::bin_extraction(name.second, _frac[name.second]);

      /////////////////////////////////////////////////////////////////////////////////////////////

    skip_data:

      if (ToLower(name.second) != "data")
        continue;

      for (Int_t j = 0; j < time_diff[name.second].size(); j++)
        _frac[name.second]->Fill(time_diff[name.second][j]); //! Filling DATA

      interference::bin_extraction(name.second, _frac[name.second]);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////

    for (const auto &name : KLOE::channName)
      _frac[name.second]->Reset("ICESM");

    /////////////////////////////////////////////////////////////////////////////////////////////

    Double_t value = 0;

    //! Sum of bins and errors for the fitting further

    b["MC sum"].clear();
    e["MC sum"].clear();

    b["MC sum"].resize(_bin_number, 0.0);
    e["MC sum"].resize(_bin_number, 0.0);

    for (const auto &name : KLOE::channName)
    {
      if (channOmit(name.second))
        continue;

      for (Int_t j = 0; j < _bin_number; j++)
      {
        // Pobieramy dt środka binu, aby wiedzieć w którym przedziale splitu jesteśmy
        Double_t dt_bin = _frac[name.second]->GetBinCenter(j + 1);

        // Magia: get_weight sam wie, czy użyć parametru pojedynczego, czy splitu
        Double_t weight = get_weight(name.second, dt_bin, xx);

        b["MC sum"][j] += weight * b[name.second][j];
        e["MC sum"][j] += std::pow(weight * e[name.second][j], 2);
      }
    }

    for (Int_t i = 0; i < _bin_number; i++)
    {
      value += std::pow(b["Data"][i] - b["MC sum"][i], 2) / (std::pow(e["Data"][i], 2) + (e["MC sum"][i]));
    }

    return value;
  };

  //! Fitting with excluded regeneration peaks
  Double_t interference::interf_chi2_excluded(const Double_t *xx)
  {
    // Exclusion of regeneration peaks - to be defined in a class' constructor

    Double_t ReFit = xx[0];
    Double_t ImFit = xx[1];
    Double_t Norm[6] = {xx[2], xx[3], xx[4], xx[5], xx[6], xx[7]};

    /////////////////////////////////////////////////////////////////////////////////////////////
    for (Int_t i = 0; i < KLOE::channNum; i++)
    {
      for (Int_t j = 0; j < time_diff[i].size(); j++)
      {
        if (time_diff[i][j] < exclusions[0] && time_diff[i][j] > exclusions[1] &&
            time_diff[i][j] < exclusions[2] && time_diff[i][j] > exclusions[3])
        {
          if (i == 0)
            _frac[i]->Fill(time_diff[i][j], fit_function(time_diff_gen[j], 0, xx));
          else
            _frac[i]->Fill(time_diff[i][j]);
        }
      }

      if (i == 0)
      {
        _frac[i]->Scale(_frac[i]->GetEntries() / _frac[i]->Integral(0, _bin_number + 1));
        interference::bin_extraction(i, _frac[i]);
      }
      else
        interference::bin_extraction(i, _frac[i]);
    }

    for (Int_t i = 0; i < KLOE::channNum - 2; i++)
      _frac[i]->Reset("ICESM");
    /////////////////////////////////////////////////////////////////////////////////////////////

    Double_t value = 0, *bin_sum, *err_sum;

    bin_sum = new Double_t[_bin_number];
    err_sum = new Double_t[_bin_number];

    for (Int_t i = 0; i < _bin_number; i++)
    {
      bin_sum[i] = 0.;
      err_sum[i] = 0.;
    }

    for (Int_t i = 0; i < KLOE::channNum - 3; i++)
      for (Int_t j = 0; j < _bin_number; j++)
      {
        bin_sum[j] += Norm[i] * b[i][j];
        err_sum[j] += std::pow(Norm[i] * e[i][j], 2);
      }

    for (Int_t i = 0; i < _bin_number; i++)
      value += std::pow(b[6][i] - bin_sum[i], 2) / (std::pow(e[6][i], 2) + err_sum[i] / 2.);

    return value;
  };

  //! Fitting 1/2 signal MC to 1/2 'DATA'
  Double_t interference::interf_chi2_mc(const Double_t *xx)
  {

    //* Fitting interference for MC splitted in half

    Double_t ReFit = xx[0];
    Double_t ImFit = xx[1];
    Double_t Norm = xx[2];

    /////////////////////////////////////////////////////////////////////////////////////////////

    for (Int_t j = 0; j < time_diff_rand_mc["Signal"].size(); j++)
      _frac["Signal"]->Fill(time_diff_rand_mc["Signal"][j], fit_function(time_diff_gen_rand_mc[j], 0, xx));

    for (Int_t j = 0; j < time_diff_rand_data["Signal"].size(); j++)
      _frac["Data"]->Fill(time_diff_rand_data["Signal"][j], fit_function(time_diff_gen_rand_data[j], 0, xx));

    _frac["Signal"]->Scale(_frac["Signal"]->GetEntries() / _frac["Signal"]->Integral(0, _bin_number + 1));
    interference::bin_extraction("Signal", _frac["Signal"]);

    _frac["Data"]->Scale(_frac["Data"]->GetEntries() / _frac["Data"]->Integral(0, _bin_number + 1));
    interference::bin_extraction("Data", _frac["Data"]);

    _frac["Signal"]->Reset("ICESM");
    _frac["Data"]->Reset("ICESM");

    /////////////////////////////////////////////////////////////////////////////////////////////

    Double_t value = 0;

    for (Int_t i = 0; i < _bin_number; i++)
    {
      b["MC sum"][i] = 0.;
      e["MC sum"][i] = 0.;
    }

    for (Int_t j = 0; j < _bin_number; j++)
    {
      b["MC sum"][j] += Norm * b["Signal"][j];
      e["MC sum"][j] += std::pow(Norm * e["Signal"][j], 2);
    }

    for (Int_t i = 0; i < _bin_number; i++)
      value += std::pow(b["Data"][i] - b["MC sum"][i], 2) / (std::pow(e["Data"][i], 2) + e["MC sum"][i]);

    return value;
  };

  //! Fitting signal MC to DATA after efficiency correction
  Double_t interference::interf_chi2_mc_data(const Double_t *xx)
  {
    Double_t ReFit = xx[0];
    Double_t ImFit = xx[1];
    std::map<TString, std::vector<Double_t>> Norm;

    for (const auto &name : KLOE::channName)
    {
      if (name.second == "Signal")
      {
        Norm[name.second].push_back(xx[2]); // Signal norm
      }

      if (name.second == "Regeneration")
      {
        Norm[name.second].push_back(tmp_norm[name.second][0]); // Left DC Wall
        Norm[name.second].push_back(tmp_norm[name.second][1]); // Left beam pipe Wall
        Norm[name.second].push_back(tmp_norm[name.second][2]); // Right beam pipe Wall
        Norm[name.second].push_back(tmp_norm[name.second][3]); // Right DC Wall
      }

      if (name.second == "Omega")
        Norm[name.second].push_back(tmp_norm[name.second][4]); // Background norms

      if (name.second == "3pi0")
        Norm[name.second].push_back(tmp_norm[name.second][5]); // Background norms

      if (name.second == "Semileptonic")
        Norm[name.second].push_back(tmp_norm[name.second][6]); // Background norms

      if (name.second == "Other")
        Norm[name.second].push_back(tmp_norm[name.second][7]); // Background norms
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    for (const auto &name : KLOE::channName)
    {
      for (Int_t j = 0; j < time_diff[name.second].size(); j++)
      {
        if (name.second == "Signal")
          _frac[name.second]->Fill(time_diff[name.second][j], Norm[name.second][0] * interf_function(time_diff[name.second][j], 0, xx)); //! Filling Signal
        else if (name.second == "Regeneration")
        {
          if (time_diff[name.second][j] < left_x_split)
            _frac[name.second]->Fill(time_diff[name.second][j], Norm[name.second][0]);
          else if (time_diff[name.second][j] > left_x_split && time_diff[name.second][j] < center_x_split)
            _frac[name.second]->Fill(time_diff[name.second][j], Norm[name.second][1]);
          else if (time_diff[name.second][j] > center_x_split && time_diff[name.second][j] < right_x_split)
            _frac[name.second]->Fill(time_diff[name.second][j], Norm[name.second][2]);
          else if (time_diff[name.second][j] > right_x_split)
            _frac[name.second]->Fill(time_diff[name.second][j], Norm[name.second][3]);
        }
        else
        {
          _frac[name.second]->Fill(time_diff[name.second][j], Norm[name.second][0]); //! Filling background
        }
      }

      //! Using correction factor and efficiency
      if (name.second == "Signal")
      {
        _frac[name.second]->Scale(_frac[name.second]->GetEntries() / _frac[name.second]->Integral(0, _bin_number + 1));

        for (Int_t k = 1; k <= _bin_number; k++)
          _frac[name.second]->SetBinContent(k, _frac[name.second]->GetBinContent(k) / corr_vals[k - 1]);

        interference::bin_extraction(name.second, _frac[name.second]);
      }

      /////////////////////////////////////////////////////////////////////////////////////////////

      if (name.second != "Data")
        continue;

      for (Int_t j = 0; j < time_diff[name.second].size(); j++)
        _frac[name.second]->Fill(time_diff[name.second][j]); //! Filling DATA

      for (const auto &name1 : KLOE::channName)
      {
        _frac[name.second]->Add(_frac[name1.second], -1.);
      }

      for (Int_t k = 1; k <= _bin_number; k++)
        _frac[name.second]->SetBinContent(k, _frac[name.second]->GetBinContent(k) / corr_vals[k - 1]);

      interference::bin_extraction(6, _frac[name.second]);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////

    for (const auto &name : KLOE::channName)
      _frac[name.second]->Reset("ICESM");

    /////////////////////////////////////////////////////////////////////////////////////////////

    Double_t value = 0;

    //! Sum of bins and errors for the fitting further

    b["MC sum"].resize(_bin_number, 0.0);
    e["MC sum"].resize(_bin_number, 0.0);

    for (Int_t j = 0; j < _bin_number; j++)
    {
      b["MC sum"][j] += Norm["Signal"][0] * b["Signal"][j];
      e["MC sum"][j] += std::pow(Norm["Signal"][0] * e["Signal"][j], 2);
    }

    for (Int_t i = 0; i < _bin_number; i++)
    {
      value += std::pow(b["Data"][i] - b["MC sum"][i], 2) / (std::pow(e["Data"][i], 2) + e["MC sum"][i]);
    }

    return value;
  };

  //! Fitting 1/2 total MC to 1/2 'DATA'
  Double_t interference::interf_chi2_bcg(const Double_t *xx)
  {

    //* Fitting interference for MC splitted in half

    Double_t ReFit = xx[0];
    Double_t ImFit = xx[1];
    Double_t Norm[6] = {xx[2], xx[3], xx[4], xx[5], xx[6], xx[7]};

    /////////////////////////////////////////////////////////////////////////////////////////////

    for (auto const &name : KLOE::channName)
    {
      for (Int_t j = 0; j < time_diff_rand_mc[name.second].size(); j++)
      {
        if (name.second == "Signal")
          _frac[name.second]->Fill(time_diff_rand_mc[name.second][j], interf_function(time_diff_gen_rand_mc[j], 0, xx));
        else
          _frac[name.second]->Fill(time_diff_rand_mc[name.second][j]);
      }

      for (Int_t j = 0; j < time_diff_rand_data[name.second].size(); j++)
      {
        if (name.second == "Signal")
          _frac_data[name.second]->Fill(time_diff_rand_data[name.second][j], interf_function(time_diff_gen_rand_data[j]));
        else
          _frac_data[name.second]->Fill(time_diff_rand_data[name.second][j]);
      }

      _frac[name.second]->Scale(_frac[name.second]->GetEntries() / _frac[name.second]->Integral(0, _bin_number + 1));
      interference::bin_extraction(name.second, _frac[name.second]);

      _frac_data[name.second]->Scale(_frac_data[name.second]->GetEntries() / _frac_data[name.second]->Integral(0, _bin_number + 1));
      _frac["Data"]->Add(_frac_data[name.second]);

      interference::bin_extraction("Data", _frac["Data"]);

      _frac[name.second]->Reset("ICESM");
      _frac_data[name.second]->Reset("ICESM");
    }

    /////////////////////////////////////////////////////////////////////////////////////////////

    Double_t value = 0;

    for (Int_t i = 0; i < _bin_number; i++)
    {
      b["MC sum"][i] = 0.;
      e["MC sum"][i] = 0.;
    }

    for (Int_t i = 0; i < KLOE::channNum; i++)
      for (Int_t j = 0; j < _bin_number; j++)
      {
        b["MC sum"][j] += Norm[i] * b[i][j];
        e["MC sum"][j] += std::pow(Norm[i] * e[i][j], 2);
      }

    for (Int_t i = 0; i < _bin_number; i++)
      value += std::pow(b["Data"][i] - b["MC sum"][i], 2) / (std::pow(e["Data"][i], 2) + e["MC sum"][i]);

    return value;
  };

  //! Auxillary function to use any of the others
  Double_t interference::interf_chi2(const Double_t *xx)
  {
    if (_mode == "split")
      return interf_chi2_split(xx);
    else if (_mode == "excluded")
      return interf_chi2_excluded(xx);
    else if (_mode == "mc")
      return interf_chi2_mc(xx);
    else if (_mode == "bcg")
      return interf_chi2_bcg(xx);
    else if (_mode == "final")
      return interf_chi2_mc_data(xx);
    else
      return -999.0;
  };

  void interference::ClearVectors()
  {
    for (const auto &name : KLOE::channName)
    {
      time_diff[name.second].clear();
      time_diff[name.second].shrink_to_fit();
    }

    time_diff_gen.clear();
    time_diff_gen.shrink_to_fit();
  }

  Double_t KLOE::interference::get_weight(TString channel, Double_t dt, const Double_t *xx)
  {
    auto it = channel_to_indices.find(channel);
    if (it == channel_to_indices.end() || it->second.empty())
      return 0.0;

    const std::vector<Int_t> &indices = it->second;

    // SCENARIUSZ A: Kanał ma tylko jeden parametr (Signal, Omega, lub grupa A_bcg)
    if (indices.size() == 1)
    {
      Int_t idx = indices[0];
      return (idx != -1) ? xx[idx] : 1.0;
    }

    // SCENARIUSZ B: Kanał ma wiele parametrów (Regeneracja ze splitem)
    if (channel == "Regeneration" && indices.size() >= 4)
    {
      if (dt < left_x_split)
        return (indices[0] != -1) ? xx[indices[0]] : 1.0;
      else if (dt < center_x_split)
        return (indices[1] != -1) ? xx[indices[1]] : 1.0;
      else if (dt < right_x_split)
        return (indices[2] != -1) ? xx[indices[2]] : 1.0;
      else
        return (indices[3] != -1) ? xx[indices[3]] : 1.0;
    }

    return 1.0;
  }
}
