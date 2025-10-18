#include <algorithm>
#include <random>

#include <TH1.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>
#include <TRatioPlot.h>

#include "interference.h"

namespace KLOE
{
	Double_t interference::interf_function(const Float_t x, Int_t check, const Double_t *par)
	{
		Double_t Value = 0;
		Double_t Epsilon = 0, RePart = 0, Dphi = 0, TauKs = 0, TauKl = 0, MassDiff = 0;
		Double_t ImPart = 0, GammaKl = 0, GammaKs = 0, Gamma = 0, DMass = 0;

		Float_t dt = x;

		// Parameters from PDG2023

		Epsilon = PhysicsConstants::mod_epsilon;
		Dphi = PhysicsConstants::phi_pm_nonCPT - PhysicsConstants::phi_00_nonCPT; // phi(+-)-phi(00) (degrees)
		TauKs = PhysicsConstants::tau_S_nonCPT * pow(10, -9);		// PDG fit not assuming CPT (s)
		TauKl = PhysicsConstants::tau_L * pow(10, -9);					// Kl mean life (s)
		MassDiff = PhysicsConstants::delta_mass_nonCPT;					// M(Kl)-M(Ks) ( (h/2pi)s-1 ):
																					// PDG fit not assuming CPT

		if (check == 0)
		{
			RePart = par[0];
			ImPart = par[1];
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
			Value = (1. + 2. * RePart) * exp(-GammaKs * abs(dt)) +
							(1. - 4. * RePart) * exp(-GammaKl * abs(dt)) -
							2. * exp(-0.5 * Gamma * abs(dt)) *
									((1. - RePart) * cos(DMass * abs(dt)) -
									 3. * ImPart * sin(DMass * abs(dt)));
		}

		return (pow(Epsilon, 2) / (2. * Gamma)) * Value * 100000;
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////

	void interference::bin_extraction(UInt_t channel, TH1 *histogram)
	{
		if (channel < KLOE::channNum)
		{
			for (Int_t i = 0; i < bin_number; i++)
				b[channel][i] = (histogram->GetBinContent(i + 1));

			for (Int_t i = 0; i < bin_number; i++)
				e[channel][i] = (histogram->GetBinError(i + 1));
		}
		else
		{
			for (Int_t i = 0; i < bin_number; i++)
				b_data[i] = (histogram->GetBinContent(i + 1));

			for (Int_t i = 0; i < bin_number; i++)
				e_data[i] = (histogram->GetBinError(i + 1));
		}
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

		Double_t ReFit = xx[0];
		Double_t ImFit = xx[1];
		Double_t Norm[9];

		Norm[0] = xx[2];	// Signal norm
		Norm[1] = xx[3];	// Left DC Wall
		Norm[2] = xx[4];	// Left beam pipe Wall
		Norm[3] = xx[5];	// Right beam pipe Wall
		Norm[4] = xx[6];	// Right DC Wall
		Norm[5] = xx[7];	// Omega norm
		Norm[6] = xx[8];	// Three norm
		Norm[7] = xx[9];	// Semi norm
		Norm[8] = xx[10]; // Other norm

		/////////////////////////////////////////////////////////////////////////////////////////////
		for (Int_t i = 0; i < KLOE::channNum; i++)
		{
			for (Int_t j = 0; j < time_diff[i].size(); j++)
			{
				if (i == 0)
					frac[i]->Fill(time_diff[i][j], interf_function(time_diff_gen[j], 0, xx)); //! Filling Signal
				else
					frac[i]->Fill(time_diff[i][j]); //! Filling background
			}

			//! Using correction factor and efficiency
			if (i == 0)
			{
				frac[i]->Scale(frac[i]->GetEntries() / frac[i]->Integral(0, bin_number + 1));

				for (UInt_t j = 0; j < bin_number; j++)
				{
					frac[i]->SetBinContent(i + 1, frac[i]->GetBinContent(i + 1) * corr_vals[i]);
				}

				interference::bin_extraction(i, frac[i]);
			}
			else
				interference::bin_extraction(i, frac[i]);
		}

		/////////////////////////////////////////////////////////////////////////////////////////////

		for (Int_t j = 0; j < time_diff_data.size(); j++)
			data->Fill(time_diff_data[j]); //! Filling DATA

		interference::bin_extraction(6, data);

		/////////////////////////////////////////////////////////////////////////////////////////////

		for (Int_t i = 0; i < KLOE::channNum; i++)
			frac[i]->Reset("ICESM");

		data->Reset("ICESM");

		/////////////////////////////////////////////////////////////////////////////////////////////

		Double_t value = 0;

		//! Sum of bins and errors for the fitting further

		for (Int_t i = 0; i < bin_number; i++)
		{
			b_mcsum[i] = 0.;
			e_mcsum[i] = 0.;
		}

		for (Int_t i = 0; i < KLOE::channNum; i++)
		{
			for (Int_t j = 0; j < bin_number; j++)
			{
				if (i == 1)
				{
					if (j + 1 < frac[i]->FindBin(left_x_split))
					{
						b_mcsum[j] += Norm[1] * b[i][j];
						e_mcsum[j] += pow(Norm[1] * e[i][j], 2);
					}
					else if (j + 1 > frac[i]->FindBin(left_x_split) && j + 1 < frac[i]->FindBin(center_x_split))
					{
						b_mcsum[j] += Norm[2] * b[i][j];
						e_mcsum[j] += pow(Norm[2] * e[i][j], 2);
					}
					else if (j + 1 > frac[i]->FindBin(center_x_split) && j + 1 < frac[i]->FindBin(right_x_split))
					{
						b_mcsum[j] += Norm[3] * b[i][j];
						e_mcsum[j] += pow(Norm[3] * e[i][j], 2);
					}
					else if (j + 1 > frac[i]->FindBin(right_x_split))
					{
						b_mcsum[j] += Norm[4] * b[i][j];
						e_mcsum[j] += pow(Norm[4] * e[i][j], 2);
					}
				}
				else if (i == 0)
				{
					b_mcsum[j] += Norm[i] * b[i][j];
					e_mcsum[j] += pow(Norm[i] * e[i][j], 2);
				}
				else
				{
					b_mcsum[j] += Norm[i + 3] * b[i][j];
					e_mcsum[j] += pow(Norm[i + 3] * e[i][j], 2);
				}
			}
		}

		for (Int_t i = 0; i < bin_number; i++)
		{
			value += pow(b_data[i] - b_mcsum[i], 2) / (pow(e_data[i], 2) + (e_mcsum[i]));
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
						frac[i]->Fill(time_diff[i][j], interf_function(time_diff_gen[j], 0, xx));
					else
						frac[i]->Fill(time_diff[i][j]);
				}
			}

			if (i == 0)
			{
				frac[i]->Scale(frac[i]->GetEntries() / frac[i]->Integral(0, bin_number + 1));
				interference::bin_extraction(i, frac[i]);
			}
			else
				interference::bin_extraction(i, frac[i]);
		}

		for (Int_t i = 0; i < KLOE::channNum - 2; i++)
			frac[i]->Reset("ICESM");
		/////////////////////////////////////////////////////////////////////////////////////////////

		Double_t value = 0, *bin_sum, *err_sum;

		bin_sum = new Double_t[bin_number];
		err_sum = new Double_t[bin_number];

		for (Int_t i = 0; i < bin_number; i++)
		{
			bin_sum[i] = 0.;
			err_sum[i] = 0.;
		}

		for (Int_t i = 0; i < KLOE::channNum - 3; i++)
			for (Int_t j = 0; j < bin_number; j++)
			{
				bin_sum[j] += Norm[i] * b[i][j];
				err_sum[j] += pow(Norm[i] * e[i][j], 2);
			}

		for (Int_t i = 0; i < bin_number; i++)
			value += pow(b[6][i] - bin_sum[i], 2) / (pow(e[6][i], 2) + err_sum[i] / 2.);

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

		for (Int_t j = 0; j < time_diff_rand_mc[0].size(); j++)
			frac[0]->Fill(time_diff_rand_mc[0][j], interf_function(time_diff_gen_rand_mc[j], 0, xx));

		for (Int_t j = 0; j < time_diff_rand_data[0].size(); j++)
			data->Fill(time_diff_rand_data[0][j], interf_function(time_diff_gen_rand_data[j]));

		frac[0]->Scale(frac[0]->GetEntries() / frac[0]->Integral(0, bin_number + 1));
		interference::bin_extraction(0, frac[0]);

		data->Scale(data->GetEntries() / data->Integral(0, bin_number + 1));
		interference::bin_extraction(6, data);

		frac[0]->Reset("ICESM");
		data->Reset("ICESM");

		/////////////////////////////////////////////////////////////////////////////////////////////

		Double_t value = 0;

		for (Int_t i = 0; i < bin_number; i++)
		{
			b_mcsum[i] = 0.;
			e_mcsum[i] = 0.;
		}

		for (Int_t j = 0; j < bin_number; j++)
		{
			b_mcsum[j] += Norm * b[0][j];
			e_mcsum[j] += pow(Norm * e[0][j], 2);
		}

		for (Int_t i = 0; i < bin_number; i++)
			value += pow(b_data[i] - b_mcsum[i], 2) / (pow(e_data[i], 2) + e_mcsum[i]);

		return value;
	};

	//! Fitting signal MC to DATA after efficiency correction
	Double_t interference::interf_chi2_mc_data(const Double_t *xx)
	{
		Double_t ReFit = xx[0];
		Double_t ImFit = xx[1];
		Double_t Norm_signal = xx[2];
		Double_t Norm[8];

		Norm[0] = tmp_norm[0]; // Left DC Wall
		Norm[1] = tmp_norm[1]; // Left beam pipe Wall
		Norm[2] = tmp_norm[2]; // Right beam pipe Wall
		Norm[3] = tmp_norm[3]; // Right DC Wall
		Norm[4] = tmp_norm[4]; // Omega norm
		Norm[5] = tmp_norm[5]; // Three norm
		Norm[6] = tmp_norm[6]; // Semi norm
		Norm[7] = tmp_norm[7]; // Other bcg norm

		/////////////////////////////////////////////////////////////////////////////////////////////
		for (Int_t i = 0; i < KLOE::channNum; i++)
		{
			for (Int_t j = 0; j < time_diff[i].size(); j++)
			{
				if (i == 0)
					frac[i]->Fill(time_diff[i][j], interf_function(time_diff[i][j], 0, xx)); //! Filling Signal
				else if (i == 1)
				{
					if (time_diff[i][j] < left_x_split)
						frac[i]->Fill(time_diff[i][j], Norm[0]);
					else if (time_diff[i][j] > left_x_split && time_diff[i][j] < center_x_split)
						frac[i]->Fill(time_diff[i][j], Norm[1]);
					else if (time_diff[i][j] > center_x_split && time_diff[i][j] < right_x_split)
						frac[i]->Fill(time_diff[i][j], Norm[2]);
					else if (time_diff[i][j] > right_x_split)
						frac[i]->Fill(time_diff[i][j], Norm[3]);
				}
				else
				{
					frac[i]->Fill(time_diff[i][j], Norm[i + 2]); //! Filling background
				}
			}

			//! Using correction factor and efficiency
			if (i == 0)
			{
				frac[i]->Scale(frac[i]->GetEntries() / frac[i]->Integral(0, bin_number + 1));

				for (Int_t k = 1; k <= bin_number; k++)
					frac[i]->SetBinContent(k, frac[i]->GetBinContent(k) / corr_vals[k - 1]);

				interference::bin_extraction(i, frac[i]);
			}
		}

		/////////////////////////////////////////////////////////////////////////////////////////////

		for (Int_t j = 0; j < time_diff_data.size(); j++)
			data->Fill(time_diff_data[j]); //! Filling DATA

		for (Int_t i = 1; i < KLOE::channNum; i++)
		{
			data->Add(frac[i], -1.);
		}

		for (Int_t k = 1; k <= bin_number; k++)
			data->SetBinContent(k, data->GetBinContent(k) / corr_vals[k - 1]);

		interference::bin_extraction(6, data);

		/////////////////////////////////////////////////////////////////////////////////////////////

		for (Int_t i = 0; i < KLOE::channNum; i++)
			frac[i]->Reset("ICESM");

		data->Reset("ICESM");

		/////////////////////////////////////////////////////////////////////////////////////////////

		Double_t value = 0;

		//! Sum of bins and errors for the fitting further

		for (Int_t i = 0; i < bin_number; i++)
		{
			b_mcsum[i] = 0.;
			e_mcsum[i] = 0.;
		}

		for (Int_t j = 0; j < bin_number; j++)
		{

			b_mcsum[j] += Norm_signal * b[0][j];
			e_mcsum[j] += pow(Norm_signal * e[0][j], 2);
		}

		for (Int_t i = 0; i < bin_number; i++)
		{
			value += pow(b_data[i] - b_mcsum[i], 2) / (pow(e_data[i], 2) + e_mcsum[i]);
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

		for (Int_t i = 0; i < KLOE::channNum; i++)
		{
			for (Int_t j = 0; j < time_diff_rand_mc[i].size(); j++)
			{
				if (i == 0)
					frac[i]->Fill(time_diff_rand_mc[i][j], interf_function(time_diff_gen_rand_mc[j], 0, xx));
				else
					frac[i]->Fill(time_diff_rand_mc[i][j]);
			}

			for (Int_t j = 0; j < time_diff_rand_data[i].size(); j++)
			{
				if (i == 0)
					frac_data[i]->Fill(time_diff_rand_data[i][j], interf_function(time_diff_gen_rand_data[j]));
				else
					frac_data[i]->Fill(time_diff_rand_data[i][j]);
			}

			frac[i]->Scale(frac[i]->GetEntries() / frac[i]->Integral(0, bin_number + 1));
			interference::bin_extraction(i, frac[i]);

			frac_data[i]->Scale(frac_data[i]->GetEntries() / frac_data[i]->Integral(0, bin_number + 1));
			data->Add(frac_data[i]);
		}

		interference::bin_extraction(6, data);

		for (Int_t i = 0; i < 6; i++)
		{
			frac[i]->Reset("ICESM");
			frac_data[i]->Reset("ICESM");
		}

		data->Reset("ICESM");

		/////////////////////////////////////////////////////////////////////////////////////////////

		Double_t value = 0;

		for (Int_t i = 0; i < bin_number; i++)
		{
			b_mcsum[i] = 0.;
			e_mcsum[i] = 0.;
		}

		for (Int_t i = 0; i < KLOE::channNum; i++)
			for (Int_t j = 0; j < bin_number; j++)
			{
				b_mcsum[j] += Norm[i] * b[i][j];
				e_mcsum[j] += pow(Norm[i] * e[i][j], 2);
			}

		for (Int_t i = 0; i < bin_number; i++)
			value += pow(b_data[i] - b_mcsum[i], 2) / (pow(e_data[i], 2) + e_mcsum[i]);

		return value;
	};

	//! Auxillary function to use any of the others
	Double_t interference::interf_chi2(const Double_t *xx)
	{
		if (mode == "split")
			return interf_chi2_split(xx);
		else if (mode == "excluded")
			return interf_chi2_excluded(xx);
		else if (mode == "mc")
			return interf_chi2_mc(xx);
		else if (mode == "bcg")
			return interf_chi2_bcg(xx);
		else if (mode == "final")
			return interf_chi2_mc_data(xx);
		else
			return -999.0;
	};

	void interference::ClearVectors()
	{
		time_diff[0].clear();
		time_diff[1].clear();
		time_diff[2].clear();
		time_diff[3].clear();
		time_diff[4].clear();
		time_diff[5].clear();
		time_diff_data.clear();
		time_diff_gen.clear();

		time_diff[0].shrink_to_fit();
		time_diff[1].shrink_to_fit();
		time_diff[2].shrink_to_fit();
		time_diff[3].shrink_to_fit();
		time_diff[4].shrink_to_fit();
		time_diff[5].shrink_to_fit();
		time_diff_data.shrink_to_fit();
		time_diff_gen.shrink_to_fit();
	}
}