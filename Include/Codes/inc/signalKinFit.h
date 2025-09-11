#pragma once

#include <TVectorD.h>
#include <TMatrixD.h>
#include <TMath.h>
#include <ConfigManager.h>
#include <charged_mom.h>
#include <neutral_mom.h>

#include <KinFitter.h>
#include <kloe_class.h>

namespace KLOE
{
	class SignalKinFit : public KinFitter, protected ChargedVtxRec<>
	{
	private:
		TVectorD
				_X,
				_C,
				_X_final,
				_L,
				_CORR,
				_X_init,
				_X_min,
				_C_min,
				_L_min,
				_C_aux,
				_L_aux,
				_X_init_min,
				_X_init_aux;

		TMatrixD
				_V,
				_D,
				_D_T,
				_V_final,
				_V_aux,
				_V_min,
				_Aux,
				_V_invert,
				_V_init;

		Double_t
				_min_value_def,
				_FUNVALMIN,
				_CHISQRMIN,
				_Chi2TriKinFit;

		Bool_t
				_isConverged;

		std::unique_ptr<Int_t[]>
				_ind_gam;

		std::vector<Float_t>
				_cluster[5],
				_bhabha_mom,
				_bhabha_mom_err,
				_bhabha_vtx,
				_Param,
				_Errors,
				_iptri_kinfit,
				_fourKnetri_kinfit,
				_neu_vtx_min,
				_gamma_mom_final[4];

		std::vector<Int_t>
				_NeuClusters,
				_g4takentri_kinfit,
				_selected;

		ConfigManager
				&_config = ConfigManager::getInstance();

		TrilaterationCode
				_recMode;

		Reconstructor
				_R;
		Solution
				_S;

	public:
		SignalKinFit(Int_t N_free, Int_t N_const, Int_t M, Int_t loopcount, Double_t chiSqrStep, Int_t jmin, Int_t jmax, ErrorHandling::ErrorLogs &logger);
		~SignalKinFit();

		void SetParameters(const std::vector<Float_t> clu[5], const std::vector<Int_t> NeuClusters, const std::vector<Float_t> bhabha_mom, const std::vector<Float_t> bhabha_mom_err, const std::vector<Float_t> bhabha_vtx)
		{
			for (Int_t i = 0; i < 5; i++)
				_cluster[i].assign(clu[i].begin(), clu[i].end());

			_NeuClusters = NeuClusters;
			_bhabha_mom = bhabha_mom;
			_bhabha_mom_err = bhabha_mom_err;
			_bhabha_vtx = bhabha_vtx;
		}

		void GetResults(Int_t &bunchnum, std::vector<Float_t> &iptri_kinfit, std::vector<Int_t> &g4takentri_kinfit, std::vector<Float_t> gamma_mom_final[4], std::vector<Float_t> &fourKnetri_kinfit, std::vector<Float_t> &neu_vtx_min, Float_t &Chi2TriKinFit)
		{
			iptri_kinfit = _iptri_kinfit;
			g4takentri_kinfit = _g4takentri_kinfit;
			for (Int_t i = 0; i < 4; i++)
				gamma_mom_final[i] = _gamma_mom_final[i];
			fourKnetri_kinfit = _fourKnetri_kinfit;
			neu_vtx_min = _neu_vtx_min;
			Chi2TriKinFit = _Chi2TriKinFit;
		}

		ErrorHandling::ErrorCodes Reconstruct();
	};

} // namespace KLOE
