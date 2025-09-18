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
			_CHISQRMIN,
			_Chi2TriKinFit;

		Bool_t
			_isConverged;

		std::vector<Float_t>
			_trackParameters[2],
			_trackParametersErr[2],
			_cluster[4],
			_chargedVtx,
			_chargedVtxErr,
			_bhabha_mom,
			_bhabha_mom_err,
			_bhabha_vtx,
			_Param,
			_Errors,
			_KchrecFit,
			_KchboostFit,
			_ipFit,
			_KnereclorFit,
			_KnerecFit,
			_photonFit[4],
			_trkFit[2];

		Int_t
			_offset;

		ConfigManager
			&_config = ConfigManager::getInstance();

	public:
		SignalKinFit(Int_t N_free, Int_t N_const, Int_t M, Int_t loopcount, Double_t chiSqrStep, ErrorHandling::ErrorLogs &logger);
		~SignalKinFit();

		void SetParameters(const std::vector<Float_t> trackParameters[2], const std::vector<Float_t> trackParametersErr[2], const std::vector<Float_t> cluster[4], const std::vector<Float_t> chargedVtx, const std::vector<Float_t> chargedVtxErr, const std::vector<Float_t> bhabha_mom, const std::vector<Float_t> bhabha_mom_err, const std::vector<Float_t> bhabha_vtx)
		{
			for (Int_t i = 0; i < 2; i++)
			{
				_trackParameters[i] = trackParameters[i];
				_trackParametersErr[i] = trackParametersErr[i];
			}

			for (Int_t i = 0; i < 4; i++)
				_cluster[i].assign(cluster[i].begin(), cluster[i].end());

			_chargedVtx = chargedVtx;
			_chargedVtxErr = chargedVtxErr;

			_bhabha_mom = bhabha_mom;
			_bhabha_mom_err = bhabha_mom_err;

			_bhabha_vtx = bhabha_vtx;
		}

		void GetResults(std::vector<Float_t> &Param, std::vector<Float_t> &Errors, std::vector<Float_t> &ParamFit, std::vector<Float_t> &ErrorsFit, std::vector<Float_t> trkFit[2], std::vector<Float_t> &KchrecFit, std::vector<Float_t> &KchboostFit, std::vector<Float_t> &ipFit, std::vector<Float_t> photonFit[4], std::vector<Float_t> &KnerecFit, std::vector<Float_t> &KnereclorFit, Float_t &Chi2SignalKinFit)
		{
			// Param.resize(_N_free);
			// Errors.resize(_N_free);
			// Param.assign(_Param.begin(), _Param.begin() + _N_free);
			// Errors.assign(_Errors.begin(), _Errors.begin() + _N_free);

			// ParamFit.resize(_N_free);
			// ErrorsFit.resize(_N_free);
			// ParamFit.assign(&_X_final.GetMatrixArray()[0], &_X_final.GetMatrixArray()[0] + _N_free);
			// ErrorsFit.assign(&_V_final.GetMatrixArray()[0], &_X_final.GetMatrixArray()[0] + _N_free);

			for (Int_t i = 0; i < 2; i++)
				trkFit[i] = _trkFit[i];

			KchrecFit = _KchrecFit;
			KchboostFit = _KchboostFit;

			ipFit = _ipFit;

			for (Int_t i = 0; i < 4; i++)
				photonFit[i] = _photonFit[i];

			KnereclorFit = _KnereclorFit;
			KnerecFit = _KnerecFit;

			Chi2SignalKinFit = _CHISQRMIN;
		}

		ErrorHandling::ErrorCodes Reconstruct();
	};

} // namespace KLOE
