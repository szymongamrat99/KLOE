#include <ErrorLogs.h>
#include <KinFitter.h>
#include <uncertainties.h>
#include <reconstructor.h>

#include <signalKinFit.h>

namespace KLOE
{
	SignalKinFit::SignalKinFit(Int_t N_free, Int_t N_const, Int_t M, Int_t loopcount, Double_t chiSqrStep, ErrorHandling::ErrorLogs &logger) : KinFitter("SignalGlobal", N_free, N_const, M, 0, loopcount, chiSqrStep, logger), _V(N_free + N_const, N_free + N_const), _D(M, N_free + N_const), _D_T(N_free + N_const, M), _V_final(N_free + N_const, N_free + N_const), _V_aux(N_free + N_const, N_free + N_const), _V_min(N_free + N_const, N_free + N_const), _Aux(M, M), _V_invert(N_free, N_free), _V_init(N_free + N_const, N_free + N_const), _X(N_free + N_const), _C(M), _X_final(N_free + N_const), _L(M), _CORR(N_free + N_const), _X_init(N_free + N_const), _X_min(N_free + N_const), _C_min(M), _L_min(M), _C_aux(M), _L_aux(M), _X_init_min(N_free + N_const), _X_init_aux(N_free + N_const), _Param(N_free + N_const), _Errors(N_free + N_const)
	{
		for (Int_t i = 0; i < 4; i++)
			_photonFit[i].resize(8);

		_ipFit.resize(3);
		_KchrecFit.resize(10);
		_KchboostFit.resize(10);
		_KnerecFit.resize(10);
		_KnereclorFit.resize(10);
		for (Int_t i = 0; i < 2; i++)
			_trkFit[i].resize(4);

		KinFitter::ConstraintSet({"PxConsvLAB",
								  "PyConsvLAB",
								  "PzConsvLAB",
								  "EnergyConsvLAB",
								  "MinvConsvNeutralKaon",
                  "MinvConsvChargedKaon",
								  "Photon1PathLAB",
								  "Photon2PathLAB",
								  "Photon3PathLAB",
								  "Photon4PathLAB"});

		gErrorIgnoreLevel = 6001;
	}

	SignalKinFit::~SignalKinFit()
	{
	}

	ErrorHandling::ErrorCodes SignalKinFit::Reconstruct()
	{
		_CHISQRMIN = 999999.;
		_isConverged = 0;

		Bool_t
			clusterEnergy = false,
			cond_clus[4] = {false, false, false, false};

		clusterEnergy = _cluster[0][4] > MIN_CLU_ENE &&
						        _cluster[1][4] > MIN_CLU_ENE &&
						        _cluster[2][4] > MIN_CLU_ENE &&
						        _cluster[3][4] > MIN_CLU_ENE;

		for (Int_t k = 0; k < 4; k++)
		{
			cond_clus[k] =
				_cluster[k][0] != 0 &&
				_cluster[k][1] != 0 &&
				_cluster[k][2] != 0;
		}

		Bool_t cond_tot = 1;

		if (cond_tot)
		{
			try
			{
				_offset = 0;

				for (Int_t i = 0; i < 4; i++)
				{
					for (Int_t j = 0; j < 5; j++)
						_Param[_offset + i * 5 + j] = _cluster[i][j];

					_Errors[_offset + i * 5] = clu_x_error(_Param[_offset + i * 5], _Param[_offset + i * 5 + 1], _Param[_offset + i * 5 + 2], _Param[_offset + i * 5 + 4]);		// cm
					_Errors[_offset + i * 5 + 1] = clu_y_error(_Param[_offset + i * 5], _Param[_offset + i * 5 + 1], _Param[_offset + i * 5 + 2], _Param[_offset + i * 5 + 4]); // cm
					_Errors[_offset + i * 5 + 2] = clu_z_error(_Param[_offset + i * 5], _Param[_offset + i * 5 + 1], _Param[_offset + i * 5 + 2], _Param[_offset + i * 5 + 4]); // cm
					_Errors[_offset + i * 5 + 3] = clu_time_error(_Param[_offset + i * 5 + 4]);																					// ns
					_Errors[_offset + i * 5 + 4] = clu_ene_error(_Param[_offset + i * 5 + 4]);																					// MeV
				}

				_offset = 20;

				for (Int_t i = 0; i < 3; i++)
				{
					_Param[_offset + i] = _chargedVtx[i];
					_Errors[_offset + i] = _chargedVtxErr[i];
				}

				_offset = 23;

				for (Int_t i = 0; i < 2; i++)
				{
					_Param[_offset + i * 3] = _trackParameters[i][0];
					_Param[_offset + i * 3 + 1] = _trackParameters[i][1];
					_Param[_offset + i * 3 + 2] = _trackParameters[i][2];

					_Errors[_offset + i * 3] = _trackParametersErr[i][0];
					_Errors[_offset + i * 3 + 1] = _trackParametersErr[i][1];
					_Errors[_offset + i * 3 + 2] = _trackParametersErr[i][2];
				}

				_offset = 29;

				for (Int_t i = 0; i < 4; i++)
				{
					_Param[_offset + i] = _bhabha_mom[i];
					_Errors[_offset + i] = _bhabha_mom_err[i];
				}

				_offset = 33;


				for (Int_t i = 0; i < 3; i++)
				{
					_Param[_offset + i] = _neuVtx[i];
					_Errors[_offset + i] = _neuVtxErr[i];
				}

				_offset = 36;

				for (Int_t i = 0; i < 3; i++)
				{
					_Param[_offset + i] = _bhabha_vtx[i];
					_Errors[_offset + i] = _bhabhaVtxErr[i];
				}

				_offset = _Param.size();

				if (_offset != _N_free + _N_const)
					throw ErrorHandling::ErrorCodes::SIGNAL_KIN_FIT;
			}
			catch (ErrorHandling::ErrorCodes &err)
			{
				return err;
			}

			KinFitter::ParameterInitialization(_Param.data(), _Errors.data());

			_CHISQRMIN = KinFitter::FitFunction();


			KinFitter::GetResults(_X_min, _V_min, _X_init_min, _V_init, _trkFit.data(), _KchrecFit, _KchboostFit, _ipFit, _photonFit.data(), _KnerecFit, _KnereclorFit);

			_isConverged = 1;
		}
		else
		{
			_isConverged = 0;
		}

		if (_isConverged)
			return ErrorHandling::ErrorCodes::NO_ERROR;
		else
			return ErrorHandling::ErrorCodes::SIGNAL_KIN_FIT;
	}
}