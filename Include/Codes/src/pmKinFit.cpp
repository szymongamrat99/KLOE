#include <ErrorLogs.h>
#include <KinFitter.h>
#include <uncertainties.h>
#include <reconstructor.h>

#include <pmKinFit.h>

namespace KLOE
{
	PMKinFit::PMKinFit(Int_t N_free, Int_t N_const, Int_t M, Int_t loopcount, Double_t chiSqrStep, ErrorHandling::ErrorLogs &logger) : KinFitter("PM", N_free, N_const, M, 0, loopcount, chiSqrStep, logger), _V(N_free + N_const, N_free + N_const), _D(M, N_free + N_const), _D_T(N_free + N_const, M), _V_final(N_free + N_const, N_free + N_const), _V_aux(N_free + N_const, N_free + N_const), _V_min(N_free + N_const, N_free + N_const), _Aux(M, M), _V_invert(N_free, N_free), _V_init(N_free + N_const, N_free + N_const), _X(N_free + N_const), _C(M), _X_final(N_free + N_const), _L(M), _CORR(N_free + N_const), _X_init(N_free + N_const), _X_min(N_free + N_const), _C_min(M), _L_min(M), _C_aux(M), _L_aux(M), _X_init_min(N_free + N_const), _X_init_aux(N_free + N_const), _Param(N_free + N_const), _Errors(N_free + N_const)
	{

		_ipFit.resize(3);
		_KchrecFit.resize(10);
		_KchboostFit.resize(10);

		for (Int_t i = 0; i < 2; i++)
			_trkFit[i].resize(4);

		KinFitter::ConstraintSet(std::vector<std::string>
                             {"energyconsvcm",
								              "minvconsvchargedkaon"});

		gErrorIgnoreLevel = 6001;
	}

	PMKinFit::~PMKinFit()
	{
	}

	ErrorHandling::ErrorCodes PMKinFit::Reconstruct()
	{
		_CHISQRMIN = 999999.;
		_isConverged = 0;

		Bool_t cond_tot = 1;

		if (cond_tot)
		{
			try
			{
				_offset = 0;

				for (Int_t i = 0; i < 3; i++)
				{
					_Param[_offset + i] = _chargedVtx[i];
					_Errors[_offset + i] = _chargedVtxErr[i];
				}

				_offset = 3;

				for (Int_t i = 0; i < 2; i++)
				{
					_Param[_offset + i * 3] = _trackParameters[i][0];
					_Param[_offset + i * 3 + 1] = _trackParameters[i][1];
					_Param[_offset + i * 3 + 2] = _trackParameters[i][2];

					_Errors[_offset + i * 3] = _trackParametersErr[i][0];
					_Errors[_offset + i * 3 + 1] = _trackParametersErr[i][1];
					_Errors[_offset + i * 3 + 2] = _trackParametersErr[i][2];
				}

				_offset += 6;

				for (Int_t i = 0; i < 4; i++)
				{
					_Param[_offset + i] = _bhabha_mom[i];
					_Errors[_offset + i] = _bhabha_mom_err[i];
				}

				_offset += 4;

				for (Int_t i = 0; i < 3; i++)
				{
					_Param[_offset + i] = _bhabha_vtx[i];
					_Errors[_offset + i] = _bhabhaVtxErr[i];
				}

				_offset += 3;

				if (_offset != _N_free + _N_const)
					throw ErrorHandling::ErrorCodes::PM_KIN_FIT;
			}
			catch (ErrorHandling::ErrorCodes &err)
			{
				return err;
			}

			KinFitter::ParameterInitialization(_Param.data(), _Errors.data());

			_CHISQRMIN = KinFitter::FitFunction();

			KinFitter::GetResults(_X_min, _V_min, _X_init_min, _V_init, _trkFit.data(), _KchrecFit, _KchboostFit, _ipFit, _phiMomFit);

			_isConverged = 1;
		}
		else
		{
			_isConverged = 0;
		}

		if (_isConverged)
			return ErrorHandling::ErrorCodes::NO_ERROR;
		else
			return ErrorHandling::ErrorCodes::PM_KIN_FIT;
	}
}