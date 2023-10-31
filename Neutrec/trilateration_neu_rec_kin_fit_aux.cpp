#include <string.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "Math/Functor.h"
#include "Math/Factory.h"
#include "Math/Minimizer.h"

#include "../../Include/Codes/reconstructor.h"
#include "../../Include/const.h"
#include "../../Include/Codes/uncertainties.h"
#include "../../Include/Codes/charged_mom.h"
#include "../../Include/Codes/neutral_mom.h"
#include "../../Include/Codes/lorentz_transf.h"
#include "../../Include/Codes/plane_intersection.h"
#include "../../Include/Codes/closest_approach.h"
#include "../../Include/Codes/kloe_class.h"
#include "../../Include/Codes/kinematic_fits.h"
#include "chain_init.C"

const UInt_t N_free = 27, N_const = 10, M = 9;
const Float_t Trf = 2.715; // ns - time of a bunch (correction)

Float_t P[N_free + N_const], P1[N_free + N_const], S[N_free]; 

Reconstructor R;
Solution S;

Int_t selected[4] = {1, 2, 3, 4};

using namespace KLOE;

Double_t trilateration_chi_square(const Double_t *x)
{
  kin_fits event(N_free, )

  for (Int_t i = 0; i < 4; i++)
		R.SetClu(i, clusters[0][i],
						 clusters[1][i],
						 clusters[2][i],
						 clusters[3][i],
						 clusters[4][i]);

	R.SetClu(4, 0., 0., 0., 0., 0.);
	R.SetClu(5, 0., 0., 0., 0., 0.);

	S = R.MySolve(selected);

	for (Int_t i = 0; i < 2; i++)
		for (Int_t j = 0; j < 4; j++)
			neu_vtx[i][j] = S.sol[i][j];
	//!
}