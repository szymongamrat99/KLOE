#ifndef MOMENTUMSMEARING_H
#define MOMENTUMSMEARING_H

#include <TRandom3.h>    // For generation of gaussian variables
#include <TMatrix.h>     // For Matrix operations
#include <TDecompChol.h> // For Cholesky Decomposition V = LLT
#include <TVector.h>     // For Vector operations
#include <TBufferJSON.h> // To store results in JSON file

#include <kloe_class.h>  // Standard KLOE class

namespace KLOE
{
  enum class CovMatrixMode
  {
    NEW_COV_MATRIX,
    MOM_SMEAR
  };

  template <typename T = Double_t>
  class MomentumSmearing : public pm00
  {
  private:
    TVectorT<T>
        _momVecMC,
        _momVecData;

    TMatrixT<T>
        _covMatrix,
        _U,
        _UT;

    Int_t
        _vecSize,
        _numberOfPoints = 0;

    CovMatrixMode
              _flag;

    void _CovMatrixToJSON();
    void _CholeskyDecomposition();

  public:
    MomentumSmearing(TVectorT<T> &momVecMC, TVectorT<T> &momVecData, TMatrixT<T> &covMatrix); // 1
    MomentumSmearing(TVectorT<T> &momVecMC, TMatrixT<T> &covMatrix);                          // 2

    void CovCalcPart(); // Usable only for constructor 1
    void GetCovMatrix();
    void CovMatrixUncertainty(Int_t numberOfSamples);


    void SmearMomentum(); // Usable only for constructor 2
    void GetSmearedMomentum(TVectorT<T> &momVecMC);

    void SetCovMatrix(TMatrixT<T> &covMatrix)
    {
      _covMatrix = covMatrix;
    }

    void SetMCVector(TVectorT<T> &momVecMC)
    {
      _momVecMC = momVecMC;
    }

    void SetDataVector(TVectorT<T> &momVecData)
    {
      _momVecData = momVecData;
    }

  };
};

#endif // !MOMENTUMSMEARING_H