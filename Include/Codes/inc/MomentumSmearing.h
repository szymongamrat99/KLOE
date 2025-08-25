#ifndef MOMENTUMSMEARING_H
#define MOMENTUMSMEARING_H

#include <TRandom3.h>    // For generation of gaussian variables
#include <TMatrix.h>     // For Matrix operations
#include <TDecompChol.h> // For Cholesky Decomposition V = LLT
#include <TVector.h>     // For Vector operations
#include <TBufferJSON.h> // To store results in JSON file
#include <vector>        // For storing multiple data vectors
#include <string>
#include <iostream>
#include <fstream>
#include "json.hpp"      // Assuming a JSON library is available, e.g., nlohmann/json

// Standard KLOE class - assuming this is part of your project
#include <kloe_class.h>  

namespace KLOE
{
  enum class CovMatrixMode
  {
    NEW_COV_MATRIX,
    MOM_SMEAR
  };

  template <typename T = Double_t>
  class MomentumSmearing : private pm00
  {
  private:
    TVectorT<T>
        _momVecMC; // Current momentum vector to be smeared

    TMatrixT<T>
        _covMatrix, // The final covariance matrix
        _U,
        _UT;

    Int_t
        _vecSize;

    CovMatrixMode
        _flag;

    // We store the full differences to allow for proper bootstrapping
    std::vector<TVectorT<T>>
        _diffVectors; 

    // Helper functions for internal use
    void _CovMatrixToJSON();
    void _CholeskyDecomposition();
    TMatrixT<T> _calculateCovarianceFromSample(const std::vector<Int_t>& indices);

  public:
    // Constructor for covariance matrix calculation mode
    MomentumSmearing(TVectorT<T> &momVecMC, TVectorT<T> &momVecData, TMatrixT<T> &covMatrix); 
    // Constructor for momentum smearing mode
    MomentumSmearing(TVectorT<T> &momVecMC, TMatrixT<T> &covMatrix);                          

    // Accumulates a new data point for covariance calculation
    void AddCovariancePoint(TVectorT<T> &momVecMC, TVectorT<T> &momVecData); 
    // Calculates the final covariance matrix from accumulated points
    void GetCovMatrix();
    // Calculates the uncertainty of the covariance matrix using bootstrap
    void CovMatrixUncertainty(Int_t numberOfSamples);

    // Applies momentum smearing to the MC vector
    void SmearMomentum(); 
    void GetSmearedMomentum(TVectorT<T> &momVecMC);

    void SetCovMatrix(TMatrixT<T> &covMatrix)
    {
      _covMatrix = covMatrix;
    }

    void SetMCVector(TVectorT<T> &momVecMC)
    {
      _momVecMC = momVecMC;
    }

  };
};

#endif // !MOMENTUMSMEARING_H