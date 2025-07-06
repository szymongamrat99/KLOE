#include <MomentumSmearing.h>

#include <SystemPaths.h>

namespace KLOE
{
  template class MomentumSmearing<Double_t>;
  template class MomentumSmearing<Float_t>;

  template <typename T>
  MomentumSmearing<T>::MomentumSmearing(TVectorT<T> &momVecMC, TVectorT<T> &momVecData, TMatrixT<T> &covMatrix) : _momVecMC(momVecMC), _momVecData(momVecData), _covMatrix(covMatrix)
  {
    if (_momVecData.GetNoElements() == _momVecMC.GetNoElements() && _momVecMC.GetNoElements() == _covMatrix.GetNcols() && _momVecMC.GetNoElements() == _covMatrix.GetNrows())
    {
      _vecSize = _momVecData.GetNoElements();
      _covMatrix.Zero();
      _flag = CovMatrixMode::NEW_COV_MATRIX;
    }
    else
    {
      std::cout << "Vectors and Matrix have different sizes!" << std::endl;
      return;
    }
  };

  template <typename T>
  MomentumSmearing<T>::MomentumSmearing(TVectorT<T> &momVecMC, TMatrixT<T> &covMatrix) : _momVecMC(momVecMC), _covMatrix(covMatrix)
  {
    if (_momVecMC.GetNoElements() == _covMatrix.GetNcols() && _momVecMC.GetNoElements() == _covMatrix.GetNrows())
    {
      _vecSize = _momVecData.GetNoElements();
      _CholeskyDecomposition();
      _flag = CovMatrixMode::MOM_SMEAR;
    }
    else
    {
      std::cout << "Vectors and covariant matrix have different sizes!" << std::endl;
      return;
    }
  };

  template <typename T>
  void MomentumSmearing<T>::CovCalcPart()
  {
    for (Int_t i = 0; i < _vecSize; i++)
    {
      for (Int_t j = 0; j < _vecSize; j++)
      {
        _covMatrix(i, j) += (_momVecMC(i) - _momVecData(i)) * (_momVecMC(j) - _momVecData(j));
      }
    }
    _numberOfPoints++;
  }

  template <typename T>
  void MomentumSmearing<T>::GetCovMatrix()
  {
    switch (_flag)
    {
    case CovMatrixMode::NEW_COV_MATRIX:
    {
      if (_numberOfPoints != 0)
      {
        for(Int_t i = 0; i < _covMatrix.GetNrows(); i++)
          for(Int_t j = 0; j < _covMatrix.GetNcols(); j++)
          {
            _covMatrix(i,j) /= _numberOfPoints;
            // _covMatrix(i,j) = _covMatrix(i,j) / sqrt(_covMatrix(i,i) * _covMatrix(j,j));
          }
      }
      else
      {
        std::cout << "Number of points = 0!" << std::endl;
        return;
      }
      break;
    }
    }
    _CholeskyDecomposition();
    _CovMatrixToJSON();

    _covMatrix.Print();
  }

  template <typename T>
  void MomentumSmearing<T>::CovMatrixUncertainty(Int_t numberOfSamples)
  {
    TRandom3 
            randGen(0);
  }

  template <typename T>
  void MomentumSmearing<T>::SmearMomentum()
  {
    TRandom3 randGen(0);
    TVectorT<T>
        gaussVec(_vecSize),
        deltaMom(_vecSize);

    gaussVec.ResizeTo(_momVecMC);
    deltaMom.ResizeTo(_momVecMC);

    for (Int_t i = 0; i < _momVecMC.GetNoElements(); i++)
    {
      gaussVec(i) = randGen.Gaus(0.0, 1.0); // Generation from normal distribution
    }

    deltaMom = _UT * gaussVec; // Smearing vector calculation

    // std::cout << "Gaussian vector: " << std::endl;
    // gaussVec.Print();
    // std::cout << "Smearing vector: " << std::endl;
    // deltaMom.Print();
    // std::cout << "Original momentum vector: " << std::endl;
    // _momVecMC.Print();
    // std::cout << "Smeared momentum vector: " << std::endl;
    // _momVecMC.Print();

  

    _momVecMC += deltaMom;

  }

  template <typename T>
  void MomentumSmearing<T>::GetSmearedMomentum(TVectorT<T> &momVecMC)
  {
    momVecMC = _momVecMC;
  }

  template <typename T>
  void MomentumSmearing<T>::_CovMatrixToJSON()
  {
    std::string json = (std::string)TBufferJSON::ToJSON(&_covMatrix);
    properties["momSmearing"]["covarianceMatrix"] = json;

    std::ofstream outfile(SystemPath::generalPropertiesPath);
    outfile << properties.dump(4);
    outfile.close();
  }

  template <typename T>
  void MomentumSmearing<T>::_CholeskyDecomposition()
  {
    // Cholesky Decomposition to get the matrix for gaussian transformation
    TDecompChol decomposition(_covMatrix);
    decomposition.Decompose();

    _U.ResizeTo(decomposition.GetU());
    _UT.ResizeTo(decomposition.GetU());
    
    _U = decomposition.GetU();
    _UT.Transpose(_U);
  }

}