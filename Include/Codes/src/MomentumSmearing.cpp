#include <MomentumSmearing.h>
#include <const.h>

namespace KLOE
{
  template class MomentumSmearing<Double_t>;
  template class MomentumSmearing<Float_t>;

  template <typename T>
  MomentumSmearing<T>::MomentumSmearing(TVectorT<T> &momVecMC, TVectorT<T> &momVecData, TMatrixT<T> &covMatrix) : _momVecMC(momVecMC), _covMatrix(covMatrix)
  {
    if (momVecData.GetNoElements() == _momVecMC.GetNoElements())
    {
      _vecSize = _momVecMC.GetNoElements();
      _covMatrix.Zero();
      _diffVectors.clear(); // Clear the vector for new data
      _flag = CovMatrixMode::NEW_COV_MATRIX;
      
      // Add the first point
      AddCovariancePoint(momVecMC, momVecData);
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
      _vecSize = _momVecMC.GetNoElements();
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
  void MomentumSmearing<T>::AddCovariancePoint(TVectorT<T> &momVecMC, TVectorT<T> &momVecData)
  {
    // Check for size consistency, although constructors already do this
    if (momVecMC.GetNoElements() != _vecSize || momVecData.GetNoElements() != _vecSize) {
        std::cout << "ERROR: Input vectors do not match class dimension." << std::endl;
        return;
    }
    TVectorT<T> diffVec = momVecMC;
    diffVec -= momVecData;
    _diffVectors.push_back(diffVec);
  }

  template <typename T>
  void MomentumSmearing<T>::GetCovMatrix()
  {
    if (_diffVectors.empty())
    {
      std::cout << "Number of points = 0!" << std::endl;
      return;
    }
    
    // Calculate the mean difference vector
    TVectorT<T> meanDiff(_vecSize);
    meanDiff.Zero();
    for (const auto& diff : _diffVectors) {
        meanDiff += diff;
    }
    meanDiff *= (1.0 / _diffVectors.size());

    // Calculate the covariance matrix using loops (no Outer)
    _covMatrix.Zero();
    for (const auto& diff : _diffVectors) {
        TVectorT<T> centeredDiff = diff;
        centeredDiff -= meanDiff;
        for(Int_t i = 0; i < _vecSize; i++) {
            for(Int_t j = 0; j < _vecSize; j++) {
                _covMatrix(i, j) += centeredDiff(i) * centeredDiff(j);
            }
        }
    }
    _covMatrix *= (1.0 / (_diffVectors.size() - 1.0));

    _CholeskyDecomposition();
    std::cout << "Calculated Covariance Matrix:" << std::endl;
    _covMatrix.Print();
  }

  template <typename T>
  TMatrixT<T> MomentumSmearing<T>::_calculateCovarianceFromSample(const std::vector<Int_t>& indices)
  {
      TMatrixT<T> cov_matrix(_vecSize, _vecSize);
      cov_matrix.Zero();
      TVectorT<T> mean_diff(_vecSize);
      mean_diff.Zero();
      
      // Calculate mean of the resampled differences
      for (const auto& index : indices) {
          mean_diff += _diffVectors[index];
      }
      mean_diff *= (1.0 / indices.size());
      
      // Calculate covariance using loops (no Outer)
      for (const auto& index : indices) {
          TVectorT<T> centered_diff = _diffVectors[index];
          centered_diff -= mean_diff;
          for(Int_t i = 0; i < _vecSize; i++) {
              for(Int_t j = 0; j < _vecSize; j++) {
                  cov_matrix(i, j) += centered_diff(i) * centered_diff(j);
              }
          }
      }
      cov_matrix *= (1.0 / (indices.size() - 1.0));
      return cov_matrix;
  }
  
  template <typename T>
  void MomentumSmearing<T>::CovMatrixUncertainty(Int_t numberOfSamples)
  {
    if (_diffVectors.size() < 2)
    {
      std::cout << "ERROR: Not enough data points (" << _diffVectors.size() << ") to perform bootstrap (need at least 2)." << std::endl;
      return;
    }
    
    TRandom3 randGen(0); // For reproducibility

    std::vector<TMatrixT<T>> bootstrap_covs;
    const int n_original_size = _diffVectors.size();

    // Perform the bootstrap resampling
    for (Int_t i = 0; i < numberOfSamples; ++i)
    {
      std::vector<Int_t> resampled_indices;
      for (Int_t j = 0; j < n_original_size; ++j)
      {
        int index = randGen.Integer(n_original_size);
        resampled_indices.push_back(index);
      }
      bootstrap_covs.push_back(_calculateCovarianceFromSample(resampled_indices));
    }

    // --- Analyze the results to find mean and uncertainty ---
    TMatrixT<T> meanCovMatrix(_vecSize, _vecSize);
    meanCovMatrix.Zero();
    TMatrixT<T> uncertaintyMatrix(_vecSize, _vecSize);
    uncertaintyMatrix.Zero();

    // Sum all bootstrap matrices
    for (const auto &matrix : bootstrap_covs)
    {
      meanCovMatrix += matrix;
    }
    meanCovMatrix *= (1.0 / numberOfSamples);

    // Calculate standard deviation using loops (no Multiply)
    for (const auto &matrix : bootstrap_covs)
    {
      TMatrixT<T> diff = matrix;
      diff -= meanCovMatrix;
      for(Int_t i = 0; i < _vecSize; i++) {
          for(Int_t j = 0; j < _vecSize; j++) {
              uncertaintyMatrix(i, j) += diff(i, j) * diff(i, j);
          }
      }
    }
    uncertaintyMatrix *= (1.0 / (numberOfSamples - 1.0));
    
    // Take the square root of each element to get standard deviation
    for (Int_t i = 0; i < _vecSize; i++)
    {
      for (Int_t j = 0; j < _vecSize; j++)
      {
        uncertaintyMatrix(i, j) = std::sqrt(uncertaintyMatrix(i, j));
      }
    }
    
    std::cout << "\n--- Bootstrap Analysis ---" << std::endl;
    std::cout << "Mean Covariance Matrix (from " << numberOfSamples << " samples):" << std::endl;
    meanCovMatrix.Print();
    std::cout << "Uncertainty Matrix (Standard Deviation of elements):" << std::endl;
    uncertaintyMatrix.Print();
    std::cout << "--------------------------" << std::endl;

    _covMatrix = meanCovMatrix; // Update the main covariance matrix with the bootstrap mean
    
    // Zapisz macierz średnią i macierz niepewności do JSON
    _CovMatrixToJSON("bootstrapMeanCovMatrix");
    
    // Tymczasowo zastąp _covMatrix macierzą niepewności aby ją zapisać
    TMatrixT<T> tempMatrix = _covMatrix;
    _covMatrix = uncertaintyMatrix;
    _CovMatrixToJSON("bootstrapUncertaintyMatrix");
    _covMatrix = tempMatrix; // Przywróć oryginalną macierz
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
    
    // Manual matrix-vector multiplication
    for (Int_t i = 0; i < _vecSize; i++) {
      deltaMom(i) = 0;
      for (Int_t j = 0; j < _vecSize; j++) {
        deltaMom(i) += _UT(i, j) * gaussVec(j);
      }
    }
    
    _momVecMC += deltaMom;
  }

  template <typename T>
  void MomentumSmearing<T>::GetSmearedMomentum(TVectorT<T> &momVecMC)
  {
    momVecMC = _momVecMC;
  }

  template <typename T>
  void MomentumSmearing<T>::_CovMatrixToJSON(const std::string& varName)
  {
    std::string json = (std::string)TBufferJSON::ToJSON(&_covMatrix);
    
    // Usuwanie niepotrzebnych znaków formatujących
    std::string cleaned_json = json;
    
    // Usuwanie \n
    size_t pos = 0;
    while ((pos = cleaned_json.find("\\n", pos)) != std::string::npos) {
        cleaned_json.replace(pos, 2, "");
        pos += 0;
    }
    
    // Usuwanie nadmiarowych spacji
    pos = 0;
    while ((pos = cleaned_json.find("  ", pos)) != std::string::npos) {
        cleaned_json.replace(pos, 2, " ");
        pos += 1;
    }
    
    // Usuwanie \" (escaped quotes)
    pos = 0;
    while ((pos = cleaned_json.find("\\\"", pos)) != std::string::npos) {
        cleaned_json.replace(pos, 2, "\"");
        pos += 1;
    }
    
    // Usuwanie okalających cudzysłowów jeśli istnieją
    if (cleaned_json.front() == '"' && cleaned_json.back() == '"') {
        cleaned_json = cleaned_json.substr(1, cleaned_json.length() - 2);
    }
    
    // Parsowanie do JSON object zamiast string
    try {
        nlohmann::json jsonObj = nlohmann::json::parse(cleaned_json);
        properties["momSmearing"][varName] = jsonObj;
    } catch (const std::exception& e) {
        std::cerr << "ERROR: Failed to parse JSON: " << e.what() << std::endl;
        std::cerr << "Raw JSON: " << json << std::endl;
        std::cerr << "Cleaned JSON: " << cleaned_json << std::endl;
        // Fallback - save as string if parsing fails
        properties["momSmearing"][varName] = cleaned_json;
    }

    std::ofstream outfile(propName); // Assuming propName is defined in your project
    if (outfile.is_open()) {
        outfile << properties.dump(4);
        outfile.close();
        std::cout << "Covariance matrix saved to JSON as '" << varName << "'" << std::endl;
    } else {
        std::cerr << "ERROR: Unable to open JSON file for writing." << std::endl;
    }
  }

  template <typename T>
  void MomentumSmearing<T>::_CholeskyDecomposition()
  {
    TDecompChol decomposition(_covMatrix);
    if (!decomposition.Decompose()) {
        std::cerr << "ERROR: Cholesky decomposition failed. Matrix is not positive-definite." << std::endl;
        return;
    }
    _U.ResizeTo(decomposition.GetU());
    _UT.ResizeTo(decomposition.GetU());
    
    _U = decomposition.GetU();
    _UT.Transpose(_U);
  }
}