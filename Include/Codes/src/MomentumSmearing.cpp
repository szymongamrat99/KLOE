#include <MomentumSmearing.h>
#include <const.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLine.h>
#include <TLegend.h>
#include <TText.h>

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

    // Plot bootstrap distributions for each matrix element
    _plotBootstrapDistributions(bootstrap_covs, meanCovMatrix, uncertaintyMatrix);

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

  template <typename T>
  void MomentumSmearing<T>::_plotBootstrapDistributions(const std::vector<TMatrixT<T>>& bootstrap_covs, const TMatrixT<T>& meanCovMatrix, const TMatrixT<T>& uncertaintyMatrix)
  {
    std::cout << "Creating bootstrap distribution plots..." << std::endl;
    
    // Labels for momentum components (can be customized)
    std::vector<std::string> componentLabels;
    for (Int_t i = 0; i < _vecSize; i++) {
      componentLabels.push_back("p" + std::to_string(i));
    }
    
    // Set ROOT plotting style
    gStyle->SetOptStat(1111); // Show statistics box
    gStyle->SetOptTitle(1);   // Show titles
    
    // Create histograms and plots for each matrix element
    for (Int_t i = 0; i < _vecSize; i++) {
      for (Int_t j = 0; j <= i; j++) { // Only lower triangular + diagonal (covariance matrix is symmetric)
        
        // Create histogram for this matrix element
        std::string histName = Form("h_cov_%d_%d", i, j);
        std::string histTitle = Form("Bootstrap Distribution of Cov(%s,%s);Covariance Value;Frequency", 
                                   componentLabels[i].c_str(), componentLabels[j].c_str());
        
        // Find min and max values for binning
        T minVal = bootstrap_covs[0](i, j);
        T maxVal = bootstrap_covs[0](i, j);
        for (const auto& matrix : bootstrap_covs) {
          T val = matrix(i, j);
          if (val < minVal) minVal = val;
          if (val > maxVal) maxVal = val;
        }
        
        // Add some margin to the range
        T range = maxVal - minVal;
        if (range == 0) range = std::abs(maxVal * 0.1); // Handle case where all values are the same
        minVal -= range * 0.1;
        maxVal += range * 0.1;
        
        // Create histogram
        TH1D* hist = new TH1D(histName.c_str(), histTitle.c_str(), 50, minVal, maxVal);
        hist->SetLineColor(kBlue);
        hist->SetFillColor(kCyan);
        hist->SetFillStyle(3001);
        
        // Fill histogram with bootstrap values
        for (const auto& matrix : bootstrap_covs) {
          hist->Fill(matrix(i, j));
        }
        
        // Create canvas
        std::string canvasName = Form("c_cov_%d_%d", i, j);
        TCanvas* canvas = new TCanvas(canvasName.c_str(), histTitle.c_str(), 800, 600);
        canvas->cd();
        
        // Draw histogram
        hist->Draw("HIST");
        
        // Add vertical lines for mean and ±1σ
        T mean = meanCovMatrix(i, j);
        T sigma = uncertaintyMatrix(i, j);
        
        TLine* meanLine = new TLine(mean, 0, mean, hist->GetMaximum());
        meanLine->SetLineColor(kRed);
        meanLine->SetLineWidth(2);
        meanLine->SetLineStyle(1);
        meanLine->Draw("same");
        
        TLine* sigmaLow = new TLine(mean - sigma, 0, mean - sigma, hist->GetMaximum());
        sigmaLow->SetLineColor(kGreen);
        sigmaLow->SetLineWidth(2);
        sigmaLow->SetLineStyle(2);
        sigmaLow->Draw("same");
        
        TLine* sigmaHigh = new TLine(mean + sigma, 0, mean + sigma, hist->GetMaximum());
        sigmaHigh->SetLineColor(kGreen);
        sigmaHigh->SetLineWidth(2);
        sigmaHigh->SetLineStyle(2);
        sigmaHigh->Draw("same");
        
        // Add legend
        TLegend* legend = new TLegend(0.65, 0.75, 0.95, 0.95);
        legend->AddEntry(hist, "Bootstrap samples", "f");
        legend->AddEntry(meanLine, Form("Mean: %.4e", (double)mean), "l");
        legend->AddEntry(sigmaLow, Form("#pm1#sigma: %.4e", (double)sigma), "l");
        legend->Draw();
        
        // Add text with statistics
        TText* text = new TText();
        text->SetNDC(true);
        text->SetTextSize(0.03);
        text->DrawText(0.15, 0.85, Form("Samples: %d", (int)bootstrap_covs.size()));
        text->DrawText(0.15, 0.80, Form("Mean: %.6e", (double)mean));
        text->DrawText(0.15, 0.75, Form("Std Dev: %.6e", (double)sigma));
        text->DrawText(0.15, 0.70, Form("RMS/Mean: %.3f%%", (double)(sigma/std::abs(mean)*100)));
        
        // Save plot
        std::string filename = Form("bootstrap_cov_%s_%s.png", 
                                  componentLabels[i].c_str(), componentLabels[j].c_str());
        canvas->SaveAs(filename.c_str());
        
        std::cout << "Saved bootstrap distribution plot: " << filename << std::endl;
        
        // Clean up
        delete canvas;
        delete hist;
        delete meanLine;
        delete sigmaLow;
        delete sigmaHigh;
        delete legend;
        delete text;
      }
    }
    
    std::cout << "Bootstrap distribution plots completed." << std::endl;
  }
}