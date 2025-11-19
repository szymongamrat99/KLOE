#include <boost/filesystem.hpp>
#include <cctype>
#include <const.h>
#include <ctime>
#include <iostream>
#include <kloe_class.h>
#include <string>
#include <vector>

namespace KLOE {
pm00::pm00(TLorentzVector *mom_list, TLorentzVector *pos_list) {
  phi_mom = mom_list[0];
  phi_pos = pos_list[0];

  kaon_mom[0] = mom_list[1];
  kaon_pos[0] = pos_list[1];

  kaon_mom[1] = mom_list[2];
  kaon_pos[1] = pos_list[2];

  pi_ch_mom[0] = mom_list[3];
  pi_ch_pos[0] = pos_list[3];

  pi_ch_mom[1] = mom_list[4];
  pi_ch_pos[1] = pos_list[4];

  pi_ne_mom[0] = mom_list[5];
  pi_ne_pos[0] = pos_list[5];

  pi_ne_mom[1] = mom_list[6];
  pi_ne_pos[1] = pos_list[6];

  photon_mom[0] = mom_list[7];
  photon_pos[0] = pos_list[7];

  photon_mom[1] = mom_list[8];
  photon_pos[1] = pos_list[8];

  photon_mom[2] = mom_list[9];
  photon_pos[2] = pos_list[9];

  photon_mom[3] = mom_list[10];
  photon_pos[3] = pos_list[10];
};

pm00::pm00() : data(NULL), mc_sum(NULL) {
  phi_mom.SetPxPyPzE(0., 0., 0., 0.);
  phi_pos.SetXYZT(0., 0., 0., 0.);

  kaon_mom[0].SetPxPyPzE(0., 0., 0., 0.);
  kaon_pos[0].SetXYZT(0., 0., 0., 0.);

  kaon_mom[1].SetPxPyPzE(0., 0., 0., 0.);
  kaon_pos[1].SetXYZT(0., 0., 0., 0.);

  pi_ch_mom[0].SetPxPyPzE(0., 0., 0., 0.);
  pi_ch_pos[0].SetXYZT(0., 0., 0., 0.);

  pi_ch_mom[1].SetPxPyPzE(0., 0., 0., 0.);
  pi_ch_pos[1].SetXYZT(0., 0., 0., 0.);

  pi_ne_mom[0].SetPxPyPzE(0., 0., 0., 0.);
  pi_ne_pos[0].SetXYZT(0., 0., 0., 0.);

  pi_ne_mom[1].SetPxPyPzE(0., 0., 0., 0.);
  pi_ne_pos[1].SetXYZT(0., 0., 0., 0.);

  photon_mom[0].SetPxPyPzE(0., 0., 0., 0.);
  photon_pos[0].SetXYZT(0., 0., 0., 0.);

  photon_mom[1].SetPxPyPzE(0., 0., 0., 0.);
  photon_pos[1].SetXYZT(0., 0., 0., 0.);

  photon_mom[2].SetPxPyPzE(0., 0., 0., 0.);
  photon_pos[2].SetXYZT(0., 0., 0., 0.);

  photon_mom[3].SetPxPyPzE(0., 0., 0., 0.);
  photon_pos[3].SetXYZT(0., 0., 0., 0.);

  _MClistPath1 = Paths::path_tmp;
  _MClistPath2 = Paths::path_tmp;
  _MClistPath3 = Paths::path_tmp;

  _DatalistPath = Paths::path_tmp;
}

void pm00::inv_mass_calc(TLorentzVector four_mom) {
  inv_mass = four_mom.Mag();
};

void pm00::angle(TLorentzVector vec1, TLorentzVector vec2) {
  angle_vec = vec1.Angle(vec2.Vect());
};

void pm00::cyl_comp(TLorentzVector vec) {
  transv = vec.Perp();
  azim_angle = vec.Phi();
};

void pm00::boost_vector(TLorentzVector four_mom) {
  boost = four_mom.BoostVector();
};

void pm00::lorentz_transf(TLorentzVector four_mom) { four_mom.Boost(-boost); };

void pm00::lorentz_transf(Double_t *boost_vec, Double_t *vec_init,
                          Double_t *vec_end) {
  TLorentzVector before(vec_init);
  TVector3 boost(boost_vec);

  before.Boost(boost);

  for (Int_t i = 0; i < 4; i++)
    vec_end[i] = before(i);
};

void pm00::lorentz_transf(TVector3 &boost_vec, TLorentzVector &vec_init,
                          TLorentzVector &vec_end) const {
  TLorentzVector before(vec_init);
  TVector3 boost(boost_vec);

  before.Boost(boost);

  vec_end = before;
};

void pm00::lorentz_transf(Float_t *boost_vec, Float_t *vec_init,
                          Float_t *vec_end) {
  TLorentzVector before(vec_init);
  TVector3 boost(boost_vec);

  before.Boost(boost);

  for (Int_t i = 0; i < 4; i++)
    vec_end[i] = before(i);
};

Double_t pm00::DeltaT(TLorentzVector *momKch, TLorentzVector *posKch,
                      TLorentzVector *momKne, TLorentzVector *posKne) {
  TVector3 boostVKch = momKch->BoostVector(), boostVKne = momKne->BoostVector();

  posKch->Boost(-boostVKch);
  posKne->Boost(-boostVKne);

  Double_t DeltaT = (posKch->T() - posKne->T()) /
                    (PhysicsConstants::tau_S_nonCPT * PhysicsConstants::cVel);

  return DeltaT;
}

std::map<Int_t, Int_t> pm00::CountRepeatingElements(std::vector<Int_t> &arr) {
  std::map<Int_t, Int_t> frequencyMap;

  for (Int_t num : arr) {
    frequencyMap[num]++;
  }

  return frequencyMap;
}

Int_t pm00::signum(Float_t value) { return value / abs(value); }

void pm00::Clear1DArray(UInt_t M, Float_t *array) {
  for (UInt_t i = 0; i < M; i++)
    array[i] = 999.;
};

void pm00::Clear2DArray(UInt_t M, UInt_t N, Float_t **array) {
  for (UInt_t i = 0; i < M; i++)
    for (UInt_t j = 0; j < N; j++)
      array[i][j] = 999.;
};

void pm00::Clear1DArray(UInt_t M, Int_t *array) {
  for (UInt_t i = 0; i < M; i++)
    array[i] = -1;
};

void pm00::Clear2DArray(UInt_t M, UInt_t N, Int_t **array) {
  for (UInt_t i = 0; i < M; i++)
    for (UInt_t j = 0; j < N; j++)
      array[i][j] = -1;
};

// Funkcja porównująca dwie nieposortowane tablice z użyciem std::map
Int_t pm00::ArrayEquality(const Int_t *array1, const Int_t *array2,
                          const Int_t &size) {
  std::map<Int_t, Int_t> numerator1, numerator2;

  // Zliczanie elementów w pierwszej tablicy
  for (Int_t i = 0; i < size; ++i) {
    ++numerator1[array1[i]];
  }

  // Zliczanie elementów w drugiej tablicy
  for (Int_t i = 0; i < size; ++i) {
    ++numerator2[array2[i]];
  }

  // Porównanie liczników i obliczenie różnic
  Int_t numberOfDifferences = 0;

  // Iteracja po elementach pierwszej mapy
  for (auto it = numerator1.begin(); it != numerator1.end(); ++it) {
    Int_t key = it->first;
    Int_t value1 = it->second;
    Int_t value2 =
        numerator2[key]; // Liczba wystąpień w drugiej tablicy (domyślnie 0)
    numberOfDifferences +=
        std::abs(value1 - value2); // Dodaj różnicę w liczbie wystąpień
  }

  std::cout << numberOfDifferences << std::endl;

  // Iteracja po elementach drugiej mapy, aby znaleźć keye, których nie było w
  // pierwszej tablicy for (auto it = numerator2.begin(); it !=
  // numerator2.end(); ++it)
  // {
  //     Int_t key = it->first;
  //     if (numerator1.find(key) == numerator1.end())
  //     { // Jeśli keya nie ma w pierwszej mapie
  //         numberOfDifferences += it->second;
  //     }
  // }

  return numberOfDifferences;
}

Int_t pm00::ListOfFiles(const std::string &mode, Int_t &maxNumOfFiles,
                        std::vector<std::string> &fileNames) {
  std::ifstream listFiles;

  try {
    if (mode == "MC1")
      listFiles.open(_MClistPath1);
    else if (mode == "MC2")
      listFiles.open(_MClistPath2);
    else if (mode == "MC3")
      listFiles.open(_MClistPath3);
    else if (mode == "Data")
      listFiles.open(_DatalistPath);
    else
      throw ErrorHandling::ErrorCodes::FILE_NOT_EXIST;
  } catch (ErrorHandling::ErrorCodes err) {
    return Int_t(err);
  }

  std::string line = "";

  try {
    if (!listFiles)
      throw ErrorHandling::ErrorCodes::FILE_NOT_EXIST;

    Int_t numOfFiles = 0;

    while (std::getline(listFiles, line) && numOfFiles < maxNumOfFiles) {
      fileNames.push_back(line);

      numOfFiles++;
    }
  } catch (ErrorHandling::ErrorCodes err) {
    return Int_t(err);
  }

  return 0;
};

Int_t pm00::neu_triangle(Float_t *TrcSumFinal, Float_t *vtxSigmaFinal,
                         Float_t Clu5Vec[4][5], Float_t *ip, Float_t *Phi4Mom,
                         Float_t *Kne4Mom, Float_t *Kne4Vec,
                         Float_t *trc) const {
  // Logging of errors
  std::string logFilename = Paths::logsCNAFDir + "TriangleIter.log";
  ErrorHandling::ErrorLogs logger(logFilename);
  // ----------------------------------------------------------------------------

  Float_t KneTotMom = 0., BetaK, CosPkD[4], DTot[4];
  Float_t A, B[4], C[4], Delta[4], lK[4][2], lGamma[4][2], lGammaFinal[4],
      trctmp[4][2], lKtrue[4];

  Float_t NeuVtxTmp[4][2][4], NeuVtxTrueClu[4][4],
      NeuVtxAvg[4] = {0.}, EneTot = 0., TrcSum = 0., vtxSigma = 0.;

  TVector3 D[4], pK(Kne4Mom[0], Kne4Mom[1], Kne4Mom[2]);

  TrcSum = 0.;
  vtxSigma = 0.;

  KneTotMom = pK.Mag();
  BetaK = KneTotMom / Kne4Mom[3];

  A = (1 - pow(BetaK, 2)) / pow(BetaK, 2);

  for (Int_t i = 0; i < 4; i++) {
    for (Int_t j = 0; j < 3; j++) {
      D[i](j) = Clu5Vec[i][j] - ip[j];
    }

    DTot[i] = D[i].Mag();

    CosPkD[i] = (D[i].Dot(pK)) / (DTot[i] * KneTotMom);

    B[i] = 2 * (DTot[i] * CosPkD[i] -
                (PhysicsConstants::cVel * Clu5Vec[i][3] / BetaK));
    C[i] = pow(PhysicsConstants::cVel * Clu5Vec[i][3], 2) - pow(DTot[i], 2);

    Delta[i] = pow(B[i], 2) - 4 * A * C[i];

    try {
      if (Delta[i] < 0.) {
        throw ErrorHandling::ErrorCodes::DELTA_LT_ZERO;
      } else if (A == 0.) {
        throw ErrorHandling::ErrorCodes::DENOM_EQ_ZERO;
      } else {
        lK[i][0] = (-B[i] - sqrt(Delta[i])) / (2. * A);
        lK[i][1] = (-B[i] + sqrt(Delta[i])) / (2. * A);

        for (Int_t j = 0; j < 2; j++) {
          for (Int_t k = 0; k < 3; k++)
            NeuVtxTmp[i][j][k] = ip[k] + (lK[i][j] * (pK(k) / KneTotMom));

          NeuVtxTmp[i][j][3] = (lK[i][j] / (PhysicsConstants::cVel * BetaK));

          lGamma[i][j] = sqrt(pow(Clu5Vec[i][0] - NeuVtxTmp[i][j][0], 2) +
                              pow(Clu5Vec[i][1] - NeuVtxTmp[i][j][1], 2) +
                              pow(Clu5Vec[i][2] - NeuVtxTmp[i][j][2], 2));

          trctmp[i][j] = Clu5Vec[i][3] - NeuVtxTmp[i][j][3] -
                         (lGamma[i][j] / PhysicsConstants::cVel);
        }

        if (abs(trctmp[i][0]) < abs(trctmp[i][1])) {
          for (Int_t j = 0; j < 4; j++)
            NeuVtxTrueClu[i][j] = NeuVtxTmp[i][0][j];
        } else {
          for (Int_t j = 0; j < 4; j++)
            NeuVtxTrueClu[i][j] = NeuVtxTmp[i][1][j];
        }

        for (Int_t j = 0; j < 3; j++)
          NeuVtxAvg[j] += Clu5Vec[i][4] * NeuVtxTrueClu[i][j];

        EneTot += Clu5Vec[i][4];
      }
    } catch (ErrorHandling::ErrorCodes err) {
      logger.getErrLog(err);

      return int(err);
    }
  }

  for (Int_t j = 0; j < 3; j++)
    NeuVtxAvg[j] = NeuVtxAvg[j] / EneTot;

  NeuVtxAvg[3] =
      sqrt(pow(NeuVtxAvg[0] - ip[0], 2) + pow(NeuVtxAvg[1] - ip[1], 2) +
           pow(NeuVtxAvg[2] - ip[2], 2)) /
      (BetaK * PhysicsConstants::cVel);

  for (Int_t i = 0; i < 4; i++) {
    lGammaFinal[i] = sqrt(pow(Clu5Vec[i][0] - NeuVtxAvg[0], 2) +
                          pow(Clu5Vec[i][1] - NeuVtxAvg[1], 2) +
                          pow(Clu5Vec[i][2] - NeuVtxAvg[2], 2));
    trc[i] = Clu5Vec[i][3] - NeuVtxAvg[3] -
             (lGammaFinal[i] / PhysicsConstants::cVel);

    vtxSigma += Clu5Vec[i][4] *
                sqrt(pow(NeuVtxTrueClu[i][0] - NeuVtxAvg[0], 2) +
                     pow(NeuVtxTrueClu[i][1] - NeuVtxAvg[1], 2) +
                     pow(NeuVtxTrueClu[i][2] - NeuVtxAvg[2], 2)) /
                EneTot;
  }

  Kne4Vec[0] = NeuVtxAvg[0];
  Kne4Vec[1] = NeuVtxAvg[1];
  Kne4Vec[2] = NeuVtxAvg[2];
  Kne4Vec[3] = NeuVtxAvg[3];

  TrcSum = trc[0] + trc[1] + trc[2] + trc[3];

  *TrcSumFinal = TrcSum;
  *vtxSigmaFinal = vtxSigma;

  return 0;
};

Int_t pm00::neu_triangle(std::vector<TLorentzVector> *Clu4Mom,
                         std::vector<TLorentzVector> *Clu4Vec, Float_t *ip,
                         TLorentzVector *Phi4Mom, TLorentzVector *Kne4Mom,
                         TLorentzVector *Kne4Vec, Float_t *trc,
                         Float_t *TrcSumFinal, Float_t *vtxSigmaFinal) const {
  // Logging of errors

  Int_t numOfPhotons = Clu4Vec->size();

  // Logging of errors
  std::string logFilename = Paths::logsCNAFDir + "TriangleIter" +
                            std::to_string(numOfPhotons) + ".log";
  ErrorHandling::ErrorLogs logger(logFilename);
  // ----------------------------------------------------------------------------

  Float_t KneTotMom = 0., BetaK, CosPkD[4], DTot[4];
  Float_t A, B[4], C[4], Delta[4], lK[4][2], lGamma[4][2], lGammaFinal[4],
      trctmp[4][2], lKtrue[4];

  Float_t NeuVtxTmp[4][2][4], NeuVtxTrueClu[4][4],
      NeuVtxAvg[4] = {0.}, EneTot = 0., TrcSum = 0., vtxSigma = 0.;

  TVector3 D[4], pK = Kne4Mom->Vect();

  TrcSum = 0.;
  vtxSigma = 0.;

  KneTotMom = pK.Mag();
  BetaK = KneTotMom / Kne4Mom->E();

  A = (1 - pow(BetaK, 2)) / pow(BetaK, 2);

  for (Int_t i = 0; i < Clu4Vec->size(); i++) {
    for (Int_t j = 0; j < 3; j++) {
      D[i](j) = (*Clu4Vec)[i][j] - ip[j];
    }

    DTot[i] = D[i].Mag();

    CosPkD[i] = (D[i].Dot(pK)) / (DTot[i] * KneTotMom);

    B[i] = 2 * (DTot[i] * CosPkD[i] -
                (PhysicsConstants::cVel * (*Clu4Vec)[i][3] / BetaK));
    C[i] = pow(PhysicsConstants::cVel * (*Clu4Vec)[i][3], 2) - pow(DTot[i], 2);

    Delta[i] = pow(B[i], 2) - 4 * A * C[i];

    try {
      if (Delta[i] < 0.) {
        throw ErrorHandling::ErrorCodes::DELTA_LT_ZERO;
      } else if (A == 0.) {
        throw ErrorHandling::ErrorCodes::DENOM_EQ_ZERO;
      } else {
        lK[i][0] = (-B[i] - sqrt(Delta[i])) / (2. * A);
        lK[i][1] = (-B[i] + sqrt(Delta[i])) / (2. * A);

        for (Int_t j = 0; j < 2; j++) {
          for (Int_t k = 0; k < 3; k++)
            NeuVtxTmp[i][j][k] = ip[k] + (lK[i][j] * (pK(k) / KneTotMom));

          NeuVtxTmp[i][j][3] = (lK[i][j] / (PhysicsConstants::cVel * BetaK));

          lGamma[i][j] = sqrt(pow((*Clu4Vec)[i][0] - NeuVtxTmp[i][j][0], 2) +
                              pow((*Clu4Vec)[i][1] - NeuVtxTmp[i][j][1], 2) +
                              pow((*Clu4Vec)[i][2] - NeuVtxTmp[i][j][2], 2));

          trctmp[i][j] = (*Clu4Vec)[i][3] - NeuVtxTmp[i][j][3] -
                         (lGamma[i][j] / PhysicsConstants::cVel);
        }

        if (abs(trctmp[i][0]) < abs(trctmp[i][1])) {
          for (Int_t j = 0; j < 4; j++)
            NeuVtxTrueClu[i][j] = NeuVtxTmp[i][0][j];
        } else {
          for (Int_t j = 0; j < 4; j++)
            NeuVtxTrueClu[i][j] = NeuVtxTmp[i][1][j];
        }

        for (Int_t j = 0; j < 3; j++)
          NeuVtxAvg[j] += (*Clu4Vec)[i][4] * NeuVtxTrueClu[i][j];

        EneTot += (*Clu4Mom)[i][3];
      }
    } catch (ErrorHandling::ErrorCodes err) {
      logger.getErrLog(err);

      return int(err);
    }
  }

  for (Int_t j = 0; j < 3; j++)
    NeuVtxAvg[j] = NeuVtxAvg[j] / EneTot;

  NeuVtxAvg[3] =
      sqrt(pow(NeuVtxAvg[0] - ip[0], 2) + pow(NeuVtxAvg[1] - ip[1], 2) +
           pow(NeuVtxAvg[2] - ip[2], 2)) /
      (BetaK * PhysicsConstants::cVel);

  for (Int_t i = 0; i < 4; i++) {
    lGammaFinal[i] = sqrt(pow((*Clu4Vec)[i][0] - NeuVtxAvg[0], 2) +
                          pow((*Clu4Vec)[i][1] - NeuVtxAvg[1], 2) +
                          pow((*Clu4Vec)[i][2] - NeuVtxAvg[2], 2));
    trc[i] = (*Clu4Vec)[i][3] - NeuVtxAvg[3] -
             (lGammaFinal[i] / PhysicsConstants::cVel);

    vtxSigma += (*Clu4Mom)[i][3] *
                sqrt(pow(NeuVtxTrueClu[i][0] - NeuVtxAvg[0], 2) +
                     pow(NeuVtxTrueClu[i][1] - NeuVtxAvg[1], 2) +
                     pow(NeuVtxTrueClu[i][2] - NeuVtxAvg[2], 2)) /
                EneTot;
  }

  (*Kne4Vec)[0] = NeuVtxAvg[0];
  (*Kne4Vec)[1] = NeuVtxAvg[1];
  (*Kne4Vec)[2] = NeuVtxAvg[2];
  (*Kne4Vec)[3] = NeuVtxAvg[3];

  *TrcSumFinal = trc[0] + trc[1] + trc[2] + trc[3];
  *vtxSigmaFinal = vtxSigma;

  return 0;
};

// chainInit is now handled by FileManager. Implementation removed from pm00 as
// part of modularization.

// Function for timer start
void pm00::startTimer() {
  _start_time = std::chrono::high_resolution_clock::now();
}
// -------------------------------------------------------------------
// Function for timer end
std::string pm00::endTimer() {
  _end_time = std::chrono::high_resolution_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::seconds>(_end_time - _start_time);

  return (std::string)Utils::elapsedTimeHMS(duration.count());
}
// -------------------------------------------------------------------

// Function for today timestamp
std::string pm00::getCurrentTimestamp() const {
  auto now = std::chrono::system_clock::now();
  std::time_t now_time = std::chrono::system_clock::to_time_t(now);
  std::tm *tm_info = std::localtime(&now_time);

  std::ostringstream oss;
  oss << std::put_time(tm_info, "%Y-%m-%d_%H%M");
  return oss.str();
}

// Function for today timestamp
std::string pm00::getCurrentDate() const {
  auto now = std::chrono::system_clock::now();
  std::time_t now_time = std::chrono::system_clock::to_time_t(now);
  std::tm *tm_info = std::localtime(&now_time);

  std::ostringstream oss;
  oss << std::put_time(tm_info, "%Y-%m-%d");
  return oss.str();
}

// Data flags inline function
void pm00::dataFlagSetter(Controls::DataType &dataType, bool &dataFlag,
                          int &mcflag, int &mctruth) {
  switch (dataType) {
  case Controls::DataType::SIGNAL_TOT: {
    dataFlag = (mcflag == 1 && (mctruth == 1 || mctruth == 2));
    break;
  }
  case Controls::DataType::SIG_BCG: {
    dataFlag = (mcflag == 1 && mctruth != 0);
    break;
  }
  case Controls::DataType::MC_DATA: {
    dataFlag = (mcflag == 0 || (mcflag == 1 && mctruth != 0));
    break;
  }
  case Controls::DataType::SIGNAL_MAX: {
    dataFlag = (mcflag == 1 && (mctruth == 1 || mctruth == 2 || mctruth == 0));
    break;
  }
  case Controls::DataType::MC_ONLY: {
    dataFlag = (mcflag == 1 && mctruth != 0);
    break;
  }
  case Controls::DataType::DATA_ONLY: {
    dataFlag = (mcflag == 0);
    break;
  }
  }
}

void pm00::ConditionInitializer(std::vector<TFormula> &formula_vector,
                                std::vector<TString> &formula_names) const {
  for (Int_t i = 0; i < formula_names.size(); i++)
    formula_vector[i] = TFormula(formula_names[i], "abs(x - [0]) > [1]");
}

void pm00::SetSingleConditionParameters(TFormula &formula,
                                        Double_t *mean_sigma) const {
  formula.SetParameters(mean_sigma);
}

void pm00::SetAllConditionParameters(
    std::vector<TFormula> &formula_vector,
    std::vector<Double_t *> &mean_sigma_vector) const {
  Int_t size_formula = formula_vector.size(),
        size_parameters = mean_sigma_vector.size();

  if (size_formula != size_parameters)
    std::cout << "Both vectors need to have the same size!" << std::endl;
  else {
    for (Int_t i = 0; i < size_formula; i++) {
      formula_vector[i].SetParameters(mean_sigma_vector[i]);
    }
  }
}

Bool_t pm00::GetSingleConditionValue(TFormula &formula, Double_t &x) const {
  return formula.Eval(x);
}

void pm00::GetAllConditionValue(std::vector<TFormula> &formula_vector,
                                std::vector<Double_t> &x,
                                std::vector<Bool_t> &values) const {
  Int_t size_formula = formula_vector.size(), size_x = x.size();

  if (size_formula != size_x)
    std::cout << "Both vectors need to have the same size!" << std::endl;
  else {
    for (Int_t i = 0; i < size_formula; i++) {
      values.push_back(formula_vector[i].Eval(x[i]));
    }
  }
}

std::string pm00::CreateDatedFolder(const std::string &basePath) {
  std::string dateStr = pm00::getCurrentDate();
  boost::filesystem::path folderPath =
      boost::filesystem::path(basePath) / dateStr;

  if (!boost::filesystem::exists(folderPath)) {
    boost::filesystem::create_directories(folderPath);
  }

  return folderPath.string();
}

Double_t pm00::WeightedAverageSymmetric(std::vector<Double_t> &values,
                                        std::vector<Double_t> &errors) {
  Double_t sumWeights = 0;
  Double_t weightedSum = 0;

  for (size_t i = 0; i < values.size(); ++i) {
    Double_t weight = 1.0 / (errors[i] * errors[i]);
    sumWeights += weight;
    weightedSum += values[i] * weight;
  }

  return weightedSum / sumWeights;
}

Double_t pm00::WAvgSymmError(std::vector<Double_t> &errors) {
  Double_t sumWeights = 0;

  for (size_t i = 0; i < errors.size(); ++i) {
    Double_t weight = 1.0 / (errors[i] * errors[i]);
    sumWeights += weight;
  }

  return sqrt(1.0 / sumWeights);
}

Double_t pm00::WeightedAverageAsymmetric(std::vector<Double_t> &values,
                                         std::vector<Double_t> &errorsLow,
                                         std::vector<Double_t> &errorsUp) {
  Double_t sumWeights = 0;
  Double_t weightedSum = 0;

  for (size_t i = 0; i < values.size(); ++i) {
    Double_t error = 0.5 * (errorsLow[i] + errorsUp[i]);
    Double_t weight = 1.0 / (error * error);

    sumWeights += weight;
    weightedSum += values[i] * weight;
  }

  return weightedSum / sumWeights;
}

Double_t pm00::WAvgAsymmError(std::vector<Double_t> &errorsLow,
                              std::vector<Double_t> &errorsUp) {
  Double_t sumWeights = 0;

  for (size_t i = 0; i < errorsLow.size(); ++i) {
    Double_t error = 0.5 * (errorsLow[i] + errorsUp[i]);
    Double_t weight = 1.0 / (error * error);

    sumWeights += weight;
  }

  return sqrt(1.0 / sumWeights);
}

Double_t pm00::TwoBodyDecayMass(Double_t M, Double_t m1, Double_t m2) {
  // M - mass of the decaying particle
  // m1, m2 - masses of the decay products

  if (M < (m1 + m2)) {
    std::cout << "Mass of the decaying particle is less than the sum of the "
                 "decay products' masses!"
              << std::endl;
    return 0.;
  }

  Double_t lambda = (pow(M, 4) + pow(m1, 4) + pow(m2, 4) -
                     2 * pow(M, 2) * pow(m1, 2) - 2 * pow(M, 2) * pow(m2, 2) -
                     2 * pow(m1, 2) * pow(m2, 2)),
           p = sqrt(lambda) / (2 * M);

  return p;
}

/**
 * @brief Extract unique run numbers from filenames in a directory using regex
 * pattern.
 * @param directory Path to the directory to search.
 * @param regex_pattern Regex pattern with a capturing group for the run number
 * (e.g. R"(_(\\d+)\\.ntu$)").
 * @return std::vector<Int_t> of unique run numbers (sorted ascending).
 */
std::set<Int_t>
KLOE::pm00::ExtractUniqueRuns(const std::string &directory,
                              const std::string &regex_pattern) {
  std::set<Int_t> unique_runs; // Set to store unique run numbers
  std::regex run_regex(
      regex_pattern); // Regex to extract run number from filename
  boost::filesystem::path dir_path(directory); // Directory path
  if (!boost::filesystem::exists(dir_path) ||
      !boost::filesystem::is_directory(dir_path)) {
    // Return empty set if directory does not exist or is not a directory
    return {};
  }
  // Iterate over all files in the directory
  for (auto &entry : boost::filesystem::directory_iterator(dir_path)) {
    if (!boost::filesystem::is_regular_file(entry.status()))
      continue; // Skip non-regular files
    std::string filename =
        entry.path().filename().string(); // Get filename as string
    std::smatch match;
    // Search for run number using regex
    if (std::regex_search(filename, match, run_regex)) {
      if (match.size() > 1) {
        try {
          Int_t run =
              std::stoi(match[1].str()); // Convert matched string to integer
          unique_runs.insert(run);       // Insert run number into set
        } catch (...) {
          // Ignore parse errors
        }
      }
    }
  }
  // Return set of unique run numbers
  return std::set<Int_t>(unique_runs.begin(), unique_runs.end());
}

void pm00::CorrectClusterTime(Float_t T0,
                              std::vector<Float_t> &cluTimeCorrected) {
  for (size_t i = 0; i < cluTimeCorrected.size(); ++i) {
    cluTimeCorrected[i] = cluTimeCorrected[i] - T0;
  }
}

ErrorHandling::ErrorCodes pm00::triangleReconstruction(
    std::vector<Int_t> g4taken_kinfit, std::vector<Float_t> cluster[5],
    std::vector<Int_t> Asscl, std::vector<Float_t> bhabha_mom,
    std::vector<Float_t> Kchboost, std::vector<Float_t> ip,
    std::vector<Float_t> &Knetriangle, std::vector<Float_t> gammatriangle[4],
    Float_t &minv4gam, std::vector<Float_t> &trcfinal,
    ErrorHandling::ErrorLogs &logger) {

  Int_t ind_gam[4];

  Bool_t cond_ene;
  Bool_t cond_clus[4];
  Float_t Clu5Vec[4][5];

  Int_t done = 0;

  ind_gam[0] = g4taken_kinfit[0];
  ind_gam[1] = g4taken_kinfit[1];
  ind_gam[2] = g4taken_kinfit[2];
  ind_gam[3] = g4taken_kinfit[3];

  cond_ene = cluster[4][Asscl[ind_gam[0]] - 1] > MIN_CLU_ENE &&
             cluster[4][Asscl[ind_gam[1]] - 1] > MIN_CLU_ENE &&
             cluster[4][Asscl[ind_gam[2]] - 1] > MIN_CLU_ENE &&
             cluster[4][Asscl[ind_gam[3]] - 1] > MIN_CLU_ENE;

  cond_clus[0] = cluster[3][Asscl[ind_gam[0]] - 1] > 0 &&
                 cluster[0][Asscl[ind_gam[0]] - 1] != 0 &&
                 cluster[1][Asscl[ind_gam[0]] - 1] != 0 &&
                 cluster[2][Asscl[ind_gam[0]] - 1] != 0;
  cond_clus[1] = cluster[3][Asscl[ind_gam[1]] - 1] > 0 &&
                 cluster[0][Asscl[ind_gam[1]] - 1] != 0 &&
                 cluster[1][Asscl[ind_gam[1]] - 1] != 0 &&
                 cluster[2][Asscl[ind_gam[1]] - 1] != 0;
  cond_clus[2] = cluster[3][Asscl[ind_gam[2]] - 1] > 0 &&
                 cluster[0][Asscl[ind_gam[2]] - 1] != 0 &&
                 cluster[1][Asscl[ind_gam[2]] - 1] != 0 &&
                 cluster[2][Asscl[ind_gam[2]] - 1] != 0;
  cond_clus[3] = cluster[3][Asscl[ind_gam[3]] - 1] > 0 &&
                 cluster[0][Asscl[ind_gam[3]] - 1] != 0 &&
                 cluster[1][Asscl[ind_gam[3]] - 1] != 0 &&
                 cluster[2][Asscl[ind_gam[3]] - 1] != 0;

  if (cond_ene == true && cond_clus[0] && cond_clus[1] && cond_clus[2] &&
      cond_clus[3]) {
    for (Int_t k = 0; k < 4; k++) {
      Clu5Vec[k][0] = cluster[0][Asscl[ind_gam[k]] - 1];
      Clu5Vec[k][1] = cluster[1][Asscl[ind_gam[k]] - 1];
      Clu5Vec[k][2] = cluster[2][Asscl[ind_gam[k]] - 1];
      Clu5Vec[k][3] = cluster[3][Asscl[ind_gam[k]] - 1];
      Clu5Vec[k][4] = cluster[4][Asscl[ind_gam[k]] - 1];
    }

    //! Using the charged part of the decay

    Float_t Knerec[9] = {};

    Knerec[0] = bhabha_mom[0] - Kchboost[0];
    Knerec[1] = bhabha_mom[1] - Kchboost[1];
    Knerec[2] = bhabha_mom[2] - Kchboost[2];
    Knerec[3] = bhabha_mom[3] - Kchboost[3];

    //!

    Float_t TrcSum = 0., vtxSigma = 0., vtxSigmaMin = 1.e6, TrcSumMin = 1.e6,
            trcsum = 0.;

    Float_t neu_vtx[4], trc[4];

    // this->neu_triangle(&TrcSum, &vtxSigma, Clu5Vec, ip.data(),
    // bhabha_mom.data(), Knerec, neu_vtx, trc);

    if (sqrt(pow(vtxSigma, 2) + pow(TrcSum, 2)) <
        sqrt(pow(vtxSigmaMin, 2) + pow(TrcSumMin, 2))) {
      vtxSigmaMin = vtxSigma;
      TrcSumMin = TrcSum;

      trcsum = TrcSumMin;

      done = 1;

      for (Int_t l = 0; l < 4; l++) {
        Knetriangle[l] = Knerec[l];

        Knetriangle[6 + l] = neu_vtx[l];

        trcfinal[l] = trc[l];

        neutral_mom(cluster[0][Asscl[ind_gam[l]] - 1],
                    cluster[1][Asscl[ind_gam[l]] - 1],
                    cluster[2][Asscl[ind_gam[l]] - 1],
                    cluster[4][Asscl[ind_gam[l]] - 1], neu_vtx,
                    gammatriangle[l].data());
      }

      minv4gam = sqrt(pow(gammatriangle[0][3] + gammatriangle[1][3] +
                              gammatriangle[2][3] + gammatriangle[3][3],
                          2) -
                      pow(gammatriangle[0][0] + gammatriangle[1][0] +
                              gammatriangle[2][0] + gammatriangle[3][0],
                          2) -
                      pow(gammatriangle[0][1] + gammatriangle[1][1] +
                              gammatriangle[2][1] + gammatriangle[3][1],
                          2) -
                      pow(gammatriangle[0][2] + gammatriangle[1][2] +
                              gammatriangle[2][2] + gammatriangle[3][2],
                          2));

      Knetriangle[4] = 0.;
      for (Int_t l = 0; l < 3; l++) {
        Knetriangle[4] += pow(Knetriangle[l], 2);
      }

      Knetriangle[5] = sqrt(pow(Knetriangle[3], 2) - Knetriangle[4]);
      Knetriangle[4] = sqrt(Knetriangle[4]);
    }
  }

  if (done == 0) {
    return ErrorHandling::ErrorCodes::TRIANGLE_REC;
  }

  return ErrorHandling::ErrorCodes::NO_ERROR;
}

ErrorHandling::ErrorCodes
pm00::triangleReconstruction(std::vector<neutralParticle> &photon, phiMeson phi,
                             kaonNeutral Kchboost, Float_t *ip,
                             kaonNeutral &Knetriangle) {
  Bool_t cond_ene;
  Bool_t cond_clus[4];
  Float_t Clu5Vec[4][5];

  Int_t done = 0;

  cond_ene = photon[0].clusterParams[4] >= MIN_CLU_ENE &&
             photon[1].clusterParams[4] >= MIN_CLU_ENE &&
             photon[2].clusterParams[4] >= MIN_CLU_ENE &&
             photon[3].clusterParams[4] >= MIN_CLU_ENE;

  cond_clus[0] = photon[0].clusterParams[0] != 0 &&
                 photon[0].clusterParams[1] != 0 &&
                 photon[0].clusterParams[2] != 0;
  cond_clus[1] = photon[1].clusterParams[0] != 0 &&
                 photon[1].clusterParams[1] != 0 &&
                 photon[1].clusterParams[2] != 0;
  cond_clus[2] = photon[2].clusterParams[0] != 0 &&
                 photon[2].clusterParams[1] != 0 &&
                 photon[2].clusterParams[2] != 0;
  cond_clus[3] = photon[3].clusterParams[0] != 0 &&
                 photon[3].clusterParams[1] != 0 &&
                 photon[3].clusterParams[2] != 0;

  if (cond_ene && cond_clus[0] && cond_clus[1] && cond_clus[2] &&
      cond_clus[3]) {
    for (Int_t k = 0; k < 4; k++) {
      Clu5Vec[k][0] = photon[k].clusterParams[0];
      Clu5Vec[k][1] = photon[k].clusterParams[1];
      Clu5Vec[k][2] = photon[k].clusterParams[2];
      Clu5Vec[k][3] = photon[k].clusterParams[3];
      Clu5Vec[k][4] = photon[k].clusterParams[4];
    }

    //! Using the charged part of the decay

    Float_t Knerec[9] = {};

    Knerec[0] = phi.fourMom[0] - Kchboost.fourMom[0];
    Knerec[1] = phi.fourMom[1] - Kchboost.fourMom[1];
    Knerec[2] = phi.fourMom[2] - Kchboost.fourMom[2];
    Knerec[3] = phi.fourMom[3] - Kchboost.fourMom[3];

    //!

    Float_t TrcSum = 0., vtxSigma = 0., vtxSigmaMin = 1.e6, TrcSumMin = 1.e6,
            trcsum = 0.;

    Float_t neu_vtx[4], trc[4];

    neu_triangle(&TrcSum, &vtxSigma, Clu5Vec, ip, phi.fourMom.data(), Knerec,
                 neu_vtx, trc);

    if (sqrt(pow(vtxSigma, 2) + pow(TrcSum, 2)) <
        sqrt(pow(vtxSigmaMin, 2) + pow(TrcSumMin, 2))) {
      vtxSigmaMin = vtxSigma;
      TrcSumMin = TrcSum;

      trcsum = TrcSumMin;

      done = 1;

      for (Int_t l = 0; l < 4; l++) {
        Knetriangle.fourMom[l] = Knerec[l];

        Knetriangle.fourPos[l] = neu_vtx[l];

        neutral_mom(photon[l].clusterParams[0], photon[l].clusterParams[1],
                    photon[l].clusterParams[2], photon[l].clusterParams[4],
                    neu_vtx, photon[l].fourMom.data());
      }

      Knetriangle.SetTotalVector();
    }
  }

  if (done == 0) {
    return ErrorHandling::ErrorCodes::TRIANGLE_REC;
  }

  return ErrorHandling::ErrorCodes::NO_ERROR;
}

ErrorHandling::ErrorCodes
pm00::NeuCluWrongCheck(std::vector<Int_t> neuclulist, phiMeson phi,
                       kaonNeutral Kchboost, std::vector<Float_t> Xcl,
                       std::vector<Float_t> Ycl, std::vector<Float_t> Zcl,
                       std::vector<Float_t> Tcl, std::vector<Float_t> Enecl) {
  TLorentzVector Kch4MomLAB, Phi4MomLAB, Kch4MomCM, Kch4MomCMKaon;
  Phi4MomLAB.SetXYZT(phi.fourMom[0], phi.fourMom[1], phi.fourMom[2],
                     phi.fourMom[3]);
  Kch4MomLAB.SetXYZT(Kchboost.fourMom[0], Kchboost.fourMom[1],
                     Kchboost.fourMom[2], Kchboost.fourMom[3]);

  TVector3 PhiBeta = -Phi4MomLAB.BoostVector();

  lorentz_transf(PhiBeta, Kch4MomLAB, Kch4MomCM);

  return ErrorHandling::ErrorCodes::NO_ERROR;
}

void pm00::PhotonPairingToPi0(std::vector<Float_t> *photonMom,
                              Int_t numOfPhotons,
                              std::vector<Int_t> &bestPairingIndex) {
  Int_t intPi0 = numOfPhotons / 2.;

  std::vector<Int_t> indices(numOfPhotons);

  std::iota(indices.begin(), indices.end(), 0);

  Float_t invMassDiffMin = 1.e6;

  do {
    std::vector<Float_t> px(intPi0, 0.), py(intPi0, 0.), pz(intPi0, 0.),
        E(intPi0, 0.), invMassDiff(intPi0, 0.);

    Float_t invMassDiffTot = 0.;

    for (Int_t i = 0; i < intPi0; i++) {
      px[i] = photonMom[indices[i * intPi0]][0] +
              photonMom[indices[i * intPi0 + 1]][0];
      py[i] = photonMom[indices[i * intPi0]][1] +
              photonMom[indices[i * intPi0 + 1]][1];
      pz[i] = photonMom[indices[i * intPi0]][2] +
              photonMom[indices[i * intPi0 + 1]][2];
      E[i] = photonMom[indices[i * intPi0]][3] +
             photonMom[indices[i * intPi0 + 1]][3];

      invMassDiff[i] =
          sqrt(pow(E[i], 2) - pow(px[i], 2) - pow(py[i], 2) - pow(pz[i], 2));

      invMassDiffTot += pow(invMassDiff[i] - PhysicsConstants::mPi0, 2);
    }

    invMassDiffTot = sqrt(invMassDiffTot);

    if (invMassDiffTot < invMassDiffMin) {
      invMassDiffMin = invMassDiffTot;
      bestPairingIndex = indices;
    }

  } while (std::next_permutation(indices.begin(), indices.end()));
}

KLOE::KaonProperTimes KLOE::pm00::CalculateKaonProperTimes(
    const std::vector<Float_t> &kaon1Mom, const std::vector<Float_t> &kaon1Pos,
    const std::vector<Float_t> &kaon2Mom, const std::vector<Float_t> &kaon2Pos,
    const std::vector<Float_t> &ipPos) {
  KaonProperTimes result;

  // Oblicz czasy dla pierwszego kaona
  CalculateSingleKaonTime(kaon1Mom, kaon1Pos, ipPos, result.kaon1TimeLAB,
                          result.kaon1TimeCM);

  // Oblicz czasy dla drugiego kaona
  CalculateSingleKaonTime(kaon2Mom, kaon2Pos, ipPos, result.kaon2TimeLAB,
                          result.kaon2TimeCM);

  // Oblicz różnice czasów
  result.deltaTimeLAB = result.kaon1TimeLAB - result.kaon2TimeLAB;
  result.deltaTimeCM = result.kaon1TimeCM - result.kaon2TimeCM;

  return result;
}

void KLOE::pm00::CalculateSingleKaonTime(const std::vector<Float_t> &kaonMom,
                                         const std::vector<Float_t> &kaonPos,
                                         const std::vector<Float_t> &ipPos,
                                         Double_t &timeLAB, Double_t &timeCM) {
  // Sprawdź rozmiary wektorów
  if (kaonMom.size() < 4 || kaonPos.size() < 3 || ipPos.size() < 3) {
    std::cerr << "ERROR: Invalid vector sizes in CalculateSingleKaonTime"
              << std::endl;
    timeLAB = -999.0;
    timeCM = -999.0;
    return;
  }

  // Oblicz pęd kaona w układzie LAB (znormalizowany przez energię)
  TVector3 kaonBetaLAB = {-kaonMom[0] / kaonMom[3], -kaonMom[1] / kaonMom[3],
                          -kaonMom[2] / kaonMom[3]};

  // Oblicz drogę kaona w układzie LAB
  Double_t kaonPathLAB =
      sqrt(pow(kaonPos[0] - ipPos[0], 2) + pow(kaonPos[1] - ipPos[1], 2) +
           pow(kaonPos[2] - ipPos[2], 2));

  // Prędkość kaona w układzie LAB
  Double_t kaonTotalBetaLAB = kaonBetaLAB.Mag();

  // Czas w układzie LAB
  if (kaonTotalBetaLAB > 0) {
    timeLAB = kaonPathLAB / kaonTotalBetaLAB;
  } else {
    timeLAB = -999.0;
    timeCM = -999.0;
    return;
  }

  // Przejście do układu CM kaona
  TLorentzVector kaon4VecLAB = {kaonPos[0] - ipPos[0], // cm
                                kaonPos[1] - ipPos[1], // cm
                                kaonPos[2] - ipPos[2], // cm
                                timeLAB};              // cm

  TLorentzVector kaon4VecKaonCM = {0., 0., 0., 0.};

  // Użyj istniejącej metody transformacji Lorentza
  lorentz_transf(kaonBetaLAB, kaon4VecLAB, kaon4VecKaonCM);

  // Czas własny w jednostkach τ_S
  timeCM = kaon4VecKaonCM.T() /
           (PhysicsConstants::cVel * PhysicsConstants::tau_S_nonCPT);
  timeLAB = timeLAB / (PhysicsConstants::cVel * PhysicsConstants::tau_S_nonCPT);
}
} // namespace KLOE

void KLOE::pm00::trilaterationReconstruction(TVectorD X, Double_t neuVtx[2][4],
                                             Bool_t neuVtxErr[2]) {
  Int_t selected[4] = {1, 2, 3, 4};

  Reconstructor R; // Reconstructor object
  Solution S;      // Solution struct

  // Setting clusters for a solution
  for (Int_t k = 0; k < 4; k++) {
    R.SetClu(k, X[k * 5], X[k * 5 + 1], X[k * 5 + 2], X[k * 5 + 3],
             X[k * 5 + 4]);

    R.SetClu(4, 0., 0., 0., 0., 0.);
    R.SetClu(5, 0., 0., 0., 0., 0.);
  }
  // -------------------------------

  S = R.MySolve(selected); // Filling up the structure

  for (Int_t i = 0; i < 2; i++) {
    neuVtxErr[i] = S.error[i];

    for (Int_t j = 0; j < 4; j++) {
      neuVtx[i][j] = S.sol[i][j];
    }
  }
}

// File Management Methods
bool KLOE::pm00::CreateCurrentLinks(const std::string &rootFilesDir) {
  try {
    // Sprawdź czy katalog root_files istnieje
    if (!boost::filesystem::exists(rootFilesDir)) {
      std::cerr << "Katalog " << rootFilesDir << " nie istnieje!" << std::endl;
      return false;
    }

    // Znajdź folder z najnowszą datą
    std::string latestDateFolder;
    std::time_t latestTime = 0;

    for (boost::filesystem::directory_iterator entry(rootFilesDir);
         entry != boost::filesystem::directory_iterator(); ++entry) {
      if (boost::filesystem::is_directory(entry->status())) {
        std::string folderName = entry->path().filename().string();
        if (IsValidDateFormat(folderName)) {
          std::time_t folderTime =
              boost::filesystem::last_write_time(entry->path());
          if (latestDateFolder.empty() || folderTime > latestTime) {
            latestDateFolder = folderName;
            latestTime = folderTime;
          }
        }
      }
    }

    if (latestDateFolder.empty()) {
      std::cerr << "Nie znaleziono folderów z datą w formacie YYYY-MM-DD!"
                << std::endl;
      return false;
    }

    // Utwórz folder current
    std::string currentPath = rootFilesDir + "/current";
    if (boost::filesystem::exists(currentPath)) {
      boost::filesystem::remove_all(currentPath);
    }
    boost::filesystem::create_directory(currentPath);

    // Ścieżka do najnowszego folderu
    std::string latestFolderPath = rootFilesDir + "/" + latestDateFolder;

    // Znajdź wszystkie pliki .root w najnowszym folderze i utwórz dla nich
    // symlinki
    int linkCount = 0;
    for (boost::filesystem::directory_iterator entry(latestFolderPath);
         entry != boost::filesystem::directory_iterator(); ++entry) {
      if (boost::filesystem::is_regular_file(entry->status())) {
        std::string fileName = entry->path().filename().string();
        // Sprawdź czy to plik .root
        if (fileName.length() > 5 &&
            fileName.substr(fileName.length() - 5) == ".root") {
          std::string sourcePath = entry->path().string();
          std::string linkPath = currentPath + "/" + fileName;

          boost::filesystem::create_symlink(sourcePath, linkPath);
          std::cout << "Utworzono link: " << fileName << " -> " << sourcePath
                    << std::endl;
          linkCount++;
        }
      }
    }

    if (linkCount == 0) {
      std::cout << "Ostrzeżenie: Nie znaleziono plików .root w folderze "
                << latestDateFolder << std::endl;
    }

    std::cout << "Pomyślnie utworzono folder current z linkami do "
              << latestDateFolder << std::endl;
    return true;
  } catch (const std::exception &e) {
    std::cerr << "Błąd podczas tworzenia linków: " << e.what() << std::endl;
    return false;
  }
}

bool KLOE::pm00::IsValidDateFormat(const std::string &folderName) {
  // Sprawdź czy string ma długość 10 (YYYY-MM-DD)
  if (folderName.length() != 10)
    return false;

  // Sprawdź czy ma myślniki w odpowiednich miejscach
  if (folderName[4] != '-' || folderName[7] != '-')
    return false;

  // Sprawdź czy pozostałe znaki to cyfry
  for (int i = 0; i < 10; i++) {
    if (i != 4 && i != 7) // Pomiń myślniki
    {
      if (!std::isdigit(folderName[i]))
        return false;
    }
  }

  // Dodatkowa walidacja miesięcy i dni
  int month = std::stoi(folderName.substr(5, 2));
  int day = std::stoi(folderName.substr(8, 2));

  if (month < 1 || month > 12)
    return false;

  if (day < 1 || day > 31)
    return false;

  return true;
}
