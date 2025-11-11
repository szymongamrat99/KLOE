#include <NeutralReconstruction.h>

namespace KLOE
{
  NeutralReconstruction::NeutralReconstruction(Int_t nPhotons) : _nPhotons(nPhotons), _nPions(nPhotons / 2.), ErrorLogs("NeutralReconstruction.log")
  {
    try
    {
      if (!_checkNPhotons())
      {
        throw std::invalid_argument("Number of photons must be a positive even integer.");
      }

      _Photons.resize(_nPhotons);
      _Pions.resize(_nPions);
      _Charged.resize(2);
    }
    catch (const std::invalid_argument &e)
    {
      std::cerr << "Error: " << e.what() << std::endl;
    }
  };

  // Method to perform photon pairing to pions
  void NeutralReconstruction::PhotonPairingToPi0()
  {
    _PhotonPairingCore();
  }

  // Overloaded method to get results directly from the code
  void NeutralReconstruction::PhotonPairingToPi0(std::vector<neutralParticle> &photons, std::vector<Int_t> &bestPairingIndex)
  {
    SetPhotonParameters(photons);

    _PhotonPairingCore();

    photons = _Photons;
    bestPairingIndex = _bestPairingIndex;
  }

  void NeutralReconstruction::PhotonPairingToPi0WithOmega(std::vector<neutralParticle> &photons, std::vector<chargedParticle> &pions, std::vector<Int_t> &bestPairingIndexNeutral, std::vector<Int_t> &bestPairingIndexOmega, neutralParticle &omega, std::vector<neutralParticle> &neutralPions)
  {
    SetPhotonParameters(photons);

    SetChargedParameters(pions);

    _PhotonPairingWithOmegaCore();

    _Pi0ReconstructionCore();
    _OmegaReconstructionCore(_bestPairingIndexOmega[0]);
    _omega.SetTotalVector();

    neutralPions[0] = _Pions[_bestPairingIndexOmega[0]];
    neutralPions[1] = _Pions[_bestPairingIndexOmega[1]];

    photons = _Photons;
    pions = _Charged;
    omega = _omega;
    bestPairingIndexNeutral = _bestPairingIndex;
    bestPairingIndexOmega = _bestPairingIndexOmega;
  }

  // Method to perform pi0 reconstruction from paired photons
  void NeutralReconstruction::Pi0Reconstruction()
  {
    _Pi0ReconstructionCore();
  }

  // Overloaded method to get results directly from the code
  void NeutralReconstruction::Pi0Reconstruction(std::vector<neutralParticle> &pions)
  {
    _Pi0ReconstructionCore();
    pions = _Pions;
  }

  // Private methods

  void NeutralReconstruction::_PhotonPairingCore()
  {
    if (_nPhotons != 4)
    {
      std::cerr << "ERROR: Photon pairing currently supports only 4 photons." << std::endl;
      return;
    }

    _invMassDiffMin = 1.e6;

    std::vector<Float_t>
        px(_nPions, 0.),
        py(_nPions, 0.),
        pz(_nPions, 0.),
        E(_nPions, 0.),
        invMass(_nPions, 0.);

    for (const auto &indices : _combinations)
    {
      std::array<Float_t, 2> invMass;
      Float_t chi2_pi0s = 0.0;

      for (Int_t i = 0; i < _nPions; i++)
      {
        Float_t px = _Photons[indices[i * 2]].fourMom[0] + _Photons[indices[i * 2 + 1]].fourMom[0];
        Float_t py = _Photons[indices[i * 2]].fourMom[1] + _Photons[indices[i * 2 + 1]].fourMom[1];
        Float_t pz = _Photons[indices[i * 2]].fourMom[2] + _Photons[indices[i * 2 + 1]].fourMom[2];
        Float_t E = _Photons[indices[i * 2]].fourMom[3] + _Photons[indices[i * 2 + 1]].fourMom[3];

        invMass[i] = sqrt(E * E - px * px - py * py - pz * pz);

        chi2_pi0s += pow(invMass[i] - PhysicsConstants::mPi0, 2);
      }

      chi2_pi0s = sqrt(chi2_pi0s);

      if (chi2_pi0s < _invMassDiffMin)
      {
        _invMassDiffMin = chi2_pi0s;
        _bestPairingIndex = std::vector<Int_t>(indices.begin(), indices.end());
      }
    }
  }

  void NeutralReconstruction::_PhotonPairingWithOmegaCore()
  {
    if (_nPhotons != 4)
    {
      std::cerr << "ERROR: Photon pairing currently supports only 4 photons." << std::endl;
      return;
    }

    Float_t mPi0Sigma = 8.5, // MeV
        mOmegaSigma = 30.0;  // MeV

    _invMassDiffMin = 1.e6;

    for (const auto &indices : _combinations)
    {
      std::array<Float_t, 2> invMass;
      Float_t chi2_pi0s = 0.0;

      for (Int_t i = 0; i < _nPions; i++)
      {
        Float_t px = _Photons[indices[i * 2]].fourMom[0] + _Photons[indices[i * 2 + 1]].fourMom[0];
        Float_t py = _Photons[indices[i * 2]].fourMom[1] + _Photons[indices[i * 2 + 1]].fourMom[1];
        Float_t pz = _Photons[indices[i * 2]].fourMom[2] + _Photons[indices[i * 2 + 1]].fourMom[2];
        Float_t E = _Photons[indices[i * 2]].fourMom[3] + _Photons[indices[i * 2 + 1]].fourMom[3];

        invMass[i] = sqrt(E * E - px * px - py * py - pz * pz);

        chi2_pi0s += pow((invMass[i] - PhysicsConstants::mPi0) / mPi0Sigma, 2);
      }

      _Pi0ReconstructionCore(std::vector<Int_t>(indices.begin(), indices.end()));

      for (Int_t pionIdx = 0; pionIdx < _nPions; pionIdx++)
      {
        _OmegaReconstructionCore(pionIdx);

        Float_t chi2_omega = pow((_omega.mass - PhysicsConstants::mOmega) / mOmegaSigma, 2);
        Float_t chi2_total = chi2_pi0s + chi2_omega;

        if (chi2_total < _invMassDiffMin)
        {
          _invMassDiffMin = chi2_total;
          _bestPairingIndex = std::vector<Int_t>(indices.begin(), indices.end());
          _bestPairingIndexOmega = {pionIdx, (pionIdx + 1) % 2};
        }
      }
    }
  }

  void NeutralReconstruction::_Pi0ReconstructionCore()
  {
    for (Int_t i = 0; i < _nPions; i++)
    {
      _Pions[i].fourMom[0] = _Photons[_bestPairingIndex[i * _nPions]].fourMom[0] + _Photons[_bestPairingIndex[i * _nPions + 1]].fourMom[0];
      _Pions[i].fourMom[1] = _Photons[_bestPairingIndex[i * _nPions]].fourMom[1] + _Photons[_bestPairingIndex[i * _nPions + 1]].fourMom[1];
      _Pions[i].fourMom[2] = _Photons[_bestPairingIndex[i * _nPions]].fourMom[2] + _Photons[_bestPairingIndex[i * _nPions + 1]].fourMom[2];
      _Pions[i].fourMom[3] = _Photons[_bestPairingIndex[i * _nPions]].fourMom[3] + _Photons[_bestPairingIndex[i * _nPions + 1]].fourMom[3];

      _Pions[i].CalculateMassFromFourMom();
      _Pions[i].CalculateTotalMomentumFromFourMom();

      _Pions[i].SetLorentzVectors();

      _Pions[i].SetTotalVector();
    }

    for (Int_t i = 0; i < _nPions; i++)
    {
      _Pions[i].openingAngle = _Pions[i].lorentzFourMom.Angle(_Pions[(i + 1) % _nPions].lorentzFourMom.Vect());
    }
  }

  void NeutralReconstruction::_Pi0ReconstructionCore(std::vector<Int_t> indices)
  {
    for (Int_t i = 0; i < _nPions; i++)
    {
      _Pions[i].fourMom[0] = _Photons[indices[i * _nPions]].fourMom[0] + _Photons[indices[i * _nPions + 1]].fourMom[0];
      _Pions[i].fourMom[1] = _Photons[indices[i * _nPions]].fourMom[1] + _Photons[indices[i * _nPions + 1]].fourMom[1];
      _Pions[i].fourMom[2] = _Photons[indices[i * _nPions]].fourMom[2] + _Photons[indices[i * _nPions + 1]].fourMom[2];
      _Pions[i].fourMom[3] = _Photons[indices[i * _nPions]].fourMom[3] + _Photons[indices[i * _nPions + 1]].fourMom[3];

      _Pions[i].CalculateMassFromFourMom();
      _Pions[i].CalculateTotalMomentumFromFourMom();

      _Pions[i].SetLorentzVectors();
    }

    for (Int_t i = 0; i < _nPions; i++)
    {
      _Pions[i].openingAngle = _Pions[i].lorentzFourMom.Angle(_Pions[(i + 1) % _nPions].lorentzFourMom.Vect());
    }
  }

  void NeutralReconstruction::_OmegaReconstructionCore(Int_t index)
  {
    for (Int_t i = 0; i < 4; i++)
    {
      _omega.fourMom[i] = _Pions[index].fourMom[i] +
                          _Charged[0].fourMom[i] +
                          _Charged[1].fourMom[i];
    }

    _omega.CalculateMassFromFourMom();
    _omega.CalculateTotalMomentumFromFourMom();
    _omega.SetLorentzVectors();
  }

  ErrorHandling::ErrorCodes NeutralReconstruction::ReconstructSixGammaVertex(const std::vector<Float_t> cluster[5], const std::vector<Int_t> &neu_clu_list, std::vector<Int_t> &bestIndices, Float_t &bestError, kaonNeutral &KnerecSix, std::vector<neutralParticle> &photonSix)
  {
    const Int_t
        photonNum = 6,
        neuCluSize = neu_clu_list.size();

    Float_t bestTotalError = 1e6;
    Bool_t solutionFound = false;

    std::vector<std::array<Float_t, 4>>
        partialSolutions;

    std::vector<Float_t>
        partialSolutionEnergy;

    // Implementation of the six-gamma vertex reconstruction algorithm
    if (neuCluSize < photonNum)
      return ErrorHandling::ErrorCodes::LESS_THAN_SIX_NEUTRAL_CLUSTERS;

    // Vector of combinations of photonNum clusters
    std::vector<std::array<Int_t, photonNum>> combinations;
    for (Int_t i1 = 0; i1 < neuCluSize - 5; i1++)
      for (Int_t i2 = i1 + 1; i2 < neuCluSize - 4; i2++)
        for (Int_t i3 = i2 + 1; i3 < neuCluSize - 3; i3++)
          for (Int_t i4 = i3 + 1; i4 < neuCluSize - 2; i4++)
            for (Int_t i5 = i4 + 1; i5 < neuCluSize - 1; i5++)
              for (Int_t i6 = i5 + 1; i6 < neuCluSize; i6++)
              {
                combinations.push_back({neu_clu_list[i1],
                                        neu_clu_list[i2],
                                        neu_clu_list[i3],
                                        neu_clu_list[i4],
                                        neu_clu_list[i5],
                                        neu_clu_list[i6]});
              }
    ////////////////////////////////////////////////////////////////////////////

    for (auto &combo : combinations)
    {
      Float_t totalEnergy = cluster[4][combo[0] - 1] + cluster[4][combo[1] - 1] +
                            cluster[4][combo[2] - 1] + cluster[4][combo[3] - 1] +
                            cluster[4][combo[4] - 1] + cluster[4][combo[5] - 1];

      Bool_t energyLimitPerCluster = cluster[4][combo[0] - 1] > MIN_CLU_ENE &&
                                     cluster[4][combo[1] - 1] > MIN_CLU_ENE &&
                                     cluster[4][combo[2] - 1] > MIN_CLU_ENE &&
                                     cluster[4][combo[3] - 1] > MIN_CLU_ENE &&
                                     cluster[4][combo[4] - 1] > MIN_CLU_ENE &&
                                     cluster[4][combo[5] - 1] > MIN_CLU_ENE,
             totalEnergyLimit = totalEnergy > 350.0 && totalEnergy < 700.0,
             condTotal = totalEnergyLimit && energyLimitPerCluster;

      // Reject clusters, which do not meet energy conditions
      if (!condTotal)
        continue;
      //////////////////////////////////////////////////////////////

      // Setting up the trilateration object
      for (Int_t i = 0; i < photonNum; i++)
        _R.SetClu(i, cluster[0][combo[i] - 1],
                  cluster[1][combo[i] - 1],
                  cluster[2][combo[i] - 1],
                  cluster[3][combo[i] - 1],
                  cluster[4][combo[i] - 1]);
      //////////////////////////////////////////////////////////////

      // Preparation of the sets of 4 photon index combinations
      std::vector<std::array<Int_t, 4>> indices;
      for (Int_t i1 = 0; i1 < photonNum - 3; i1++)
        for (Int_t i2 = i1 + 1; i2 < photonNum - 2; i2++)
          for (Int_t i3 = i2 + 1; i3 < photonNum - 1; i3++)
            for (Int_t i4 = i3 + 1; i4 < photonNum; i4++)
            {
              indices.push_back({i1 + 1, i2 + 1, i3 + 1, i4 + 1});
            }
      ////////////////////////////////////////////////////////////////////////////

      Float_t totalErrorTmp = 0.;
      std::vector<std::array<Float_t, 4>>
          partialSolutionsTmp;

      std::vector<Float_t>
          partialSolutionEnergyTmp;

      // Calculation of trilateration solutions for every possible combination
      for (auto &ind : indices)
      {
        _S = _R.MySolve(ind.data());

        // Check validity of solutions and calculate errors
        Bool_t isValid[2] = {!_S.error[0], !_S.error[1]},
               anyValid = isValid[0] || isValid[1];
        Float_t errTmp[2] = {999999., 999999.};

        // Reject event if no solution is valid
        if (!anyValid)
          continue;
        ////////////////////////////////////////////////////

        for (Int_t i = 0; i < 2; i++)
        {
          if (isValid[i])
            errTmp[i] = _R.ResidualErrTot(_S.sol[i]);
        }

        if (errTmp[0] < errTmp[1])
        {
          // Save a partial solution in the vector
          totalErrorTmp += errTmp[0];
          partialSolutionsTmp.push_back({_S.sol[0][0], _S.sol[0][1], _S.sol[0][2], _S.sol[0][3]});
          Float_t energy = 0.;

          for (Int_t j = 0; j < 4; j++)
            energy += cluster[4][combo[ind[j] - 1] - 1];

          partialSolutionEnergyTmp.push_back(energy);
        }
        else if (errTmp[1] < errTmp[0])
        {
          // Save a partial solution in the vector
          totalErrorTmp += errTmp[1];
          partialSolutionsTmp.push_back({_S.sol[1][0], _S.sol[1][1], _S.sol[1][2], _S.sol[1][3]});
          Float_t energy = 0.;

          for (Int_t j = 0; j < 4; j++)
            energy += cluster[4][combo[ind[j] - 1] - 1];

          partialSolutionEnergyTmp.push_back(energy);
        }
      }
      //////////////////////////////////////////////////////
      // Choice of the solution
      if (totalErrorTmp < bestTotalError)
      {
        solutionFound = true;
        bestTotalError = totalErrorTmp;
        bestIndices = {combo[0], combo[1], combo[2], combo[3], combo[4], combo[5]};
        partialSolutions = partialSolutionsTmp;           // Update to the best partial solutions
        partialSolutionEnergy = partialSolutionEnergyTmp; // Update to the best partial solution energies
      }
    }

    if (!solutionFound)
      return ErrorHandling::ErrorCodes::NO_VALID_SIX_GAMMA_SOLUTION;

    // Fill the results of algorithm
    bestError = bestTotalError;
    WeightedMeanVertex(partialSolutions, partialSolutionEnergy, KnerecSix.fourPos);

    KnerecSix.fourMom = {0., 0., 0., 0.};

    for (Int_t i = 0; i < photonNum; i++)
    {
      Int_t ind = bestIndices[i] - 1;

      neutral_mom(cluster[0][ind],
                  cluster[1][ind],
                  cluster[2][ind],
                  cluster[4][ind],
                  KnerecSix.fourPos.data(),
                  photonSix[i].fourMom.data());

      photonSix[i].fourMomFilled = true;

      photonSix[i].clusterParams[0] = cluster[0][ind];
      photonSix[i].clusterParams[1] = cluster[1][ind];
      photonSix[i].clusterParams[2] = cluster[2][ind];
      photonSix[i].clusterParams[3] = cluster[3][ind];
      photonSix[i].clusterParams[4] = cluster[4][ind];

      photonSix[i].SetTotalVectorPhoton();

      KnerecSix.fourMom[0] += photonSix[i].fourMom[0];
      KnerecSix.fourMom[1] += photonSix[i].fourMom[1];
      KnerecSix.fourMom[2] += photonSix[i].fourMom[2];
      KnerecSix.fourMom[3] += photonSix[i].fourMom[3];
    }

    KnerecSix.SetTotalVector();
    /////////////////////////////////////////////////////////////

    return ErrorHandling::ErrorCodes::NO_ERROR;
  }
}
