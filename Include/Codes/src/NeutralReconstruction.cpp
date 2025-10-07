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

  void NeutralReconstruction::PhotonPairingToPi0WithOmega(std::vector<neutralParticle> &photons, std::vector<chargedParticle> &pions, std::vector<Int_t> &bestPairingIndexNeutral, std::vector<Int_t> &bestPairingIndexOmega, neutralParticle &omega)
  {
    SetPhotonParameters(photons);
    SetChargedParameters(pions);

    _PhotonPairingWithOmegaCore();
    _Pi0ReconstructionCore();
    _OmegaReconstructionCore(_bestPairingIndexOmega[0]);

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
    std::vector<Int_t> indices(_nPhotons);

    std::iota(indices.begin(), indices.end(), 0);

    _invMassDiffMin = 1.e6;

    do
    {
      std::vector<Float_t>
          px(_nPions, 0.),
          py(_nPions, 0.),
          pz(_nPions, 0.),
          E(_nPions, 0.),
          invMass(_nPions, 0.);

      Float_t invMassDiffTot = 0.;

      for (Int_t i = 0; i < _nPions; i++)
      {
        px[i] = _Photons[indices[i * _nPions]].fourMom[0] + _Photons[indices[i * _nPions + 1]].fourMom[0];
        py[i] = _Photons[indices[i * _nPions]].fourMom[1] + _Photons[indices[i * _nPions + 1]].fourMom[1];
        pz[i] = _Photons[indices[i * _nPions]].fourMom[2] + _Photons[indices[i * _nPions + 1]].fourMom[2];
        E[i] = _Photons[indices[i * _nPions]].fourMom[3] + _Photons[indices[i * _nPions + 1]].fourMom[3];

        invMass[i] = sqrt(pow(E[i], 2) - pow(px[i], 2) - pow(py[i], 2) - pow(pz[i], 2));

        invMassDiffTot += pow(invMass[i] - mPi0, 2);
      }

      invMassDiffTot = sqrt(invMassDiffTot);

      if (invMassDiffTot < _invMassDiffMin)
      {
        _invMassDiffMin = invMassDiffTot;
        _bestPairingIndex = indices;
      }

    } while (std::next_permutation(indices.begin(), indices.end()));
  }

  void NeutralReconstruction::_PhotonPairingWithOmegaCore()
  {
    std::vector<Int_t> indices(_nPhotons);

    std::iota(indices.begin(), indices.end(), 0);

    Double_t mPi0Sigma = 0.0005, // MeV
        mOmegaSigma = 0.13;      // MeV

    _invMassDiffMin = 1.e6;

    do
    {
      std::vector<Float_t>
          px(_nPions, 0.),
          py(_nPions, 0.),
          pz(_nPions, 0.),
          E(_nPions, 0.),
          invMass(_nPions, 0.);

      Float_t invMassDiffTot = 0.,
              invMassDiffTotOmega[2] = {0.};

      for (Int_t i = 0; i < _nPions; i++)
      {
        px[i] = _Photons[indices[i * _nPions]].fourMom[0] + _Photons[indices[i * _nPions + 1]].fourMom[0];
        py[i] = _Photons[indices[i * _nPions]].fourMom[1] + _Photons[indices[i * _nPions + 1]].fourMom[1];
        pz[i] = _Photons[indices[i * _nPions]].fourMom[2] + _Photons[indices[i * _nPions + 1]].fourMom[2];
        E[i] = _Photons[indices[i * _nPions]].fourMom[3] + _Photons[indices[i * _nPions + 1]].fourMom[3];

        invMass[i] = sqrt(pow(E[i], 2) - pow(px[i], 2) - pow(py[i], 2) - pow(pz[i], 2));

        invMassDiffTot += pow(invMass[i] - mPi0, 2);
      }

      invMassDiffTot = invMassDiffTot / pow(mPi0Sigma, 2);

      _Pi0ReconstructionCore(indices);

      for (Int_t i = 0; i < _nPions; i++)
      {
        _OmegaReconstructionCore(i);
        invMassDiffTotOmega[i] = invMassDiffTot +
                                 pow(_omega.mass - mOmega, 2) / pow(mOmegaSigma, 2);
      }

      if (invMassDiffTotOmega[0] < invMassDiffTotOmega[1])
      {
        invMassDiffTot = invMassDiffTotOmega[0];
        _bestPairingIndexOmega = {0, 1};
      }
      else
      {
        invMassDiffTot = invMassDiffTotOmega[1];
        _bestPairingIndexOmega = {1, 0};
      }

      if (invMassDiffTot < _invMassDiffMin)
      {
        _invMassDiffMin = invMassDiffTot;
        _bestPairingIndex = indices;
      }

    } while (std::next_permutation(indices.begin(), indices.end()));
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
}
