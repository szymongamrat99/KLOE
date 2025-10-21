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

    _omega.SetTotalVector();

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
        mOmegaSigma = 30.0;    // MeV

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
