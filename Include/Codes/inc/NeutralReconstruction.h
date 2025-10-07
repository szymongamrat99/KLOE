#pragma once

#include <kloe_class.h>

namespace KLOE
{
  class NeutralReconstruction : public pm00, public ErrorHandling::ErrorLogs
  {
  public:
    NeutralReconstruction(Int_t nPhotons);
    NeutralReconstruction() : ErrorLogs("NeutralReconstruction.log") {};
    void SetPhotonParameters(const std::vector<neutralParticle> &photons)
    {
      try
      {
        if (photons.size() != _nPhotons)
        {
          throw ErrorHandling::ErrorCodes::INVALID_PHOTON_NUMBER;
        }

        _Photons = photons;

        // Check if four mom is filled
        for (auto &photon : _Photons)
        {
          if (!photon.fourMomFilled)
          {
            throw ErrorHandling::ErrorCodes::FOUR_MOM_NOT_FILLED;
          }
        }
      }
      catch (ErrorHandling::ErrorCodes &e)
      {
        ErrorLogs::getErrLog(e);
        return;
      }
    }
    void SetChargedParameters(const std::vector<chargedParticle> &pions)
    {
      try
      {
        // Check if four mom is filled
        for (auto &pion : pions)
        {
          if (!pion.fourMomFilled)
          {
            throw ErrorHandling::ErrorCodes::FOUR_MOM_NOT_FILLED;
          }
        }

        _Charged = pions;
        _Charged[0].SetLorentzVectors();
        _Charged[1].SetLorentzVectors();
      }
      catch (ErrorHandling::ErrorCodes &e)
      {
        ErrorLogs::getErrLog(e);
        return;
      }
    }
    void GetInvMassMinimalDifference(Float_t &invMassDiffMin) const
    {
      invMassDiffMin = _invMassDiffMin;
    }
    void GetBestPairingIndex(std::vector<Int_t> &bestPairingIndex) const
    {
      bestPairingIndex = _bestPairingIndex;
    }
    void GetPions(std::vector<neutralParticle> &pions) const
    {
      pions = _Pions;
    }
    void GetPhotons(std::vector<neutralParticle> &photons) const
    {
      photons = _Photons;
    }

    void SetNumberOfPhotons(Int_t nPhotons)
    {
      try
      {
        _nPhotons = nPhotons;
        if (!_checkNPhotons())
        {
          throw std::invalid_argument("Number of photons must be a positive even integer.");
        }
        _nPions = _nPhotons / 2.;
        _Photons.resize(_nPhotons);
        _Pions.resize(_nPions);
      }
      catch (const std::invalid_argument &e)
      {
        std::cerr << "Error: " << e.what() << std::endl;
      }
    }

    void PhotonPairingToPi0();
    void PhotonPairingToPi0(std::vector<neutralParticle> &photons, std::vector<Int_t> &bestPairingIndex);

    void PhotonPairingToPi0WithOmega(std::vector<neutralParticle> &photons, std::vector<chargedParticle> &pions, std::vector<Int_t> &bestPairingIndexNeutral, std::vector<Int_t> &bestPairingIndexCharged, neutralParticle &omega);

    void Pi0Reconstruction();
    void Pi0Reconstruction(std::vector<neutralParticle> &pions);

  private:
    Int_t _nPhotons, _nPions; // Number of photons and pions in the event
    Float_t _invMassDiffMin;
    std::vector<neutralParticle>
        _Photons, // Vector of photons for each event
        _Pions;   // Vector of pions for each event

    neutralParticle _omega;          // Omega meson reconstructed from pions

    std::vector<chargedParticle> _Charged; // Vector of charged particles for each event

    std::vector<Int_t> _bestPairingIndex,
        _bestPairingIndexOmega; // Best pairing index for photons to pions

    Bool_t _checkNPhotons() const
    {
      return (_nPhotons > 0 && _nPhotons % 2 == 0);
    }

    void _PhotonPairingCore();
    void _PhotonPairingWithOmegaCore();
    void _Pi0ReconstructionCore();
    void _Pi0ReconstructionCore(std::vector<Int_t> indices);
    void _OmegaReconstructionCore(Int_t index);
  };
}