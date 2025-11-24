#pragma once

#include <kloe_class.h>

namespace KLOE
{
  class NeutralReconstruction : public pm00, public ErrorHandling::ErrorLogs
  {
  public:
    NeutralReconstruction(Int_t nPhotons);
    NeutralReconstruction() : ErrorLogs("NeutralReconstruction.log"), _nPhotons(4), _nPions(2) {};
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
            std::cout << "DEBUG: Photon four momentum not filled." << std::endl;
            throw ErrorHandling::ErrorCodes::FOUR_MOM_NOT_FILLED;
          }
        }
      }
      catch (ErrorHandling::ErrorCodes &e)
      {
        std::cout << "DEBUG: Exception caught in SetPhotonParameters" << std::endl;
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
            std::cout << "DEBUG: Charged particle four momentum not filled." << std::endl;
            throw ErrorHandling::ErrorCodes::FOUR_MOM_NOT_FILLED;
          }
        }

        _Charged = pions;
        _Charged[0].SetLorentzVectors();
        _Charged[1].SetLorentzVectors();
      }
      catch (ErrorHandling::ErrorCodes &e)
      {
        std::cout << "DEBUG: Exception caught in SetChargedParameters" << std::endl;
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

    void PhotonPairingToPi0WithOmega(std::vector<neutralParticle> &photons, std::vector<chargedParticle> &pions, std::vector<Int_t> &bestPairingIndexNeutral, std::vector<Int_t> &bestPairingIndexCharged, neutralParticle &omega, std::vector<neutralParticle> &neutralPions);

    void Pi0Reconstruction();
    void Pi0Reconstruction(std::vector<neutralParticle> &pions);

    ErrorHandling::ErrorCodes ReconstructSixGammaVertex(
        const std::vector<Float_t> cluster[5],
        const std::vector<Int_t> &neu_clu_list,
        std::vector<Int_t> &bestIndices,
        Float_t &bestError,
        kaonNeutral &KnerecSix, 
        std::vector<neutralParticle> &photonFourMomSix);

    ErrorHandling::ErrorCodes ReconstructSixGammaVertexWithFourTaken(const std::vector<Float_t> cluster[5], const std::vector<Int_t> &neu_clu_list, const std::vector<Int_t> &g4taken, std::vector<Int_t> &bestIndices, Float_t &bestError, kaonNeutral &KnerecSix, std::vector<neutralParticle> &photonSix);

  private:
    Int_t _nPhotons, _nPions; // Number of photons and pions in the event
    Float_t _invMassDiffMin;
    std::vector<neutralParticle>
        _Photons, // Vector of photons for each event
        _Pions;   // Vector of pions for each event

    neutralParticle _omega; // Omega meson reconstructed from pions

    std::vector<chargedParticle> _Charged; // Vector of charged particles for each event

    // ✅ TYLKO 3 unikalne kombinacje dla 4 fotonów:
    const std::array<std::array<int, 4>, 3> _combinations = {{
        {0, 1, 2, 3}, // π⁰₁=(γ₀,γ₁), π⁰₂=(γ₂,γ₃)
        {0, 2, 1, 3}, // π⁰₁=(γ₀,γ₂), π⁰₂=(γ₁,γ₃)
        {0, 3, 1, 2}  // π⁰₁=(γ₀,γ₃), π⁰₂=(γ₁,γ₂)
    }};

    std::vector<Int_t> _bestPairingIndex,
        _bestPairingIndexOmega; // Best pairing index for photons to pions

    Reconstructor _R;
    Solution _S;
    std::vector<Int_t> _selected;

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