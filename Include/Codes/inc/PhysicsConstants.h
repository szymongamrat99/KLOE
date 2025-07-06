// PhysicsConstants.h
#pragma once

#include <string>
#include <nlohmann/json.hpp>

/**
 * @brief Manages and provides access to experiment-specific configuration.
 * This class combines compile-time constants with runtime-loaded parameters
 * from a JSON file. All runtime parameters are loaded into strongly-typed
 * member variables during construction for easy and type-safe access.
 * This class is NOT a Singleton; it's designed for Dependency Injection.
 */
class PhysicsConstants
{
public:
  // --- Compile-Time Constants (Known at compilation, always available) ---
  // These values are embedded directly into the executable.
  // They are ideal for fundamental, unchanging physical constants.
  static constexpr double cVel = 29.9792458;       ///< Speed of light in cm/ns
  static constexpr double eleCh = 1.602176634E-19; ///< Electron charge in C
  static constexpr double hBar = 6.582119569E-34;  ///< Reduced Planck constant in MeV*s

  // --- Populate PhysicalConstants members from the loaded JSON ---

  // Particles' masses
  static constexpr double mPhi = 1019.461;     // Mass of phi in MeV/c^2
  static constexpr double mPi0 = 134.9768;     // Mass of pi0 in MeV
  static constexpr double mPiCh = 139.57039;   // Mass of pi+/- in MeV
  static constexpr double mMuon = 105.6583755; // Mass of muon in MeV
  static constexpr double mElec = 0.510998950; // Mass of electron in MeV
  static constexpr double mOmega = 782.66;     // Mass of omega in MeV

  // Branching ratios
  // Phi
  static constexpr double br_phi_kskl = 0.339;
  static constexpr double br_phi_omegapi0 = 4.7E-5; // Branching ratio for phi -> omega pi0

  // K-short
  static constexpr double br_ks_pi0pi0 = 0.3069;       // Branching ratio
  static constexpr double br_ks_pippim = 0.6920;       // Branching ratio for K-short -> pi+ pi-
  static constexpr double br_ks_pippimgamma = 1.79E-3; // Branching ratio for K-short -> pi+ pi- gamma
  static constexpr double br_ks_piele = 7.04E-4;       // Branching ratio for K-short -> pi+ e-
  static constexpr double br_ks_pimu = 4.56E-4;        // Branching ratio for K-short -> pi+ mu-

  // K-long
  static constexpr double br_kl_pi0pi0 = 8.64E-4;   // Branching ratio for K-long -> pi0 pi0
  static constexpr double br_kl_pippim = 1.967E-3;  // Branching ratio for K-long -> pi+ pi-
  static constexpr double br_kl_pippimpi0 = 0.1254; // Branching ratio for K-long -> pi+ pi- pi0
  static constexpr double br_kl_3pi0 = 0.1952;      // Branching ratio for K-long -> 3 pi0
  static constexpr double br_kl_piele = 0.4055;     // Branching ratio for K-long -> pi+ e-
  static constexpr double br_kl_pimu = 0.2704;      // Branching ratio for K-long -> pi+ mu-

  static constexpr double tau_S_CPT = 0.8954E-1;      // CPT kaon S lifetime in ns
  static constexpr double delta_mass_CPT = 0.5293E10; // CPT delta mass in hbar s^-1
  static constexpr double Im_CPT = -0.002;            // CPT imaginary part of epsilon in degrees
  static constexpr double phi_pm_CPT = 43.51;
  static constexpr double phi_00_CPT = 43.52; // CPT phi00 parameter in degrees
  static constexpr double TRF = 2.715;        // Time of DAFNE bunch in ns

  // --- Constructor ---
  /**
   * @brief Constructor for PhysicsConstants.
   * Loads runtime configuration parameters from the specified JSON file
   * and assigns them to strongly-typed member variables.
   * @param jsonFilePath The full path to the JSON file containing runtime configuration.
   * @throws std::runtime_error if the JSON file cannot be opened or parsed,
   * or if critical parameters are missing/invalid.
   */
  explicit PhysicsConstants(const std::string &jsonFilePath);

  // --- Public Getters for Runtime Configuration (loaded from JSON) ---
  // These methods provide type-safe access to values populated during construction.
  double getKaonMass() const { return mK0; }
  double getTauSNonCPT() const { return tau_S_nonCPT; } // K-short lifetime in ns
  double getTauL() const { return tau_L; }             // K-long lifetime in ns
  double getDeltaMassNonCPT() const { return delta_mass_nonCPT; } // Delta mass in hbar s^-1
  double getModulusEpsilon() const { return mod_epsilon; } // Modulus of Epsilon
  double getRe() const { return Re; } // Real part of Epsilon
  double getImNonCPT() const { return Im_nonCPT; } // Imaginary part of Epsilon (non-CPT)
  double getPhiPMNonCPT() const { return phi_pm_nonCPT; } // Phi PM (non-CPT)
  double getPhi00NonCPT() const { return phi_00_nonCPT; } // Phi 00 (non-CPT)

  // Add more getters for all runtime parameters...

private:
  // --- Private Member Variables for Runtime Configuration ---
  // These will be populated from the JSON file during construction.
  // Constants from PDG
  double mK0;
  double tau_S_nonCPT;
  double tau_L;
  double delta_mass_nonCPT;
  double mod_epsilon;
  double Re;
  double Im_nonCPT;
  double phi_pm_nonCPT;
  double phi_00_nonCPT;

  // Internal JSON object (only used during construction)
  nlohmann::json m_rawConfigData;

  /**
   * @brief Private helper to load the JSON file into m_rawConfigData.
   * @param filePath The path to the JSON configuration file.
   * @throws std::runtime_error if the file cannot be opened or parsed.
   */
  void loadJsonFile(const std::string &filePath);

  /**
   * @brief Private helper to populate member variables from the loaded JSON data.
   * @throws std::runtime_error if critical parameters are missing or have incorrect types.
   */
  void populateMembersFromJson();

  // Helper to safely get value from JSON, with error reporting
  template <typename T>
  T getValue(const std::string &jsonPath, const std::string &paramName) const;

  // Helper to get value from JSON with default, with error reporting
  template <typename T>
  T getValueOrDefault(const std::string &jsonPath, const std::string &paramName, T defaultValue) const;
};