#ifndef KLOE_CLASS_H
#define KLOE_CLASS_H

#include <fstream>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <set>
#include <regex>
#include <boost/filesystem.hpp>

#include <TMath.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TVector3.h>
#include <TVectorD.h>
#include <TFormula.h>

#include <HypothesisCodes.h>
#include <reconstructor.h>
#include <ErrorLogs.h>
#include <neutral_mom.h>
#include <MainMenu.h>
#include <const.h>

/**
 * @brief KLOE namespace
 */
namespace KLOE
{
  struct chargedParticle
  {
    chargedParticle() : fourMom(4, 0.),
                        fourPos(4, 0.),
                        trackParams(3, 0.),
                        lorentzFourMom(0., 0., 0., 0.),
                        lorentzFourPos(0., 0., 0., 0.),
                        mass(0.),
                        totalMomentum(0.) {};

    std::vector<Float_t> fourMom; /*!< 4-momentum of the charged particle */
    std::vector<Float_t> fourPos; /*!< 4-momentum of the charged particle */

    std::vector<Float_t> trackParams; /*!< Track parameters of the charged particle */

    Float_t mass;                 /*!< Mass of the charged particle */
    Float_t totalMomentum;        /*!< Total momentum of the charged particle */
    Float_t openingAngle;         /*!< Opening angle between two charged particles */
    Bool_t fourMomFilled = false; /*!< Flag to check if four momentum is filled */

    TLorentzVector lorentzFourMom; /*!< 4-momentum of the charged particle */
    TLorentzVector lorentzFourPos; /*!< 4-momentum of the charged particle */

    void SetLorentzVectors()
    {
      lorentzFourMom.SetPxPyPzE(fourMom[0], fourMom[1], fourMom[2], fourMom[3]);
      lorentzFourPos.SetXYZT(fourPos[0], fourPos[1], fourPos[2], fourPos[3]);
    };

    void FillFourMom(Float_t &px, Float_t &py, Float_t &pz, Float_t &E)
    {
      fourMom[0] = px;
      fourMom[1] = py;
      fourMom[2] = pz;
      fourMom[3] = E;

      fourMomFilled = true;
    };

    void CalculateMassFromFourMom()
    {
      mass = sqrt(pow(fourMom[3], 2) - (pow(fourMom[0], 2) + pow(fourMom[1], 2) + pow(fourMom[2], 2)));
    };

    void CalculateTotalMomentumFromFourMom()
    {
      totalMomentum = sqrt(pow(fourMom[0], 2) + pow(fourMom[1], 2) + pow(fourMom[2], 2));
    };
  };

  struct neutralParticle
  {
    neutralParticle() : fourMom(4, 0.),
                        fourPos(4, 0.),
                        total(8, 0.),
                        clusterParams(5, 0.),
                        lorentzFourMom(0., 0., 0., 0.),
                        lorentzFourPos(0., 0., 0., 0.),
                        mass(0.),
                        totalMomentum(0.) {};

    std::vector<Float_t> fourMom; /*!< 4-momentum of the neutral particle */
    std::vector<Float_t> fourPos; /*!< 4-momentum of the neutral particle */
    std::vector<Float_t> total;   /*!< Total vector of the neutral particle */

    std::vector<Float_t> clusterParams; /*!< Track parameters of the neutral particle */

    Float_t path = 0;          /*!< Path of the neutral particle */
    Float_t timeOfFlight = 0;  /*!< Time of the neutral particle */
    Float_t mass = 0;          /*!< Mass of the neutral particle */
    Float_t totalMomentum = 0; /*!< Total momentum of the neutral particle */
    Float_t openingAngle = 0;  /*!< Opening angle between two neutral particles */

    Bool_t fourMomFilled = false; /*!< Flag to check if four momentum is filled */

    TLorentzVector lorentzFourMom; /*!< 4-momentum of the neutral particle */
    TLorentzVector lorentzFourPos; /*!< 4-momentum of the neutral particle */

    void calculatePath(Float_t *neuVtx)
    {
      path = sqrt(pow(fourPos[0] - neuVtx[0], 2) + pow(fourPos[1] - neuVtx[1], 2) + pow(fourPos[2] - neuVtx[2], 2)); // cm
    };

    void calculateTimeOfFlightPhoton()
    {
      timeOfFlight = path / PhysicsConstants::cVel;
    }

    void FillFourMom(Float_t &px, Float_t &py, Float_t &pz, Float_t &E)
    {
      fourMom[0] = px;
      fourMom[1] = py;
      fourMom[2] = pz;
      fourMom[3] = E;

      fourMomFilled = true;
    };

    void SetTotalVectorPhoton()
    {
      total[0] = fourMom[0];
      total[1] = fourMom[1];
      total[2] = fourMom[2];
      total[3] = fourMom[3];
      total[4] = clusterParams[0];
      total[5] = clusterParams[1];
      total[6] = clusterParams[2];
      total[7] = clusterParams[3];
    };

    void SetTotalVector()
    {
      total.resize(9);

      total[0] = fourMom[0];
      total[1] = fourMom[1];
      total[2] = fourMom[2];
      total[3] = fourMom[3];
      total[4] = sqrt(pow(fourMom[0], 2) + pow(fourMom[1], 2) + pow(fourMom[2], 2));
      total[5] = sqrt(pow(fourMom[3], 2) - pow(total[4], 2));
      total[6] = fourPos[0];
      total[7] = fourPos[1];
      total[8] = fourPos[2];
    }

    void SetLorentzVectors()
    {
      lorentzFourMom.SetPxPyPzE(fourMom[0], fourMom[1], fourMom[2], fourMom[3]);
      lorentzFourPos.SetXYZT(fourPos[0], fourPos[1], fourPos[2], fourPos[3]);
    };

    void CalculateMassFromFourMom()
    {
      mass = sqrt(pow(fourMom[3], 2) - (pow(fourMom[0], 2) + pow(fourMom[1], 2) + pow(fourMom[2], 2)));
    };

    void CalculateTotalMomentumFromFourMom()
    {
      totalMomentum = sqrt(pow(fourMom[0], 2) + pow(fourMom[1], 2) + pow(fourMom[2], 2));
    };
  };

  struct kaonNeutral
  {
    kaonNeutral() : fourMom(4, 0.),
                    fourPos(4, 0.),
                    vtxPos(3, 0.),
                    total(10, 0.),
                    lorentzFourMom(0., 0., 0., 0.),
                    lorentzFourPos(0., 0., 0., 0.),
                    mass(0.),
                    totalMomentum(0.) {};

    std::vector<Float_t> fourMom; /*!< 4-momentum of the charged particle */
    std::vector<Float_t> fourPos; /*!< 4-momentum of the charged particle */
    std::vector<Float_t> vtxPos;  /*!< 4-momentum of the charged particle */
    std::vector<Float_t> total;   /*!< 9-vector of the charged particle */

    Float_t totalMomentum = 0; /*!< Total momentum of the kaon */
    Float_t beta = 0;          /*!< Beta of the kaon */
    Float_t gamma = 0;         /*!< Gamma of the kaon */
    Float_t mass = 0;          /*!< Mass of the kaon */
    Float_t path = 0;          /*!< Path of the kaon */
    Float_t lifetimeLAB = 0;   /*!< Lifetime of the kaon in LAB*/
    Float_t lifetimeCM = 0;    /*!< Lifetime of the kaon in CM*/

    TLorentzVector lorentzFourMom; /*!< 4-momentum of the charged particle */
    TLorentzVector lorentzFourPos; /*!< 4-momentum of the charged particle */

    void calculatePath(Float_t *ip)
    {
      path = sqrt(pow(fourPos[0] - ip[0], 2) + pow(fourPos[1] - ip[1], 2) + pow(fourPos[2] - ip[2], 2)); // cm
    };

    void calculateBeta()
    {
      calculateTotalMomentum();

      beta = totalMomentum / fourMom[3];
    };

    void calculateTotalMomentum()
    {
      totalMomentum = sqrt(pow(fourMom[0], 2) + pow(fourMom[1], 2) + pow(fourMom[2], 2));
    }

    void calculateLifetimeLAB()
    {
      calculateBeta();

      lifetimeLAB = path / (beta * PhysicsConstants::cVel);
    };

    void SetTotalVector()
    {
      total[0] = fourMom[0];
      total[1] = fourMom[1];
      total[2] = fourMom[2];
      total[3] = fourMom[3];
      total[4] = sqrt(pow(fourMom[0], 2) + pow(fourMom[1], 2) + pow(fourMom[2], 2));
      total[5] = sqrt(pow(fourMom[3], 2) - pow(total[4], 2));
      total[6] = fourPos[0];
      total[7] = fourPos[1];
      total[8] = fourPos[2];

      calculateLifetimeLAB();
      total[9] = lifetimeLAB;
      fourPos[3] = lifetimeLAB;
    };

    void SetPositionAndMomentumFromTotal()
    {
      fourMom[0] = total[0];
      fourMom[1] = total[1];
      fourMom[2] = total[2];
      fourMom[3] = total[3];
      fourPos[0] = total[6];
      fourPos[1] = total[7];
      fourPos[2] = total[8];
      fourPos[3] = total[9];
    };

    void SetLorentzVectors()
    {
      lorentzFourMom.SetPxPyPzE(fourMom[0], fourMom[1], fourMom[2], fourMom[3]);
      lorentzFourPos.SetXYZT(fourPos[0], fourPos[1], fourPos[2], fourPos[3]);
    };

    void CalculateDerivedQuantities()
    {
      totalMomentum = lorentzFourMom.P();
      beta = totalMomentum / lorentzFourMom.E();
      gamma = 1. / sqrt(1. - beta * beta);
      mass = lorentzFourMom.M();
    };
  };

  struct phiMeson
  {
    phiMeson() : fourMom(4, 0.),
                 vtxPos(3, 0.),
                 total(7, 0.) {};

    std::vector<Float_t> fourMom; /*!< 4-momentum of the charged particle */
    std::vector<Float_t> vtxPos;  /*!< 3-momentum of the charged particle */
    std::vector<Float_t> total;   /*!< 7-vector of the charged particle */

    void SetTotalVector()
    {
      total[0] = fourMom[0];
      total[1] = fourMom[1];
      total[2] = fourMom[2];
      total[3] = fourMom[3];
      total[4] = vtxPos[0];
      total[5] = vtxPos[1];
      total[6] = vtxPos[2];
    };
  };

  struct KaonProperTimes
  {
    Double_t kaon1TimeLAB = 0.;
    Double_t kaon1TimeCM = 0.;
    Double_t kaon2TimeLAB = 0.;
    Double_t kaon2TimeCM = 0.;
    Double_t deltaTimeLAB = 0.;
    Double_t deltaTimeCM = 0.;
  };
  /**
   * @class General pm00 class for KLOE analysis. Includes most fundamental functions like timestamp, datestamp, array clearing functions, etc.
   * @author @szymongamrat99
   */
  class pm00
  {
  private:
    std::string
        _MClistPath1,  /*!< First path for MC list */
        _MClistPath2,  /*!< Second path for MC list */
        _MClistPath3,  /*!< Third path for MC list */
        _DatalistPath; /*!< Path for data list */

    std::chrono::high_resolution_clock::time_point
        _start_time, /*!< Start time of the counter. */
        _end_time;   /*!< End time of the counter. */

  public:
    TLorentzVector
        phi_mom,       /*!< 4-momentum of \f$ \phi \f$ meson*/
        phi_pos,       /*!< 4-position of \f$ \phi \f$ meson*/
        kaon_mom[2],   /*!< 4-momenta of both kaons*/
        kaon_pos[2],   /*!< 4-positions of both kaons*/
        pi_ch_mom[2],  /*!< 4-momenta of both charged \f$ \pi \f$ mesons*/
        pi_ch_pos[2],  /*!< 4-positions of both charged \f$ \pi \f$ mesons*/
        pi_ne_mom[2],  /*!< 4-momenta of both \f$ \pi^{0} \f$ mesons*/
        pi_ne_pos[2],  /*!< 4-positions of both \f$ \pi^{0} \f$ mesons*/
        photon_mom[4], /*!< 4-momenta of four photons*/
        photon_pos[4]; /*!< 4-positions of four photons*/

    UInt_t
        bin_number; /*!< number of histogram bins*/

    Double_t
        x_min,      /*!< minimal value of histogram's x range*/
        x_max,      /*!< maximal value of histogram's x range*/
        y_min,      /*!< minimal value of histogram's y range*/
        y_max,      /*!< maximal value of histogram's y range*/
        inv_mass,   /*!< invariant mass of a particle*/
        angle_vec,  /*!< angle between two vectors*/
        transv,     /*!< polar angle of a vector (with respect to z-axis)*/
        azim_angle; /*!< azimuthal angle of a vector (in x-y plane)*/

    TVector3
        boost; /*!< boost vector between the frames*/

    std::vector<chargedParticle>
        pionCh;

    std::vector<neutralParticle>
        photon,
        pionNe;

    neutralParticle
        omega;

    kaonNeutral
        Kchrec,    /*!< Charged kaon reconstructed from pions*/
        Kchboost,  /*!< Charged kaon fixed with lorentz boost*/
        Knerec,    /*!< Neutral kaon reconstructed from Photons*/
        Knereclor, /*!< Neutral kaon reconstructed with Kchboost*/
        KnerecCMPhi;

    phiMeson
        phi; /*!< Phi meson */

    std::vector<Float_t>
        bhabha_mom,
        trackParameters,
        ip;

    std::vector<std::vector<Float_t>>
        cluster;

    std::vector<TH1 *> frac, frac_data;
    TH1 *data, *mc_sum, *data_sub, *mc_sub;

    /**
     * @brief Constructor of pm00 class.
     * @param mom_list Pointer to the list of TLorentzVectors with momenta
     * @param pos_list Pointer to the list of TLorentzVectors with positions
     */
    pm00(TLorentzVector *mom_list, TLorentzVector *pos_list);

    /**
     * @brief Default constructor of pm00 class.
     */
    pm00();

    /**
     * @brief Calculation of the invariant mass using the Utils::properties of TLorentzVector.
     * @param four_mom Four momentum of the particle.
     */
    void inv_mass_calc(TLorentzVector four_mom);
    void angle(TLorentzVector vec1, TLorentzVector vec2);
    void cyl_comp(TLorentzVector vec);

    void boost_vector(TLorentzVector four_mom);
    void lorentz_transf(TLorentzVector four_mom);
    void lorentz_transf(Float_t *, Float_t *, Float_t *);
    void lorentz_transf(Double_t *, Double_t *, Double_t *);

    void lorentz_transf(TVector3 &, TLorentzVector &, TLorentzVector &) const;

    void trilaterationReconstruction(TVectorD X, Double_t neuVtx[2][4], Bool_t neuVtxErr[2]);

    /**
     * @brief Sign of a value.
     * @param value Floating number, for which the sign will be found.
     * @returns Sign of the number:
     * - +1: number is positive
     * - -1: number is negative
     */
    Int_t signum(Float_t value);

    /**
     * @brief Difference of Kaons' lifetimes in their CM frames. The function uses both momenta and positions to get them in kaon's CM frames from scratch.
     * @param momKch Pointer to the TLorentzVector of \f$ K\to\pi^{+}\pi^{-} \f$ 4-momentum
     * @param posKch Pointer to the TLorentzVector of \f$ K\to\pi^{+}\pi^{-} \f$ 4-position
     * @param momKne Pointer to the TLorentzVector of \f$ K\to\pi^{0}\pi^{0} \f$ 4-momentum
     * @param posKne Pointer to the TLorentzVector of \f$ K\to\pi^{0}\pi^{0} \f$ 4-position
     * @returns \f$ \Delta t = t_{K\to\pi^{+}\pi^{-}} - t_{K\to\pi^{0}\pi^{0}} [\tau_S]\f$, which is a main observable of the signal channel.
     */
    Double_t DeltaT(TLorentzVector *momKch, TLorentzVector *posKch, TLorentzVector *momKne, TLorentzVector *posKne);

    /**
     * @brief Function to clear the 1D static / dynamic arrays.
     * @param M length of the array
     * @param array pointer to the table of elements
     */
    void Clear1DArray(UInt_t M, Float_t *array);
    /**
     * @brief Function to clear the 1D static / dynamic arrays.
     * @param M length of the array
     * @param array pointer to the table of elements
     */
    void Clear1DArray(UInt_t M, Int_t *array);
    /// @overload

    /**
     * @brief Function to clear the 2D static / dynamic arrays.
     * @param M number of rows
     * @param N number of columns
     * @param array pointer to the matrix elements
     */
    void Clear2DArray(UInt_t M, UInt_t N, Float_t **array);
    /**
     * @brief Function to clear the 2D static / dynamic arrays.
     * @param M number of rows
     * @param N number of columns
     * @param array pointer to the matrix elements
     */
    void Clear2DArray(UInt_t M, UInt_t N, Int_t **array);

    /// @overload

    /**
     * @brief Method to start the counter.
     */
    void startTimer();
    /**
     * @brief Method to stop the counter.
     * @returns The std::string with time delta between end and start of the counter.
     */
    std::string endTimer();

    /**
     * @brief Method to get a current timestamp.
     * @returns The std::string with the current timestamp in the format: yyyy-MM-dd_HHmm
     */
    std::string getCurrentTimestamp() const;

    /**
     * @brief Method to get a current datestamp.
     * @returns The std::string with the current timestamp in the format: yyyy-MM-dd
     */
    std::string getCurrentDate() const;

    /**
     * @brief Method to count the repeating elements in one table.
     * @param arr std::vector<Int_t> where the elements to check are stored
     * @returns The map<Int_t, Int_t> how many times each element is repeated.
     */
    std::map<Int_t, Int_t> CountRepeatingElements(std::vector<Int_t> &arr);

    /**
     * @brief Method check the equality of two tables up to the given size. Equality means, that all the elements are the same for both tables, order does not matter.
     * @param table1, table2 the tables to check the equality.
     * @param size number of elements in both tables.
     * @returns the number, how many elements are the same.
     */
    Int_t ArrayEquality(const Int_t *table1, const Int_t *table2, const Int_t &size);

    /**
     * @brief A function to get the list of filenames described in the subfiles. Needs to be used only in special cases.
     * @param mode Mode of operation:
     * - MC1: to get the list from all_phys
     * - MC2: to get the list from all_phys2
     * - MC3: to get the list from all_phys3
     * - Data: to get the list from data
     * @param maxNumOfFiles number of files to fetch into the analysis
     * @param fileNames address to the vector of filenames to be used in the analysis
     * @returns The function returns:
     * - 0: no errors
     * - ErrorHandling::ErrorCodes::FILE_NOT_EXIST
     */
    Int_t ListOfFiles(const std::string &mode, Int_t &maxNumOfFiles, std::vector<std::string> &fileNames);

    /**
     * @brief Reconstruction of neutral vertex using the triangle method
     * @param TrcSumFinal pointer to variable to get the TrcSum after the method application $$\sum^{4}_{i=1} (T_{cl,i} - \frac{d_{\gamma,i}}{c} - \frac{d_{K}}{c\beta_{K}})$$
     * @param vtxSigmaFinal pointer variable to get the average error of the neutral vertex reconstruction (over clusters' energies) $$\frac{\sum^{4}_{i=1}E_{cl,i}\times\sqrt{(X_{neu,i} - X_{neu,avg})^2 + (Y_{neu,i} - Y_{neu,avg})^2 + (Z_{neu,i} - Z_{neu,avg})^2}}{\sum^{4}_{i=1}E_{cl,i}}$$
     * @param Clu5Vec matrix of clusters' variables:
     * - first index: choice of cluster
     * - second index:
     *  + 0-3: ($X_{cl} [cm], Y_{cl} [cm], Z_{cl} [cm], T_{cl} [ns]$)
     *  + 4: $E_{cl} [MeV]$
     * @param ip pointer to the array of interaction point's spatial coordinates: ($X_{IP}, Y_{IP}, Z_{IP}$)
     * @param Phi4Mom pointer to the array of $\phi$ meson's 4-momentum: ($p_{x}, p_{y}, p_{z}, \sqrt{s}$)
     * @param Kne4Mom pointer to the array of kaon's 4-momentum: ($p_{x}, p_{y}, p_{z}, E_{K}$)
     * @param Kne4Vec pointer to the array of kaon's 4-position: ($x_{neu,avg}, y_{neu,avg}, z_{neu,avg}, t_{K}$)
     * @param trc pointer to the array of $$(T_{cl,i} - \frac{d_{\gamma,i}}{c} - \frac{d_{K}}{c\beta_{K}})$$ for each cluster
     * @returns The function returns:
     * - 0: no errors
     * - ErrorHandling::ErrorCodes::DELTA_LT_ZERO
     * - ErrorHandling::ErrorCodes::DENOM_EQ_ZERO
     */
    Int_t neu_triangle(Float_t *TrcSumFinal, Float_t *vtxSigmaFinal, Float_t Clu5Vec[4][5], Float_t *ip, Float_t *Phi4Mom, Float_t *Kne4Mom, Float_t *Kne4Vec, Float_t *trc) const;

    /**
     * @brief Reconstruction of neutral vertex using the triangle method
     * @param TrcSumFinal pointer to variable to get the TrcSum after the method application $$\sum^{4}_{i=1} (T_{cl,i} - \frac{d_{\gamma,i}}{c} - \frac{d_{K}}{c\beta_{K}})$$
     * @param vtxSigmaFinal pointer variable to get the average error of the neutral vertex reconstruction (over clusters' energies) $$\frac{\sum^{4}_{i=1}E_{cl,i}\times\sqrt{(X_{neu,i} - X_{neu,avg})^2 + (Y_{neu,i} - Y_{neu,avg})^2 + (Z_{neu,i} - Z_{neu,avg})^2}}{\sum^{4}_{i=1}E_{cl,i}}$$
     * @param Clu4Mom std::vector of TLorentzVector objects for 4-momentum:
     * - first index: choice of cluster
     * - vectors' content: ($p_{x,cl} [MeV/c], p_{y,cl} [MeV/c], p_{z,cl} [Mev/c], E_{cl} [MeV]$)
     * @param Clu4Vec std::vector of TLorentzVector objects for 4-position:
     * - first index: choice of cluster
     * - vectors' content: ($X_{cl} [cm], Y_{cl} [cm], Z_{cl} [cm], T_{cl} [ns]$)
     * @param ip pointer to the array of interaction point's spatial coordinates: ($X_{IP}, Y_{IP}, Z_{IP}$)
     * @param Phi4Mom pointer to the array of $\phi$ meson's 4-momentum: ($p_{x}, p_{y}, p_{z}, \sqrt{s}$)
     * @param Kne4Mom pointer to the array of kaon's 4-momentum: ($p_{x}, p_{y}, p_{z}, E_{K}$)
     * @param Kne4Vec pointer to the array of kaon's 4-position: ($x_{neu,avg}, y_{neu,avg}, z_{neu,avg}, t_{K}$)
     * @param trc pointer to the array of $$(T_{cl,i} - \frac{d_{\gamma,i}}{c} - \frac{d_{K}}{c\beta_{K}})$$ for each cluster
     * @returns The function returns:
     * - 0: no errors
     * - ErrorHandling::ErrorCodes::DELTA_LT_ZERO
     * - ErrorHandling::ErrorCodes::DENOM_EQ_ZERO
     */
    Int_t neu_triangle(std::vector<TLorentzVector> *Clu4Mom, std::vector<TLorentzVector> *Clu4Vec, Float_t *ip, TLorentzVector *Phi4Mom, TLorentzVector *Kne4Mom, TLorentzVector *Kne4Vec, Float_t *trc, Float_t *TrcSumFinal, Float_t *vtxSigmaFinal) const;
    /// @overload

    ErrorHandling::ErrorCodes NeuCluWrongCheck(std::vector<Int_t> neuclulist, phiMeson phi, kaonNeutral Kchboost, std::vector<Float_t> Xcl, std::vector<Float_t> Ycl, std::vector<Float_t> Zcl, std::vector<Float_t> Tcl, std::vector<Float_t> Enecl);

    static ErrorHandling::ErrorCodes triangleReconstruction(std::vector<Int_t> g4taken_kinfit, std::vector<Float_t> cluster[5], std::vector<Int_t> Asscl, std::vector<Float_t> bhabha_mom, std::vector<Float_t> Kchboost, std::vector<Float_t> ip, std::vector<Float_t> &Knetriangle, std::vector<Float_t> gammatriangle[4], Float_t &minv4gam, std::vector<Float_t> &trcfinal, ErrorHandling::ErrorLogs &logger);

    ErrorHandling::ErrorCodes triangleReconstruction(std::vector<neutralParticle> &photon, phiMeson phi, kaonNeutral Kchboost, Float_t *ip, kaonNeutral &Knetriangle);

    /**
     * @brief Method to initialize a TChain with the list of files. Defined per branch, which files (Prod2ntu, Prod2root, Old analysis) are to be taken into account.
     * @param chain_init address to the externally initialized TChain object
     * @param dataTypeOpt address to the variable, which stores Controls::DataType value
     * @param firstData first data file index
     * @param lastData last data file index
     * @param firstMC first MC file index
     * @param lastMC last MC file index
     * @param logger address to the ErrorHandling::ErrorLogs object to log the errors / infos during analysis operation
     */
    // void chainInit(TChain &chain_init, Controls::DataType &dataTypeOpt, UInt_t &firstData, UInt_t &lastData, UInt_t &firstMC, UInt_t &lastMC, ErrorHandling::ErrorLogs &logger, Int_t csFlag);

    void checkFilesList();

    /**
     * @brief Method to set the data type, which will be used during the analysis
     * @param dataType address to the variable, which stores Controls::DataType value
     * @param dataFlag address to variable where result will be stored
     * @param mcflag the flag, if an event comes from MC simulation:
     * - 0: data event
     * - 1: MC event
     * @param mctruth the indicator of a decay channel:
     * - 0: \f$ K_{S}K_{L}\to \pi^{+}\pi^{-}\pi^{0}\pi^{0} \f$ with error
     * - 1: \f$ K_{S}K_{L}\to \pi^{+}\pi^{-}\pi^{0}\pi^{0} \f$ after cuts
     * - 2: \f$ K_{S}K_{L}\to \pi^{+}\pi^{-}\pi^{0}\pi^{0} \f$ rejected by cuts
     * - 3: regeneration
     * - 4: \f$ \omega\pi^{0}\to\pi^{+}\pi^{-}\pi^{0}\pi^{0} \f$
     * - 5: \f$ K_{S}K_{L}\to \pi^{+}\pi^{-}3\pi^{0} \f$
     * - 6: \f$ K_{S}K_{L}\to \pi^{0} \pi^{\pm}l^{\mp}\nu \f$
     * - 7: Other background
     * - 8: \f$ K_{S}K_{L}\to \pi^{+}\pi^{-}\pi^{+}\pi^{-} \f$
     */
    void dataFlagSetter(Controls::DataType &dataType, bool &dataFlag, int &mcflag, int &mctruth);

    void ConditionInitializer(std::vector<TFormula> &formula_vector, std::vector<TString> &formula_names) const;

    void SetSingleConditionParameters(TFormula &formula, Double_t *mean_sigma) const;

    void SetAllConditionParameters(std::vector<TFormula> &formula_vector, std::vector<Double_t *> &mean_sigma_vector) const;

    Bool_t GetSingleConditionValue(TFormula &formula, Double_t &x) const;

    void GetAllConditionValue(std::vector<TFormula> &formula_vector, std::vector<Double_t> &x, std::vector<Bool_t> &values) const;

    static Double_t MomMinAngleFunction(Double_t *x, Double_t *p)
    {
      Double_t
          numerator = 0,
          denominator = 0,
          value = 0;

      numerator = sqrt(pow(p[0] - p[4] * cos(p[5]), 2) + pow(p[1] - p[4] * sin(p[5]), 2)) +
                  sqrt(pow(p[2] - p[4] * cos(M_PI + p[5]), 2) + pow(p[3] - p[4] * sin(M_PI + p[5]), 2));
      denominator = sqrt(pow(p[0], 2) + pow(p[1], 2)) + sqrt(pow(p[2], 2) + pow(p[3], 2));

      if (denominator != 0 && !TMath::IsNaN(denominator) && !TMath::IsNaN(numerator))
        value = numerator / denominator;
      else
        value = 999.;

      return value;
    }

    std::string CreateDatedFolder(const std::string &basePath);

    Double_t WeightedAverageSymmetric(std::vector<Double_t> &values, std::vector<Double_t> &errors);
    Double_t WAvgSymmError(std::vector<Double_t> &errors);

    Double_t WeightedAverageAsymmetric(std::vector<Double_t> &values, std::vector<Double_t> &errorsLow, std::vector<Double_t> &errorsUp);
    Double_t WAvgAsymmError(std::vector<Double_t> &errorsLow, std::vector<Double_t> &errorsUp);

    Double_t TwoBodyDecayMass(Double_t M, Double_t m1, Double_t m2);

    /**
     * @brief Extract unique run numbers from filenames in a directory.
     * @param directory Path to the directory with files.
     * @param pattern Regex pattern to extract run numbers (e.g. R"(data_stream42_(\\d+)\\.ntu)").
     * @returns Set of unique run numbers found in filenames.
     */
    std::set<Int_t> ExtractUniqueRuns(const std::string &directory, const std::string &pattern);

    // -------------------------------------------------------------------------------
    static HypothesisCode StringToHypothesisCode(const std::string &str)
    {
      static const std::map<std::string, HypothesisCode> enumMap = {
          {"SIGNAL", HypothesisCode::SIGNAL},
          {"OMEGAPI", HypothesisCode::OMEGAPI},
          {"FOUR_PI", HypothesisCode::FOUR_PI}};

      auto it = enumMap.find(str);
      if (it != enumMap.end())
      {
        return it->second;
      }

      std::cout << "Wrong analysis code!" << std::endl;

      return HypothesisCode::INVALID_VALUE;
    };

    static std::string HypothesisCodeToString(HypothesisCode code)
    {
      static const std::map<HypothesisCode, std::string> enumMap = {
          {HypothesisCode::SIGNAL, "SIGNAL"},
          {HypothesisCode::OMEGAPI, "OMEGAPI"},
          {HypothesisCode::FOUR_PI, "FOUR_PI"}};

      auto it = enumMap.find(code);
      if (it != enumMap.end())
      {
        return it->second;
      }

      std::cout << "Wrong analysis code!" << std::endl;

      return "INVALID_VALUE";
    };

    // -------------------------------------------------------------------------------
    static TrilaterationCode StringToTrilaterationCode(const std::string &str)
    {
      static const std::map<std::string, TrilaterationCode> enumMap = {
          {"TWO_PI0", TrilaterationCode::TWO_PI0},
          {"THREE_PI0", TrilaterationCode::THREE_PI0}};

      auto it = enumMap.find(str);
      if (it != enumMap.end())
      {
        return it->second;
      }

      std::cout << "Wrong analysis code!" << std::endl;

      return TrilaterationCode::INVALID_VALUE;
    };

    void CorrectClusterTime(Float_t T0, std::vector<Float_t> &cluTimeCorrected);

    void PhotonPairingToPi0(std::vector<Float_t> *photon, Int_t numOfPhotons, std::vector<Int_t> &bestPairingIndex);

    /**
     * @brief Oblicza czasy własne kaonów i ich różnice
     * @param kaon1Mom Wektor pędu pierwszego kaona [px, py, pz, E]
     * @param kaon1Pos Wektor położenia pierwszego kaona [x, y, z]
     * @param kaon2Mom Wektor pędu drugiego kaona [px, py, pz, E]
     * @param kaon2Pos Wektor położenia drugiego kaona [x, y, z]
     * @param phiMom Wektor pędu układu phi [px, py, pz, E]
     * @param ipPos Wektor położenia IP [x, y, z]
     * @return Struktura z czasami własnymi w układzie LAB i CM oraz różnicami
     */
    KaonProperTimes CalculateKaonProperTimes(const std::vector<Float_t> &kaon1Mom,
                                             const std::vector<Float_t> &kaon1Pos,
                                             const std::vector<Float_t> &kaon2Mom,
                                             const std::vector<Float_t> &kaon2Pos,
                                             const std::vector<Float_t> &ipPos);

    /**
     * @brief Oblicza czas własny pojedynczego kaona
     * @param kaonMom Wektor pędu kaona [px, py, pz, E]
     * @param kaonPos Wektor położenia kaona [x, y, z]
     * @param phiMom Wektor pędu układu phi [px, py, pz, E]
     * @param ipPos Wektor położenia IP [x, y, z]
     * @param timeLAB Referencja na czas w układzie LAB (wyjście)
     * @param timeCM Referencja na czas w układzie CM (wyjście)
     */
    void CalculateSingleKaonTime(const std::vector<Float_t> &kaonMom,
                                 const std::vector<Float_t> &kaonPos,
                                 const std::vector<Float_t> &ipPos,
                                 Double_t &timeLAB,
                                 Double_t &timeCM);

    // File Management Methods
    /**
     * @brief Tworzy linki symboliczne w folderze current dla najnowszych plików
     * @param rootFilesDir Ścieżka do katalogu root_files
     * @return true jeśli operacja się udała, false w przeciwnym razie
     */
    bool CreateCurrentLinks(const std::string &rootFilesDir = "root_files");

    /**
     * @brief Sprawdza czy string ma format daty YYYY-MM-DD
     * @param folderName Nazwa folderu do sprawdzenia
     * @return true jeśli nazwa ma format daty, false w przeciwnym razie
     */
    bool IsValidDateFormat(const std::string &folderName);

    // Weighted mean
    void WeightedMeanVertex(const std::vector<std::array<Float_t, 4>> &values, const std::vector<Float_t> &weights, std::vector<Float_t> &weightedMean)
    {
      if (values.size() != weights.size() || values.empty())
        return;

      weightedMean.resize(4, 0.0);

      for (size_t i = 0; i < 4; ++i)
      {
        Double_t numerator = 0.0;
        Double_t denominator = 0.0;

        for (size_t j = 0; j < values.size(); ++j)
        {
          numerator += values[j][i] * weights[j];
          denominator += weights[j];
        }

        if (denominator != 0)
          weightedMean[i] = numerator / denominator;
        else
          weightedMean[i] = 0.0;
      }
    };
  };
}

// Forward declaration for FileManager
namespace KLOE
{
  class FileManager;
}

#endif