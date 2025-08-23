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

/**
 * @brief KLOE namespace
 */
namespace KLOE
{
    
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
         * @brief Calculation of the invariant mass using the properties of TLorentzVector.
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
        std::string getCurrentTimestamp();

        /**
         * @brief Method to get a current datestamp.
         * @returns The std::string with the current timestamp in the format: yyyy-MM-dd
         */
        std::string getCurrentDate();

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
        Int_t neu_triangle(Float_t *TrcSumFinal, Float_t *vtxSigmaFinal, Float_t Clu5Vec[4][5], Float_t *ip, Float_t *Phi4Mom, Float_t *Kne4Mom, Float_t *Kne4Vec, Float_t *trc);

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
        Int_t neu_triangle(std::vector<TLorentzVector> *Clu4Mom, std::vector<TLorentzVector> *Clu4Vec, Float_t *ip, TLorentzVector *Phi4Mom, TLorentzVector *Kne4Mom, TLorentzVector *Kne4Vec, Float_t *trc, Float_t *TrcSumFinal, Float_t *vtxSigmaFinal);
        /// @overload

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
        HypothesisCode StringToHypothesisCode(const std::string &str)
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
    };
}

// Forward declaration for FileManager
namespace KLOE
{
    class FileManager;
}

#endif