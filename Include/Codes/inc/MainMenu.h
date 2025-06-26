#ifndef MAIN_MENU_H
#define MAIN_MENU_H

#include <iostream>
#include <map>

#include <TString.h>

/**
 * @namespace Controls
 * @brief A namespace for Controls in the project
 */
namespace Controls
{
  /**
   * @enum Controls::MainMenu
   * The options for Main Menu of the project
   */

  /**
   * @var Controls::MainMenu::GEN_VARS
   * Part of the analysis to get:
   * - generated variables from Monte Carlo
   * - mctruth splitted files
   * - expected Monte Carlo generated cluster choice
   * @var Controls::MainMenu::KCHREC_NO_BOOST
   * Reconstruction of charged decay pm without boost method
   * @var Controls::MainMenu::KCHREC_BOOST
   * Use of boost method to improve results of pm reconstruction
   * @var Controls::MainMenu::IP_EV_BY_EV
   * Reconstruction of IP event by event using:
   * - plane intersection method
   * - point of closes approach method
   * @var Controls::MainMenu::OMEGA_REC
   * Reconstruction of \f$(\omega \pi^{0})\f$ decay for the rejection improvement
   * @var Controls::MainMenu::REGEN_REJ
   * Methods to reject regeneration channel from the sample
   * @var Controls::MainMenu::KNEREC_TRILAT
   * Reconstruction of 00 using trilateration method
   * @var Controls::MainMenu::KNEREC_TRIANGLE
   * Reconstruction of 00 using triangle method
   * @var Controls::MainMenu::EFF_SIG_TO_BCG
   * Efficiency of signal with respect to background, signal to background ratio, etc.
   * @var Controls::MainMenu::KIN_FITS
   * Global kinematic fit under signal hypothesis.
   * @var Controls::MainMenu::TRANSF_TO_CM
   * Transformation to CM of all variables in the analysis.
   * @var Controls::MainMenu::CPV_NORM
   * Final CP normalization to get the CP violation parameters.
   * @var Controls::MainMenu::PLOTS
   * Auxiliary plots of the analysis.
   * @var Controls::MainMenu::EXIT
   * Exit from the analysis chain.
   * @var Controls::MainMenu::OPT_TOT
   * Total number of options (for loop limits).
   */
  enum class MainMenu
  {
    GEN_VARS = 1,
    KCHREC_NO_BOOST = 2,
    KCHREC_BOOST = 3,
    IP_EV_BY_EV = 4,
    OMEGA_REC = 5,
    REGEN_REJ = 6,
    KNEREC_TRILAT = 7,
    KNEREC_TRIANGLE = 8,
    EFF_SIG_TO_BCG = 9,
    KIN_FITS = 10,
    TRANSF_TO_CM = 11,
    CPV_NORM = 12,
    PLOTS = 13,
    COVMATRIX = 14,
    EXIT = 15,

    OPT_TOT = 16
  };

  enum class MainMenuControlSample
  {
    KCHREC = 1,
    COV_MATRIX = 2,
    CORR_FACTOR = 3,
    EXIT = 4,

    OPT_TOT = 5
  };

  /**
   * @enum Controls::KchRecMenu
   * The options for \f$ K\to \pi^{+} \pi^{-} \f$ reconstruction.
   */

  /**
   * @var Controls::KchRecMenu::K_MASS
   * The vertex and tracks fulfilling the best kaon mass constraint are taken.
   * @var Controls::KchRecMenu::KSKL
   * The vertex closest to the IP and fulfilling the best kaon mass is taken as \f$ K_S \f$, then another one is taken as \f$ K_L \f$.
   * @var Controls::KchRecMenu::CLOSEST
   * The vertex closest to the IP is taken (no kaon mass condition).
   * @var Controls::KchRecMenu::BOOST
   * Reconstruction using Lorentz boost method.
   * @var Controls::KchRecMenu::EXIT
   * Exit from the analysis chain.
   * @var Controls::KchRecMenu::OPT_TOT
   * Total number of options (for loop limits).
   */
  enum class KchRecMenu
  {
    K_MASS = 1,
    KSKL = 2,
    CLOSEST = 3,
    BOOST = 4,
    EXIT = 5,

    OPT_TOT = 6
  };

    /**
   * @enum Controls::NeutRecMenu
   * The options for \f$ K\to \pi^{0} \pi^{0} \f$ reconstruction.
   */

  /**
   * @var Controls::NeutRecMenu::BARE_TRILATERATION
   * The trilateration used without any corrections.
   * @var Controls::NeutRecMenu::TRILATERATION_KIN_FIT
   * The trilateration used along with the kinematic fit to improve the efficiency of proper vertices choice.
   * @var Controls::NeutRecMenu::TRIANGLE
   * The reconstruction of neutral vertex using triangle method, but on clusters chosen with trilateration.
   * @var Controls::NeutRecMenu::COMP_OF_MET
   * Comparison of methods using plots, statistics, etc.
   * @var Controls::NeutRecMenu::EXIT
   * Exit from the analysis chain.
   * @var Controls::NeutRecMenu::OPT_TOT
   * Total number of options (for loop limits).
   */
  enum class NeutRecMenu
  {
    BARE_TRILATERATION = 1,
    TRILATERATION_KIN_FIT = 2,
    TRIANGLE = 3,
    COMP_OF_MET = 4,
    EXIT = 5,

    OPT_TOT = 6
  };

  enum class CPFitMenu
  {
    HALF_SIGNAL_MC = 1,
    HALF_SIG_BCG_MC = 2,
    MC_DATA = 3,
    EXIT = 4,

    OPT_TOT = 5
  };

  enum class GenVars
  {
    GEN_VARS = 1,
    SPLIT_CHANN = 2,
    EXIT = 3,

    OPT_TOT = 4
  };

  enum class OmegaRec
  {
    OMEGA_REC = 1,
    OMEGA_CUTS = 2,
    PLOTS = 3,
    EXIT = 4,

    OPT_TOT = 5
  };

  enum class Regen
  {
    REGEN_REJEC_TEST = 1,
    PLOTS = 2,
    EXIT = 3,

    OPT_TOT = 4
  };

  enum class Efficiency
  {
    EFF_SCAN = 1,
    EFF_DIST = 2,
    EXIT = 3,

    OPT_TOT = 4
  };

  enum class Plots
  {
    PLOTS_GENERAL = 1,
    PLOTS_NORM = 2,
    EXIT = 3,

    OPT_TOT = 4
  };

  enum class CovMatrix
  {
    USING_MC_DATA = 1,
    USING_CONTROL_SAMPLE = 2,
    EXIT = 3,

    OPT_TOT = 4
  };

  enum class DataType
  {
    SIGNAL_TOT = 1,
    SIG_BCG = 2,
    MC_DATA = 3,
    SIGNAL_MAX = 4,

    OPT_TOT = 5
  };

  template <typename T>
  inline std::istream &operator>>(std::istream &is, T &opt)
  {
    int a;
    is >> a;
    opt = static_cast<T>(a); // Rzutowanie z int na typ docelowy

    return is;
  }

  class Menu
  {
  private:
    std::map<int, TString> MenuOpt;
    std::vector<TString> fileOpt;

    const int ChooseMenu;
    std::vector<TString> MenuName = {"KchRec Menu", "NeutRec Menu", "Data Type", "Analysis file", "Final CPV Fit", "Generated Variables", "Omega-pi0 Reconstruction Menu", "Efficiency Check Menu", "Regeneration Rejection Menu", "Plots Menu", "Main Menu", "Covariant Matrix Determination", "Main Menu"};

    TString ChooseOpt = "Choose the option: ";

  public:
    Menu(int ChooseMenu) : ChooseMenu(ChooseMenu)
    {
      switch (ChooseMenu)
      {
      case 0:
      {
        MenuOpt[int(KchRecMenu::K_MASS)] = Form("%d. Using Kaon mass as a benchmark.", int(KchRecMenu::K_MASS));
        MenuOpt[int(KchRecMenu::KSKL)] = Form("%d. Finding KS and then KL.", int(KchRecMenu::KSKL));
        MenuOpt[int(KchRecMenu::CLOSEST)] = Form("%d. Charged decay closest to the Bhabha IP.", int(KchRecMenu::CLOSEST));
        MenuOpt[int(KchRecMenu::BOOST)] = Form("%d. Correction of reconstructed momentum with boost method.", int(KchRecMenu::BOOST));
        MenuOpt[int(KchRecMenu::EXIT)] = Form("%d. Exit.", int(KchRecMenu::EXIT));

        break;
      }
      case 1:
      {
        MenuOpt[int(NeutRecMenu::BARE_TRILATERATION)] = Form("%d. Bare trilateration.", int(NeutRecMenu::BARE_TRILATERATION));
        MenuOpt[int(NeutRecMenu::TRILATERATION_KIN_FIT)] = Form("%d. Trilateration with kinematic fit.", int(NeutRecMenu::TRILATERATION_KIN_FIT));
        MenuOpt[int(NeutRecMenu::TRIANGLE)] = Form("%d. Trilateration with triangle.", int(NeutRecMenu::TRIANGLE));
        MenuOpt[int(NeutRecMenu::COMP_OF_MET)] = Form("%d. Comparison of methods, plotting, etc.", int(NeutRecMenu::COMP_OF_MET));
        MenuOpt[int(NeutRecMenu::EXIT)] = Form("%d. Exit.", int(NeutRecMenu::EXIT));

        break;
      }
      case 2:
      {
        MenuOpt[int(DataType::SIGNAL_TOT)] = Form("%d. Total MC signal.", int(DataType::SIGNAL_TOT));
        MenuOpt[int(DataType::SIG_BCG)] = Form("%d. Signal + background MC.", int(DataType::SIG_BCG));
        MenuOpt[int(DataType::MC_DATA)] = Form("%d. MC + Data.", int(DataType::MC_DATA));
        MenuOpt[int(DataType::SIGNAL_MAX)] = Form("%d. MC - without cuts and errors.", int(DataType::SIGNAL_MAX));

        break;
      }
      case 4:
      {
        MenuOpt[int(CPFitMenu::HALF_SIGNAL_MC)] = Form("%d. Only Signal; 1/2 MC vs. 1/2 'Data'", int(CPFitMenu::HALF_SIGNAL_MC));
        MenuOpt[int(CPFitMenu::HALF_SIG_BCG_MC)] = Form("%d. Signal + Bcg; 1/2 MC vs. 1/2 'Data'", int(CPFitMenu::HALF_SIG_BCG_MC));
        MenuOpt[int(CPFitMenu::MC_DATA)] = Form("%d. Signal + Bcg; MC vs. Data", int(CPFitMenu::MC_DATA));
        MenuOpt[int(CPFitMenu::EXIT)] = Form("%d. Exit.", int(CPFitMenu::EXIT));

        break;
      }
      case 5:
      {
        MenuOpt[int(GenVars::GEN_VARS)] = Form("%d. Generated variables to channel mapper.", int(GenVars::GEN_VARS));
        MenuOpt[int(GenVars::SPLIT_CHANN)] = Form("%d. Split stream into decay channels.", int(GenVars::SPLIT_CHANN));
        MenuOpt[int(GenVars::EXIT)] = Form("%d. Exit.", int(GenVars::EXIT));

        break;
      }
      case 6:
      {
        MenuOpt[int(OmegaRec::OMEGA_REC)] = Form("%d. Reconstruction of Omega channel based on the kinematic fit.", int(OmegaRec::OMEGA_REC));
        MenuOpt[int(OmegaRec::OMEGA_CUTS)] = Form("%d. Adjustment of cuts to reject Omega channel.", int(OmegaRec::OMEGA_CUTS));
        MenuOpt[int(OmegaRec::PLOTS)] = Form("%d. Plots, comparisons, etc.", int(OmegaRec::PLOTS));
        MenuOpt[int(OmegaRec::EXIT)] = Form("%d. Exit.", int(OmegaRec::EXIT));

        break;
      }
      case 7:
      {
        MenuOpt[int(Efficiency::EFF_SCAN)] = Form("%d. Scan of the efficiency for the range of chosen cuts.", int(Efficiency::EFF_SCAN));
        MenuOpt[int(Efficiency::EFF_DIST)] = Form("%d. Distribution of the efficiency for the range of chosen cuts.", int(Efficiency::EFF_DIST));

        MenuOpt[int(Efficiency::EXIT)] = Form("%d. Exit.", int(Efficiency::EXIT));

        break;
      }
      case 8:
      {
        MenuOpt[int(Regen::REGEN_REJEC_TEST)] = Form("%d. Test of the Regeneration rejection variables.", int(Regen::REGEN_REJEC_TEST));
        MenuOpt[int(Regen::PLOTS)] = Form("%d. Plotting the results.", int(Regen::PLOTS));

        MenuOpt[int(Regen::EXIT)] = Form("%d. Exit.", int(Regen::EXIT));

        break;
      }
      case 9:
      {
        MenuOpt[int(Plots::PLOTS_GENERAL)] = Form("%d. General plots of the analysis.", int(Plots::PLOTS_GENERAL));
        MenuOpt[int(Plots::PLOTS_NORM)] = Form("%d. Plots for normalization constants.", int(Plots::PLOTS_NORM));

        MenuOpt[int(Plots::EXIT)] = Form("%d. Exit.", int(Plots::EXIT));

        break;
      }
      case 11:
      {
        MenuOpt[int(CovMatrix::USING_MC_DATA)] = Form("%d. Using MC Generated vs. Reconstructed.", int(CovMatrix::USING_MC_DATA));
        MenuOpt[int(CovMatrix::USING_CONTROL_SAMPLE)] = Form("%d. Using control samples.", int(CovMatrix::USING_CONTROL_SAMPLE));

        MenuOpt[int(CovMatrix::EXIT)] = Form("%d. Exit.", int(CovMatrix::EXIT));

        break;
      }
      case 12:
      {
        MenuOpt[int(MainMenuControlSample::KCHREC)] = Form("%d. Charged Kaon reconstruction.", int(MainMenuControlSample::KCHREC));
        MenuOpt[int(MainMenuControlSample::COV_MATRIX)] = Form("%d. Covariant Matrix Determination", int(MainMenuControlSample::COV_MATRIX));
        MenuOpt[int(MainMenuControlSample::CORR_FACTOR)] = Form("%d. Correction Factor Determination", int(MainMenuControlSample::CORR_FACTOR));

        MenuOpt[int(MainMenuControlSample::EXIT)] = Form("%d. Exit.", int(MainMenuControlSample::EXIT));

        break;
      }
      }
    };

    Menu(int ChooseMenu, std::vector<TString> fileNames) : ChooseMenu(ChooseMenu), fileOpt(fileNames)
    {
      switch (ChooseMenu)
      {
      case 0:
      {
        MenuOpt[int(KchRecMenu::K_MASS)] = Form("%d. Using Kaon mass as a benchmark.", int(KchRecMenu::K_MASS));
        MenuOpt[int(KchRecMenu::KSKL)] = Form("%d. Finding KS and then KL.", int(KchRecMenu::KSKL));
        MenuOpt[int(KchRecMenu::CLOSEST)] = Form("%d. Charged decay closest to the Bhabha IP.", int(KchRecMenu::CLOSEST));
        MenuOpt[int(KchRecMenu::BOOST)] = Form("%d. Correction of reconstructed momentum with boost method.", int(KchRecMenu::BOOST));
        MenuOpt[int(KchRecMenu::EXIT)] = Form("%d. Exit.", int(KchRecMenu::EXIT));
        break;
      }
      case 1:
      {
        MenuOpt[int(NeutRecMenu::BARE_TRILATERATION)] = Form("%d. Bare trilateration.", int(NeutRecMenu::BARE_TRILATERATION));
        MenuOpt[int(NeutRecMenu::TRILATERATION_KIN_FIT)] = Form("%d. Trilateration with kinematic fit.", int(NeutRecMenu::TRILATERATION_KIN_FIT));
        MenuOpt[int(NeutRecMenu::TRIANGLE)] = Form("%d. Trilateration with triangle.", int(NeutRecMenu::TRIANGLE));
        MenuOpt[int(NeutRecMenu::COMP_OF_MET)] = Form("%d. Comparison of methods, plotting, etc.", int(NeutRecMenu::COMP_OF_MET));
        MenuOpt[int(NeutRecMenu::EXIT)] = Form("%d. Exit.", int(NeutRecMenu::EXIT));

        break;
      }
      case 2:
      {
        MenuOpt[int(DataType::SIGNAL_TOT)] = Form("%d. Total MC signal.", int(DataType::SIGNAL_TOT));
        MenuOpt[int(DataType::SIG_BCG)] = Form("%d. Signal + background MC.", int(DataType::SIG_BCG));
        MenuOpt[int(DataType::MC_DATA)] = Form("%d. MC + Data.", int(DataType::MC_DATA));
        MenuOpt[int(DataType::SIGNAL_MAX)] = Form("%d. MC - without cuts and errors.", int(DataType::SIGNAL_MAX));

        break;
      }
      case 4:
      {
        MenuOpt[int(CPFitMenu::HALF_SIGNAL_MC)] = Form("%d. Only Signal; 1/2 MC vs. 1/2 'Data'", int(CPFitMenu::HALF_SIGNAL_MC));
        MenuOpt[int(CPFitMenu::HALF_SIG_BCG_MC)] = Form("%d. Signal + Bcg; 1/2 MC vs. 1/2 'Data'", int(CPFitMenu::HALF_SIG_BCG_MC));
        MenuOpt[int(CPFitMenu::MC_DATA)] = Form("%d. Signal + Bcg; MC vs. Data", int(CPFitMenu::MC_DATA));
        MenuOpt[int(CPFitMenu::EXIT)] = Form("%d. Exit.", int(CPFitMenu::EXIT));

        break;
      }
      case 5:
      {
        MenuOpt[int(GenVars::GEN_VARS)] = Form("%d. Generated variables to channel mapper.", int(GenVars::GEN_VARS));
        MenuOpt[int(GenVars::SPLIT_CHANN)] = Form("%d. Split stream into decay channels.", int(GenVars::SPLIT_CHANN));
        MenuOpt[int(GenVars::EXIT)] = Form("%d. Exit.", int(GenVars::EXIT));

        break;
      }
      case 6:
      {
        MenuOpt[int(OmegaRec::OMEGA_REC)] = Form("%d. Reconstruction of Omega channel based on the kinematic fit.", int(OmegaRec::OMEGA_REC));
        MenuOpt[int(OmegaRec::OMEGA_CUTS)] = Form("%d. Adjustment of cuts to reject Omega channel.", int(OmegaRec::OMEGA_CUTS));
        MenuOpt[int(OmegaRec::PLOTS)] = Form("%d. Plots, comparisons, etc.", int(OmegaRec::PLOTS));
        MenuOpt[int(OmegaRec::EXIT)] = Form("%d. Exit.", int(OmegaRec::EXIT));

        break;
      }
      case 7:
      {
        MenuOpt[int(Efficiency::EFF_SCAN)] = Form("%d. Scan of the efficiency for the range of chosen cuts.", int(Efficiency::EFF_SCAN));
        MenuOpt[int(Efficiency::EFF_DIST)] = Form("%d. Distribution of the efficiency for the range of chosen cuts.", int(Efficiency::EFF_DIST));

        MenuOpt[int(Efficiency::EXIT)] = Form("%d. Exit.", int(Efficiency::EXIT));

        break;
      }
      case 8:
      {
        MenuOpt[int(Regen::REGEN_REJEC_TEST)] = Form("%d. Test of the Regeneration rejection variables.", int(Regen::REGEN_REJEC_TEST));
        MenuOpt[int(Regen::PLOTS)] = Form("%d. Plotting the results.", int(Regen::PLOTS));

        MenuOpt[int(Regen::EXIT)] = Form("%d. Exit.", int(Regen::EXIT));

        break;
      }
      case 9:
      {
        MenuOpt[int(Plots::PLOTS_GENERAL)] = Form("%d. General plots of the analysis.", int(Plots::PLOTS_GENERAL));
        MenuOpt[int(Plots::PLOTS_NORM)] = Form("%d. Plots for normalization constants.", int(Plots::PLOTS_NORM));

        MenuOpt[int(Plots::EXIT)] = Form("%d. Exit.", int(Plots::EXIT));

        break;
      }
      case 11:
      {
        MenuOpt[int(CovMatrix::USING_MC_DATA)] = Form("%d. Using MC Generated vs. Reconstructed.", int(CovMatrix::USING_MC_DATA));
        MenuOpt[int(CovMatrix::USING_CONTROL_SAMPLE)] = Form("%d. Using control samples.", int(CovMatrix::USING_CONTROL_SAMPLE));

        MenuOpt[int(CovMatrix::EXIT)] = Form("%d. Exit.", int(CovMatrix::EXIT));

        break;
      }
      case 12:
      {
        MenuOpt[int(MainMenuControlSample::KCHREC)] = Form("%d. Charged Kaon reconstruction.", int(MainMenuControlSample::KCHREC));
        MenuOpt[int(MainMenuControlSample::COV_MATRIX)] = Form("%d. Covariant Matrix Determination", int(MainMenuControlSample::COV_MATRIX));
        MenuOpt[int(MainMenuControlSample::CORR_FACTOR)] = Form("%d. Correction Factor Determination", int(MainMenuControlSample::CORR_FACTOR));

        MenuOpt[int(MainMenuControlSample::EXIT)] = Form("%d. Exit.", int(MainMenuControlSample::EXIT));

        break;
      }
      }
    };

    void InitMenu() { std::cout << MenuName[ChooseMenu] << std::endl; };
    void EndMenu() { std::cout << ChooseOpt; };
    void ShowOpt()
    {
      switch (ChooseMenu)
      {
      case 0:
      {
        for (int i = 1; i < int(KchRecMenu::OPT_TOT); i++)
          std::cout << MenuOpt[i] << std::endl;

        break;
      }
      case 1:
      {
        for (int i = 1; i < int(NeutRecMenu::OPT_TOT); i++)
          std::cout << MenuOpt[i] << std::endl;

        break;
      }
      case 2:
      {
        for (int i = 1; i < int(DataType::OPT_TOT); i++)
        {
          std::cout << MenuOpt[i] << std::endl;
        }

        break;
      }
      case 3:
      {
        int iter = 0;
        for (const TString &fileName : fileOpt)
        {
          std::cout << iter << ". " << fileName << std::endl;
          iter++;
        }
      }
      case 4:
      {
        for (int i = 1; i < int(CPFitMenu::OPT_TOT); i++)
        {
          std::cout << MenuOpt[i] << std::endl;
        }

        break;
      }
      case 5:
      {
        for (int i = 1; i < int(GenVars::OPT_TOT); i++)
        {
          std::cout << MenuOpt[i] << std::endl;
        }

        break;
      }
      case 6:
      {
        for (int i = 1; i < int(OmegaRec::OPT_TOT); i++)
        {
          std::cout << MenuOpt[i] << std::endl;
        }

        break;
      }
      case 7:
      {
        for (int i = 1; i < int(Efficiency::OPT_TOT); i++)
        {
          std::cout << MenuOpt[i] << std::endl;
        }

        break;
      }
      case 8:
      {
        for (int i = 1; i < int(Regen::OPT_TOT); i++)
        {
          std::cout << MenuOpt[i] << std::endl;
        }

        break;
      }
      case 9:
      {
        for (int i = 1; i < int(Plots::OPT_TOT); i++)
        {
          std::cout << MenuOpt[i] << std::endl;
        }

        break;
      }
      case 11:
      {
        for (int i = 1; i < int(CovMatrix::OPT_TOT); i++)
        {
          std::cout << MenuOpt[i] << std::endl;
        }

        break;
      }
      case 12:
      {
        for (int i = 1; i < int(MainMenuControlSample::OPT_TOT); i++)
        {
          std::cout << MenuOpt[i] << std::endl;
        }

        break;
      }
      }
    }
  };
}

#endif