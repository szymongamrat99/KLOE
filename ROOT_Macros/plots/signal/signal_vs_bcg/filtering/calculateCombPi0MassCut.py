import ROOT
import csv
import os
import json
from modules import utilities as utils
import math
import numpy as np

ROOT.gInterpreter.Declare('''
#include <TError.h>
#include <cstring>

namespace KloePyRoot {
  void SelectiveErrorHandler(int level, Bool_t abort, const char *location, const char *msg)
  {
    if (level >= kError && location && msg)
    {
      const bool noisyListCleanup =
          ((std::strstr(location, "TList::Clear") != nullptr) ||
           (std::strstr(location, "THashList::Delete") != nullptr)) &&
          (std::strstr(msg, "already deleted") != nullptr);

      if (noisyListCleanup)
        return;
    }

    DefaultErrorHandler(level, abort, location, msg);
  }

  void EnableSelectiveErrorHandler()
  {
    SetErrorHandler(SelectiveErrorHandler);
  }
}
''')

ROOT.KloePyRoot.EnableSelectiveErrorHandler()


def find_repo_root(start_dir):
  """Walk up the directory tree and find repository root containing Include/Codes/inc."""
  current = os.path.abspath(start_dir)
  while True:
    if os.path.isdir(os.path.join(current, "Include", "Codes", "inc")):
      return current
    parent = os.path.dirname(current)
    if parent == current:
      raise RuntimeError("Could not locate repository root with Include/Codes/inc")
    current = parent

# Enable ROOT's implicit multi-threading to speed up the processing of large datasets.
ROOT.EnableImplicitMT()

# Ensure ROOT's interpreter can resolve project headers included below.
repo_root = find_repo_root(os.path.dirname(__file__))
ROOT.gInterpreter.AddIncludePath(os.path.join(repo_root, "Include", "Codes", "inc"))

# Include necessary headers for constants and error logging.
ROOT.gInterpreter.Declare('#include <const.h>')
ROOT.gInterpreter.Declare('#include <ErrorLogs.h>')
ROOT.gInterpreter.Declare('#include <kloe_class.h>')
ROOT.gInterpreter.Declare('#include <double_gaus.h>')

logger = ROOT.ErrorHandling.ErrorLogs("log/")

ROOT.KLOE.setGlobalStyle()
ROOT.Utils.InitializeVariables(logger)

pm00 = ROOT.KLOE.pm00()

rootFileDatedFolder = 'root_files/' + str(pm00.getCurrentDate())

if not os.path.exists(rootFileDatedFolder):
  os.makedirs(rootFileDatedFolder)

channName = dict(ROOT.KLOE.channName)

configPath = "config/variables.json"
with open(configPath, "r") as configFile:
  config = json.load(configFile)
#######################################################################################

# config1D = utils.parse_histogram_configs(str(ROOT.Paths.histogramConfig1DPath))
# config2D = utils.parse_histogram_configs(str(ROOT.Paths.histogramConfig2DPath))

#######################################################################################

filePaths = config["filePathsMC"]

filesMC = {}
for key, chann in channName.items():
  if chann == 'Signal': #chann != 'Data' and chann != 'MC sum' and chann != 'pi+pi-pi+pi-':
    filesMC[key] = "mk0*.root"

filesChain = {}
for key, file in filesMC.items():
  sublist = []
  for path in filePaths:
    sublist.append(path + "/" + file)
  filesChain[key] = sublist
  
# Initialize RDataFrame
rdf = {}
for key, file in filesChain.items():
  rdf[key] = ROOT.RDataFrame("h1", file)
  rdf[key] = rdf[key].Define("mctruth_int", "mctruth == 0 ? " + str(key) + " : mctruth")

columnsToSave = ["mcflag", "errorcode", "cutApplied", "mctruth_int", "mctruth", "Chi2SignalKinFit", "Chi2OmegaKinFit", "Chi2TriKinFit", "Bx", "By", "Bz", "Bpx", "Bpy", "Bpz", "Broots", "Cotv1", "Cotv2", "Curv1", "Curv2", "Phiv1", "Phiv2", "CotvSmeared1", "CotvSmeared2", "CurvSmeared1", "CurvSmeared2", "PhivSmeared1", "PhivSmeared2", "Qmiss", "TrcSum", "minv4gam", "bestErrorSixGamma", "CotvMC", "CurvMC", "PhivMC", "Kchboost", "KchboostFit", "Kchmc", "Kchrec", "KchrecFit", "Knemc", "Knerec", "KnerecFit", "KnereclorFit", "KnetriKinFit", "gammaMomTriKinFit1", "gammaMomTriKinFit2", "gammaMomTriKinFit3", "gammaMomTriKinFit4", "gammaMomTriangle1", "gammaMomTriangle2", "gammaMomTriangle3", "gammaMomTriangle4", "ip", "ipFit", "ipOmegaFit", "ipTriKinFit", "ipmc", "neuVtxTriKinFit", "omega", "omegaFit", "phiOmegaFit", "photonFit1", "photonFit2", "photonFit3", "photonFit4", "photonOmegaFit1", "photonOmegaFit2", "photonOmegaFit3", "photonOmegaFit4", "pi01", "pi02", "pi01Fit", "pi02Fit", "pi0Omega1", "pi0Omega2", "pi0OmegaFit1", "pi0OmegaFit2", "pullsSignalFit", "trcfinal", "trk1", "trk2", "trk1Fit", "trk2Fit", "trk1MC", "trk2MC", "trkOmegaFit1", "trkOmegaFit2", "KnerecSix"]

hist1D = {}

# Filtered only to proper events
for key, frame in rdf.items():
  if key == 1:
    rdfErrorsSignal = frame.Filter("mctruth_int == -1")

  rdf[key] = frame.Filter("mctruth_int == " + str(key))

# Calculation of the errors w.r.t. the signal events

totSignalEventsNoErrors = int(rdf[1].Count().GetValue())
totSignalEvents = totSignalEventsNoErrors + int(rdfErrorsSignal.Count().GetValue())

totSignalEventsErrors = {}

errorCodes = {
    "NO_VTX_WITH_TWO_TRACKS":                   300,
    "LESS_THAN_FOUR_NEUTRAL_CLUSTERS":           301,
    "LESS_THAN_SIX_NEUTRAL_CLUSTERS":            302,
    "NO_VTX_WITH_OPPOSITE_TRACKS":               303,
    "LESS_THAN_FOUR_CLUSTERS_WITH_GOOD_ENERGY":  304,
    "NO_TWO_VTX_WITH_TWO_TRACKS":                305,
    "CHARGED_KAON_MASS_PRE":                     306,
    "TRILATERATION_KIN_FIT":                     307,
    "NO_VALID_SIX_GAMMA_SOLUTION":               308,
    "TRIANGLE_REC":                              309,
    "SIGNAL_KIN_FIT":                            310,
    "OMEGA_KIN_FIT":                             311,
}

reversedErrorCodes = {v: k for k, v in errorCodes.items()}

rdf_sig = rdf[1].Filter("Chi2SignalKinFit < 30").Define("Pi0Mass1", "pi01Fit[5]").Define("Pi0Mass2", "pi02Fit[5]")

modelComb = ROOT.RDF.TH1DModel(
  "hCombpi0Mass",
  "Combined pi0 mass error;#sqrt{(m_{#pi^{0}_{1}}-m_{#pi^{0}})^{2} + (m_{#pi^{0}_{2}}-m_{#pi^{0}})^{2}} [MeV/c^{2}];Counts",
  100,
  0,
  15,
)

modelPi0 = ROOT.RDF.TH1DModel(
  "hPi0Mass",
  "Pi0 mass error; m^{inv,fit}_{2#gamma} [MeV/c^{2}];Counts",
  100,
  80,
  190,
)

modelPi02D = ROOT.RDF.TH2DModel(
  "hPi0Mass2D",
  "Pi0 mass error; m^{inv,fit}_{2#gamma,1} [MeV/c^{2}]; m^{inv,fit}_{2#gamma,2} [MeV/c^{2}]",
  100,
  80,
  190,
  100,
  80,
  190,
)

hist1 = rdf_sig.Histo1D(modelPi0, "Pi0Mass1")
hist2 = rdf_sig.Histo1D(modelPi0, "Pi0Mass2")

hist2d = rdf_sig.Histo2D(modelPi02D, "Pi0Mass1", "Pi0Mass2")

f1 = ROOT.TF1('f_double_1', "[0] * exp(-0.5*((x-[1])/[2])^2) + [3] * exp(-0.5*((x-[1])/[4])^2)", 80.0, 190.0, 5)
f1.SetParNames('A1', 'mean', 'sigma1', 'A2', 'sigma2')
f1.SetParameters(hist1.GetMaximum()*0.7, hist1.GetMean(), hist1.GetRMS(), hist1.GetMaximum()*0.3, hist1.GetRMS() * 4)
f1.SetParLimits(0, 0.0, hist1.GetMaximum() * 2)
f1.SetParLimits(2, 0.1, 30.0)
f1.SetParLimits(3, 0.0, hist1.GetMaximum() * 2)
f1.SetParLimits(4, 0.1, 30.0)

f2 = ROOT.TF1('f_double_2', "[0] * exp(-0.5*((x-[1])/[2])^2) + [3] * exp(-0.5*((x-[1])/[4])^2)", 80.0, 190.0, 5)
f2.SetParNames('A1', 'mean', 'sigma1', 'A2', 'sigma2')
f2.SetParameters(hist2.GetMaximum()*0.7, hist2.GetMean(), hist2.GetRMS(), hist2.GetMaximum()*0.3, hist2.GetRMS() * 4)
f2.SetParLimits(0, 0.0, hist2.GetMaximum() * 2)
f2.SetParLimits(2, 0.1, 30.0)
f2.SetParLimits(3, 0.0, hist2.GetMaximum() * 2)
f2.SetParLimits(4, 0.1, 30.0)

hist1.Fit(f1, "SRI")
hist2.Fit(f2, "SRI")

f1part1 = ROOT.TF1('f_double_1_part1', "[0] * exp(-0.5*((x-[1])/[2])^2)", 80.0, 190.0, 3)
f1part1.SetParNames('A1', 'mean', 'sigma1')
f1part1.SetParameters(f1.GetParameter(0), f1.GetParameter(1), f1.GetParameter(2))
f1part1.SetLineColor(ROOT.kGreen)
f1part1.SetLineStyle(ROOT.kDashed)

f1part2 = ROOT.TF1('f_double_1_part2', "[0] * exp(-0.5*((x-[1])/[2])^2)", 80.0, 190.0, 3)
f1part2.SetParNames('A2', 'mean', 'sigma2')
f1part2.SetParameters(f1.GetParameter(3), f1.GetParameter(1), f1.GetParameter(4))
f1part2.SetLineColor(ROOT.kBlue)
f1part2.SetLineStyle(ROOT.kDashed)

f2part1 = ROOT.TF1('f_double_2_part1', "[0] * exp(-0.5*((x-[1])/[2])^2)", 80.0, 190.0, 3)
f2part1.SetParNames('A1', 'mean', 'sigma1')
f2part1.SetParameters(f2.GetParameter(0), f2.GetParameter(1), f2.GetParameter(2))
f2part1.SetLineColor(ROOT.kGreen)
f2part1.SetLineStyle(ROOT.kDashed)

f2part2 = ROOT.TF1('f_double_2_part2', "[0] * exp(-0.5*((x-[1])/[2])^2)", 80.0, 190.0, 3)
f2part2.SetParNames('A2', 'mean', 'sigma2')
f2part2.SetParameters(f2.GetParameter(3), f2.GetParameter(1), f2.GetParameter(4))
f2part2.SetLineColor(ROOT.kBlue)
f2part2.SetLineStyle(ROOT.kDashed)

print("Pi0Mass1 fit parameters:")
print("A1 =", f1.GetParameter(0))
print("mean =", f1.GetParameter(1))
print("sigma1 =", f1.GetParameter(2))
print("A2 =", f1.GetParameter(3))
print("sigma2 =", f1.GetParameter(4))

print("Pi0Mass2 fit parameters:")
print("A1 =", f2.GetParameter(0))
print("mean =", f2.GetParameter(1))
print("sigma1 =", f2.GetParameter(2))
print("A2 =", f2.GetParameter(3))
print("sigma2 =", f2.GetParameter(4))

sigmas1 = [f1.GetParameter(2), f1.GetParameter(4)]
sigmasErr1 = [f1.GetParError(2), f1.GetParError(4)]
sigmas2 = [f2.GetParameter(2), f2.GetParameter(4)]
sigmasErr2 = [f2.GetParError(2), f2.GetParError(4)]

minSigma1 = min(sigmas1)
minSigmaErr1 = sigmasErr1[sigmas1.index(minSigma1)]
minSigma2 = min(sigmas2)
minSigmaErr2 = sigmasErr2[sigmas2.index(minSigma2)]

tailSigma1 = max(sigmas1)
tailSigmaErr1 = sigmasErr1[sigmas1.index(tailSigma1)]
tailSigma2 = max(sigmas2)
tailSigmaErr2 = sigmasErr2[sigmas2.index(tailSigma2)]

print("Minimum mean for Pi0Mass1: ", f1.GetParameter(1), "+/-", f1.GetParError(1))
print("Minimum sigma for Pi0Mass1: ", minSigma1, "+/-", minSigmaErr1)
print("Minimum mean for Pi0Mass2: ", f2.GetParameter(1), "+/-", f2.GetParError(1))
print("Minimum sigma for Pi0Mass2: ", minSigma2, "+/-", minSigmaErr2)

meanPi0Mass1 = f1.GetParameter(1)
meanPi0Mass2 = f2.GetParameter(1)

rho = hist2d.GetCorrelationFactor()
print("Correlation factor between Pi0Mass1 and Pi0Mass2: ", rho)

## Normalized radius

dx = "(Pi0Mass1 - {})/{}".format(meanPi0Mass1, minSigma1)
dy = "(Pi0Mass2 - {})/{}".format(meanPi0Mass2, minSigma2)

rho_factor = "1 / (1 - {})".format(rho**2)

rnorm = "sqrt({}*(pow({}, 2) + pow({}, 2) - 2*{}*{}*{}))".format(rho_factor, dx, dy, rho, dx, dy)

## Square cut

u = "((Pi0Mass1 - {}) + (Pi0Mass2 - {})) / sqrt(2)".format(meanPi0Mass1, meanPi0Mass2)
v = "((Pi0Mass1 - {}) - (Pi0Mass2 - {})) / sqrt(2)".format(meanPi0Mass1, meanPi0Mass2)

varu = "0.5 * ({}^2 + {}^2 + 2*{}*{}*{})".format(minSigma1, minSigma2, rho, minSigma1, minSigma2)
varv = "0.5 * ({}^2 + {}^2 - 2*{}*{}*{})".format(minSigma1, minSigma2, rho, minSigma1, minSigma2)

sigmau = "sqrt({})".format(varu)
sigmav = "sqrt({})".format(varv)

# Numeric values for quick logging (the string expressions above are used in RDF filters).
sigmau_val = math.sqrt(0.5 * (minSigma1**2 + minSigma2**2 + 2 * rho * minSigma1 * minSigma2))
sigmav_val = math.sqrt(0.5 * (minSigma1**2 + minSigma2**2 - 2 * rho * minSigma1 * minSigma2))

print("sigmau expression:", sigmau)
print("sigmav expression:", sigmav)
print("sigmau value =", sigmau_val)
print("sigmav value =", sigmav_val)

squareCut = "std::abs({}) < 3 * {} && std::abs({}) < 3 * {}".format(u, sigmau, v, sigmav)

## Display the rectangle

u_max = 3 * math.sqrt(0.5 * (minSigma1**2 + minSigma2**2 + 2 * rho * minSigma1 * minSigma2))
v_max = 3 * math.sqrt(0.5 * (minSigma1**2 + minSigma2**2 - 2 * rho * minSigma1 * minSigma2))

x_coords = []
y_coords = []

for u_sign, v_sign in [(1,1), (-1,1), (-1,-1), (1,-1), (1,1)]:
  u_cut = u_sign * u_max
  v_cut = v_sign * v_max

  x_coords.append(meanPi0Mass1 + (u_cut + v_cut) / math.sqrt(2))
  y_coords.append(meanPi0Mass2 + (u_cut - v_cut) / math.sqrt(2))

x_coords = np.array(x_coords)
y_coords = np.array(y_coords)

cut_rectangle = ROOT.TPolyLine(5, x_coords, y_coords)
cut_rectangle.SetLineColor(ROOT.kRed)
cut_rectangle.SetLineWidth(3)
cut_rectangle.SetLineStyle(2) # Linia przerywana


fitParameters1 = ROOT.TPaveText(0.6, 0.65, 1.0, 0.85, "NDC")
fitParameters1.SetBorderSize(1)
fitParameters1.SetFillColor(ROOT.kWhite)
fitParameters1.SetTextSize(0.03)
fitParameters1.AddText("#mu(m): {:.2f} #pm {:.2f} MeV/c^{{2}}".format(meanPi0Mass1, f1.GetParError(1)))
fitParameters1.AddText("#sigma(m)_{{core}}: {:.2f} #pm {:.2f} MeV/c^{{2}}".format(minSigma1, minSigmaErr1))
fitParameters1.AddText("#sigma(m)_{{tail}}: {:.2f} #pm {:.2f} MeV/c^{{2}}".format(tailSigma1, tailSigmaErr1))

fitParameters2 = ROOT.TPaveText(0.6, 0.65, 1.0, 0.85, "NDC")
fitParameters2.SetBorderSize(1)
fitParameters2.SetFillColor(ROOT.kWhite)
fitParameters2.SetTextSize(0.03)
fitParameters2.AddText("#mu(m): {:.2f} #pm {:.2f} MeV/c^{{2}}".format(meanPi0Mass2, f2.GetParError(1)))
fitParameters2.AddText("#sigma(m)_{{core}}: {:.2f} #pm {:.2f} MeV/c^{{2}}".format(minSigma2, minSigmaErr2))
fitParameters2.AddText("#sigma(m)_{{tail}}: {:.2f} #pm {:.2f} MeV/c^{{2}}".format(tailSigma2, tailSigmaErr2))


fitParameters = ROOT.TPaveText(0.6, 0.65, 1.0, 0.85, "NDC")
fitParameters.SetBorderSize(1)
fitParameters.SetFillColor(ROOT.kWhite)
fitParameters.SetTextSize(0.03)
fitParameters.AddText("#mu(m)_{{core,1}}: {:.2f} #pm {:.2f} MeV/c^{{2}}".format(meanPi0Mass1, f1.GetParError(1)))
fitParameters.AddText("#sigma(m)_{{core,1}}: {:.2f} #pm {:.2f} MeV/c^{{2}}".format(minSigma1, minSigmaErr1))
fitParameters.AddText("#mu(m)_{{core,2}}: {:.2f} #pm {:.2f} MeV/c^{{2}}".format(meanPi0Mass2, f2.GetParError(1)))
fitParameters.AddText("#sigma(m)_{{core,2}}: {:.2f} #pm {:.2f} MeV/c^{{2}}".format(minSigma2, minSigmaErr2))
fitParameters.AddText("Correlation factor: {:.2f}".format(rho))


rdf_sig = rdf_sig.Define("CombPi0MassError", rnorm)

histComb = rdf_sig.Histo1D(modelComb, "CombPi0MassError")

TLine = ROOT.TLine(3., 0, 3., histComb.GetMaximum())
TLine.SetLineColor(ROOT.kRed)
TLine.SetLineStyle(ROOT.kDashed)

c1 = ROOT.TCanvas("c1", "c1", 800, 800)

histComb.SetStats(0)
histComb.SetLineColor(ROOT.kBlack)
histComb.SetTitle("Combined pi0 mass error")
histComb.GetXaxis().SetTitle("R_{mass}^{#pi^{0}} [MeV/c^{2}]")
histComb.GetYaxis().SetTitle("Counts")
histComb.Draw()
TLine.Draw("same")
fitParameters.Draw("same")
c1.SaveAs(rootFileDatedFolder + "/pi0MassCutSignal.svg")

c2 = ROOT.TCanvas("c2", "c2", 800, 800)
hist1.SetStats(0)
hist1.SetLineColor(ROOT.kBlack)
hist1.SetTitle("Pi0 mass error for signal events")
hist1.GetXaxis().SetTitle("m^{inv,fit}_{2#gamma,1} [MeV/c^{2}]")
hist1.GetYaxis().SetTitle("Counts")
hist1.Draw()
f1.Draw("same")
f1part1.Draw("same")
f1part2.Draw("same")
fitParameters1.Draw("same")
c2.SaveAs(rootFileDatedFolder + "/pi0Mass1.svg")

c3 = ROOT.TCanvas("c3", "c3", 800, 800)
hist2.SetStats(0)
hist2.SetLineColor(ROOT.kBlack)
hist2.SetTitle("Pi0 mass error for signal events")
hist2.GetXaxis().SetTitle("m^{inv,fit}_{2#gamma,2} [MeV/c^{2}]")
hist2.GetYaxis().SetTitle("Counts")
hist2.Draw()
f2.Draw("same")
f2part1.Draw("same")
f2part2.Draw("same")
fitParameters2.Draw("same")
c3.SaveAs(rootFileDatedFolder + "/pi0Mass2.svg")

c4 = ROOT.TCanvas("c4", "c4", 800, 800)
hist2d.SetStats(0)
hist2d.SetTitle("Pi0 mass error for signal events")
hist2d.GetXaxis().SetTitle("m^{inv,fit}_{2#gamma,1} [MeV/c^{2}]")
hist2d.GetYaxis().SetTitle("m^{inv,fit}_{2#gamma,2} [MeV/c^{2}]")
hist2d.Draw("COLZ")
cut_rectangle.Draw("same")
c4.SaveAs(rootFileDatedFolder + "/pi0Mass2D.svg")

