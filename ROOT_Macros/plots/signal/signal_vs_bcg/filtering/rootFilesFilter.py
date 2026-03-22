import ROOT
import csv
import os
import json
from modules import utilities as utils

# Enable ROOT's implicit multi-threading to speed up the processing of large datasets.
ROOT.EnableImplicitMT()

# Include necessary headers for constants and error logging.
ROOT.gInterpreter.Declare('#include <const.h>')
ROOT.gInterpreter.Declare('#include <ErrorLogs.h>')
ROOT.gInterpreter.Declare('#include <kloe_class.h>')

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
    filesMC[key] = "mk0*" + chann + "*.root"

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


tmpiter = 0
lastCount = 0
lastCode = 0

for code in sorted(errorCodes.values()):
    if tmpiter == 0:
      lastCount = totSignalEvents
    else:
      lastCount = totSignalEventsErrors[lastCode]

    tmpEvents = int(rdfErrorsSignal.Filter("errorcode == {}".format(code)).Count().GetValue())

    if tmpEvents > 0:
      totSignalEventsErrors[code] = lastCount - tmpEvents
      lastCode = code
      tmpiter = 1

# Filter for good clusters
numGoodClustersAll = rdf[1].Filter("goodClustersTriKinFitSize == 4").Count().GetValue()
numGoodClustersThreeOrAll = rdf[1].Filter("goodClustersTriKinFitSize >= 3").Count().GetValue()

# Print all the events and percentages after each error code
print("Total signal events (signal + errors): {}".format(totSignalEvents))

csvOutputPath = os.path.join(rootFileDatedFolder, "signal_error_summary.csv")
csvRows = []

for code, count in sorted(totSignalEventsErrors.items()):
    percentage = round((float(count) / float(totSignalEvents)) * 100, 2)
    print("Error code {} ({}): {} events ({}%)".format(reversedErrorCodes[code], code, count, percentage))
    csvRows.append([
      reversedErrorCodes[code],
      code,
      count,
      "{:.2f}".format(percentage)
    ])

print("Number of events with at least 3 good clusters: {}".format(numGoodClustersThreeOrAll))
print("Number of events with exactly 4 good clusters: {}".format(numGoodClustersAll))

with open(csvOutputPath, "wb") as csvFile:
  csvWriter = csv.writer(csvFile)
  csvWriter.writerow(["error_name", "error_code", "events", "percentage"])
  csvWriter.writerow(["TOTAL_SIGNAL_EVENTS", "", totSignalEvents, "100.00"])
  csvWriter.writerows(csvRows)

print("CSV saved to {}".format(csvOutputPath))

