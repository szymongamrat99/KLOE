# KLOE Physics Analysis Framework (KLSPM00)

Particle physics analysis framework for the KLOE experiment, built on ROOT/RooFit and CERNLIB (Fortran) libraries. It processes `TChain`-based ROOT data through a set of kinematic-fitting and reconstruction stages (charged/neutral kaon reconstruction, kinematic fits for the signal, trilateration and omega hypotheses, covariance-matrix/efficiency studies, regeneration analysis, etc.) and produces ROOT output files, CSV summaries, and log reports.

## Overview

The executable `KLSPM00` reads KLOE ROOT ntuples via a `TChain`, then either:
- runs the **initial analysis** (`InitAnalysis_main`) that pre-processes raw files and prepares reconstructed variables, or
- runs the **interactive analysis menu**, which sequentially invokes the sub-analyses (`GenVars`, `KchRec`, `Neutrec`, `OmegaRec`, `CPFit`, `CovMatrix`, `Regen`, ...) declared in [Include/klspm00.hpp](Include/klspm00.hpp).

Behavior is driven by JSON configuration files under [Subanalysis/Properties/](Subanalysis/Properties) (analysis flags, cut limits, histogram definitions, PDG data, file paths), validated against JSON schemas.

## Repository Layout

```
app/                    Main entry point (CPVAnalysis.cpp) and its documentation
Include/                Shared headers and core framework code
  klspm00.hpp           Declarations of the *_main() entry points for each sub-analysis
  Codes/                Core framework: KinFitter engine, constraints, config/error/data managers
  MainMenuHandler/       Interactive menu loop and input handling
  FortranAnalysis/       Legacy Fortran/CERNLIB analysis code
Subanalysis/            One subdirectory per analysis stage (each builds its own static library)
  InitialAnalysis/       Pre-processing entry point (initanalysis_full.cpp)
  KchRec/, Neutrec/, OmegaRec/   Charged/neutral/omega reconstruction
  CPFit/, GeneratedVars/, CovarianceMatrix/, RegenerationAnalysis/, EfficiencyAnalysis/, ...
  Properties/            JSON configuration: analysis flags, cut limits, histogram configs, schemas, PDG data
ROOT_Macros/             Standalone ROOT macros for plotting/inspection
tests/                   Unit tests (CTest)
build/                   CMake build directory (generated, contains bin/KLSPM00)
execute_analysis.sh      Build + run a single analysis job
run_parallel.sh          Launch multiple jobs in parallel via GNU `parallel`
copy_libs.sh             Deploy built binary/libraries to the experiment software area
clean.sh                 `make clean` in the build directory
```

## Prerequisites

- CMake >= 3.13
- C++17 compiler and Fortran compiler (project mixes CXX and Fortran)
- [ROOT](https://root.cern) built with: `Core RIO Tree Physics MathCore Matrix TMVA Hist Graf RooFit RooFitCore`
- Boost (`system filesystem unit_test_framework`)
- CURL, LAPACK, BLAS

## Build

```bash
mkdir -p build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$(nproc)
```

The resulting binary is placed at `build/bin/KLSPM00`. Predefined VS Code tasks are also available: `build-release`, `build-debug`, `rebuild-debug`, `clean` (see [.vscode](.vscode)).

Debug builds enable `-Wall -Wextra -Wpedantic -fsanitize=address,undefined` (see [CMakeLists.txt](CMakeLists.txt)).

## Running

### Single job

```bash
export ANALYSIS_CONFIG_FILE=/path/to/analysis_config.json
./execute_analysis.sh <num_proc> <file_list.txt>
```

This rebuilds the project (Release) and runs `./bin/KLSPM00 <file_list.txt>` from the `build/` directory. `<file_list.txt>` lists the input ROOT files for the `TChain`.

### Interactive mode

When run without `initialAnalysisExec` forced in the config, `KLSPM00` prompts for:
1. Measurement (1) vs. Control Sample (2)
2. First/last file index (validated against `rootFiles.lastFileMax` in the config)
3. Data type (signal, background, MC, ...)

then executes the corresponding sub-analysis chain. See [app/CPVAnalysis.md](app/CPVAnalysis.md) for the full parameter/flow description.

### Standalone CP fit and systematic scans

Build the standalone executable with:

```bash
cmake --build build --target cpfit_standalone -j2
```

Set the same environment variables as for the main application:
`WORKDIR`, `KLOE_PROPERTIES_PATH`, `KLOE_LOGS_PATH`, `KLOE_DBV26_DK0`,
`KLOE_DBV26_MK0`, `ANALYSIS_CONFIG_FILE`, and `ROOT_RESULTS_PATH`. These are
used by the framework while libraries load, even if a particular input is not
used by the fit. `WORKDIR` must point to the repository. The files
`KLOE_PROPERTIES_PATH/paths-extensions.json` and
`KLOE_PROPERTIES_PATH/pdg_api/pdg_cache/pdg_cache_{S011,S012,S013}_2025.json`
must exist. Standalone initializes constants from this local PDG cache, as in
the menu, and does not download missing cache files.

The fit accepts any number of ROOT files containing the same tree:

```bash
build/Subanalysis/CPFit/cpfit_standalone \
  Subanalysis/CPFit/config/config.json result.json TREE_NAME file1.root file2.root
```

`result.json` must not already exist. The executable runs in batch mode, does
not produce plots or update global `properties.json`, and writes initialization
logs under `cpfit_logs/` next to the result. It exits nonzero on input,
configuration, or minimizer failure. Python (including a future Jupyter
notebook) can pass a list of paths without PyROOT:

```python
import subprocess

paths = ["file1.root", "file2.root"]
subprocess.run(
    ["build/Subanalysis/CPFit/cpfit_standalone",
     "Subanalysis/CPFit/config/config.json", "result.json", "TREE_NAME", *paths],
    check=True,
)
```

Edit a copy of [the CP-fit configuration](Subanalysis/CPFit/config/config.json)
without recompiling:

- `eventWeights.apply.signal`, `.semileptonic`, `.threePi0`, and
  `.regeneration` independently enable event corrections. A disabled correction
  uses weight 1; the corresponding Signal or Regeneration file is not opened.
  Signal corrections are configured independently for `ks_pm`, `kl_pm`,
  `ks_00`, and `kl_00` under `eventWeights.signalComponents`. Each enabled
  component contributes `ratio(R_ch) * correction(Delta_t_tri)` to its mode
  sum; if `ratioFile` is empty, the ratio factor is 1. The final Signal weight
  is `(ks_pm + kl_pm) * (ks_00 + kl_00)`, where each term denotes its complete
  ratio-times-correction weight. A mode with no enabled components contributes
  the identity factor 1. Set `useCorrectionSpline` per component to choose
  the configured `correctionSpline`; `false` uses ROOT
  histogram linear interpolation (`TH1::Interpolate`). Relative errors are
  always evaluated from the configured correction histogram. A correction
  spline evaluated outside its knot range (or at a non-finite coordinate)
  contributes `1`; ratio spline behavior is unchanged. The current
  `ks_pm` and `kl_pm` entries select their correction splines. The `ks_pm` and
  `kl_pm` components are enabled by default; `ks_00` and `kl_00` are disabled
  until their correction files are available; ratio data is optional.
  Correction and ratio
  file paths are relative to `Subanalysis/CPFit/` and also accept absolute
  paths. `signalMaxRelativeError` is checked independently for each enabled
  correction histogram; only a high-error component's histogram factor is set
  to 1, so that component still contributes its ratio factor (or 1 if omitted)
  and other components remain unchanged. `signalWeights.largeRelativeErrorReplacedByComponent`
  reports the number of substituted components. `regenerationShapeCorrection.weightFunctionsFile`
  is relative to `Subanalysis/RegenerationAnalysis/results/` and also accepts
  absolute paths.
- `eventWeights.crossChannelMeanWeights.threePi0FromKsPm` and
  `.semileptonicFromKs00` optionally multiply the existing 3pi0 and
  Semileptonic weights by the corresponding component's event-weighted mean
  over selected events of that same channel. Each 3pi0 event is weighted with
  `ks_pm`, and each Semileptonic event with `ks_00`, using that event's own
  `R_ch` and `Delta_t_tri`; the mean is taken after applying the configured
  relative-error substitution. The component's files are required even if its
  direct Signal weighting is disabled. Both switches default to `false`;
  calculated means are written to `crossChannelMeanWeights`.
- `diagnostics.signalWeightPlots.enabled` enables Signal correction-weight
  diagnostics for both standalone and menu runs. `weightBins` and `deltaTBins`
  set the histogram bin counts (1–1000); the weight axis is chosen from the
  data and the Δt axis uses `deltaT.rangeX`. Output alongside the result JSON:
  `_signal_weights_1d.pdf`, `_signal_weights_dt.pdf`, and
  `_signal_weights.root`. The 1D curves show the applied total and the two
  mode factors (`ks_pm + kl_pm`, `ks_00 + kl_00`), each including configured
  ratio factors. The 2D
  plot shows event counts versus `Dtboostlor` and applied total weight, with
  the mean per Δt bin overlaid; the ROOT file also contains the profile.
  When the entire event weight falls back to 1, both mode factors are recorded
  as 0 and `signal_weight_fallback` records the fallback. These are correction weights before interference, luminosity,
  branching-ratio and fitted normalization factors. ROOT underflow/overflow
  bins retain events outside the configured Δt range.
- `regenerationExclusion.enabled` controls the regeneration veto for both data
  and MC; it is disabled by default. The current veto uses the spherical means
  and sigmas together with `beamPipeBoundary` and `sigmaMultiplier`. The
  `cylindrical` settings are not used by this veto.
- `cuts.ellipse`, `cuts.trackAngle`, and `cuts.closeVertex` have `enabled`
  switches and numerical thresholds. Copy these sections when migrating an
  older CP-fit JSON. `deltaT` configures the histogram and fit range.
- `regenerationShapeCorrection.correctionFile` and `corrections` do not control
  the current CP-fit regeneration weights; the fitted functions from
  `weightFunctionsFile` are used instead.

For the RegenerationAnalysis fit, [fit_config.json](Subanalysis/RegenerationAnalysis/config/fit_config.json)
has independent `generalSettings.signalWeightComponents` switches for `ks_pm`,
`kl_pm`, `ks_00`, and `kl_00`. The Signal weight is the product of the
`ks_pm + kl_pm` mode sum and the `ks_00 + kl_00` mode sum; an empty mode uses
the identity factor 1. Each component has `correctionHistogram`,
`correctionSpline`, and `useCorrectionSpline`; the switch chooses spline
evaluation or histogram linear interpolation. Relative errors are always
measured from the histogram. A correction spline evaluated outside its knot
range (or at a non-finite coordinate) contributes `1`; ratio spline behavior
is unchanged. Each enabled component's `ratio(R_ch) * correction(Delta_t_tri)`
is added to its mode sum; an empty `ratioFile` skips ratio loading and uses 1.
The two `00` components are off and have empty correction paths until their
inputs are available.
`signalWeightRelativeErrorLimit`
uses a relative fraction (for example, `0.1` for 10%; `0` disables the error
cut). The error limit is checked independently for each enabled component;
only that component's correction histogram factor is set to 1, preserving its
ratio spline contribution and the other components. A non-finite/non-positive
composed weight is still replaced with 1. `isSignalWeighted` separately
controls multiplication by the interference weight.
`generalSettings.regenerationNormalization` selects how regeneration weights
are derived: `"luminosity"` (default) uses `(Data - MC background) / Regeneration`
with the configured `lumiFactor`; `"fit"` keeps Regeneration unchanged and
fits one background scale per histogram, `Data ~ Regeneration + a_B * MC background`
(bin-by-bin weighted least squares; background = all channels except
Regeneration). `generalSettings.backgroundFitBins` chooses the bins entering
the fit: `"nonRegeneration"` (default) excludes bins inside the regeneration
regions, `"all"` uses every bin. With `backgroundFitNonNegativeCap` (default `true`), `a_B` is capped at
`min_i(Data_i / Background_i)` over all bins so that `Data - a_B * Background`
has no negative bins; set it to `false` to remove the cap. The weight is
`(Data - a_B * Background) / Regeneration`; with
`divideWeightByBackgroundScale` (default `false`) it is additionally divided by
`a_B`, assuming `a_B` measures only the luminosity mismatch.
`applyBackgroundScaleToMCSum` (default `true`) scales the background channels of
the final 1D histograms and `MC sum` by `a_B` (Regeneration too, when the weights
are divided by `a_B`), so Data/MC plots and chi2 are consistent with the fit; the
time-difference and 2D histograms are not rescaled. If the fit fails, that
histogram falls back to the luminosity method. Fitted scales are printed.
`generalSettings.crossChannelMeanWeights.threePi0FromKsPm` and
`.semileptonicFromKs00` optionally scale the 3pi0 and Semileptonic histograms
by the corresponding component's mean weight over selected events in that
same channel. Each event's own `R_ch` and `Delta_t_tri` are used for
`ratioSpline(R_ch) * correction(Delta_t_tri)`, using the configured histogram
interpolation or correction spline and the same
relative-error substitution. The required component files must be configured
even if that component is disabled for direct Signal weighting. Both switches
default to `false`.

The JSON result contains the complete configuration, initialized physics
constants, branching-ratio corrections, luminosity factor, fit status, chi-square,
parameter errors and correlations, selected counts and fitted integrals by
channel, and Signal weight sum, sum of squares and replacement counts.
Disabling event weights does **not** disable the independent branching-ratio or
luminosity corrections. The Signal template is normalized to the sum of its
event weights before the fitted Signal amplitude is applied; that amplitude
alone is therefore not a direct measurement of the average weight.

For a one-parameter scan, use a **new** output directory:

```bash
python3 Subanalysis/CPFit/app/scan_cuts.py \
  --runner build/Subanalysis/CPFit/cpfit_standalone \
  --config Subanalysis/CPFit/config/config.json \
  --output-dir ellipse_scan --tree TREE_NAME \
  --key cuts.ellipse.nSigma --values '[2.0, 2.5, 3.0]' file1.root file2.root
```

The key may also name a boolean switch, for example
`--key eventWeights.apply.signal --values '[false, true]'`. Each variant runs
in a separate process and rereads all input events. The directory contains a
manifest, `summary.csv`, and, for each variant, its full configuration,
`result.json`, and `run.log`. Failed variants remain in the summary with a
nonzero return code; the scanner also returns nonzero. The scanner varies one
key at a time; it does not calculate systematic uncertainties automatically.

### Parallel batch jobs

```bash
./run_parallel.sh <analysis_type> <date> <first_idx> <last_idx> <num_proc> <config_file>
```

Uses GNU `parallel` to run `execute_analysis.sh` for a range of job indices in the background (`nohup`), writing logs to `parallel_logs_<type>_<first>_<last>/` and `nohup_<type>_<first>_<last>.log`.

## Configuration

Analysis behavior is controlled by JSON files in [Subanalysis/Properties/](Subanalysis/Properties):
- `properties.json` — global flags (`analysisCode`, `covMatrixType`, `initialAnalysisExec`, `trilaterationReconstructionMode`, momentum-smearing matrices, ...)
- `analysis_config*.json` — per-mode/per-dataset configs (Data/MC, Signal/3pi0/Semileptonic), validated against `schemas/`
- `cut-limits.json` / `cut-limits-final.json` — selection cut values
- `histogram_conf/histogram1D.csv`, `histogram2D.csv` — histogram definitions parsed by `KLOE::Histograms::LoadHistogramConfigs1D/2D()`
- `root-files.json`, `paths-extensions.json` — input/output file path resolution
- `pdg_api/` — cached PDG particle data

Select the active config via the `ANALYSIS_CONFIG_FILE` environment variable before running.

## Core Framework

- **KinFitter** ([Include/Codes/src/KinFitter.cpp](Include/Codes/src/KinFitter.cpp)) — Newton–Raphson based kinematic fitting engine shared by all fit hypotheses; produces fitted parameters `_X` and `_chiSqr`.
- **Constraints** (`ConstraintsSignal`, `ConstraintsTrilateration`, `ConstraintsOmega`, ...) — hypothesis-specific constraint sets derived from `ConstraintsBase`.
- **Analysis entry points** ([Subanalysis/InitialAnalysis/src/initanalysis_full.cpp](Subanalysis/InitialAnalysis/src/initanalysis_full.cpp)) — runs `SignalKinFit::Reconstruct()` → `TrilerationKinFit::Reconstruct()` → `OmegaKinFit::Reconstruct()` and collects results.
- **Error handling/logging** — `ErrorHandling::ErrorLogs`, used throughout via try/catch around risky operations (matrix inversion, file I/O, fits); logs are written under [log/](log/).

See [.github/copilot-instructions.md](.github/copilot-instructions.md) for detailed coding conventions, ROOT cross-version (6.24 → 6.36+) compatibility notes, and common pitfalls.

## Testing

```bash
cd build
ctest
```

Tests live in [tests/](tests) (e.g. `test_charged_vtx_rec.cpp`) and are registered via CTest when the `tests` subdirectory is enabled in [CMakeLists.txt](CMakeLists.txt).

The standalone CP-fit integration test creates synthetic ROOT inputs in a
temporary directory. After building `cpfit_standalone` and setting the
environment and local PDG cache described above, run:

```bash
python3 tests/test_cpfit_standalone.py
```

## Debugging

```bash
# Memory leak detection
valgrind --leak-check=full ./bin/KLSPM00

# AddressSanitizer / UndefinedBehaviorSanitizer (Debug build)
cmake -DCMAKE_BUILD_TYPE=Debug ..
export ASAN_OPTIONS="abort_on_error=1:detect_leaks=1"

# Segfault with full stack trace
gdb --args ./bin/KLSPM00
```

## Deployment

`copy_libs.sh` copies the built binary and all `*.so` libraries to the shared experiment software area (`/opt/exp_software/kloe/users/.../KLSPM00/`). It runs automatically at the end of `execute_analysis.sh` when executed on host `ui-tier1.cr.cnaf.infn.it`.

Czynnik świetlności MC → Data ustawia `luminosityFactor` na górze `config.json` (domyślnie 1.07906 przy braku klucza); trafia do JSON wyniku.

`signalNormalizationInRangeOnly` (domyślnie `false`): `true` normalizuje Signal do sumy wag tylko zdarzeń z Δt w zakresie fitu (`deltaT.rangeX`), zamiast wszystkich wybranych zdarzeń.

`branchingRatioCorrection`: obiekt `{kanał: true/false}` (kanały: Signal, Regeneration, Omega, 3pi0, Semileptonic, Other; brak klucza = włączone); `false` daje czynnik 1 i zerowy błąd BR dla tego kanału. Pojedyncze `true`/`false` działa dla wszystkich kanałów. JSON zawiera `branchingRatioCorrectionApplied` per kanał.

### Diagnostyka fitu CP

Plik wynikowy JSON zawiera sekcję `diagnostics`:
- `parameters.<nazwa>.atLowerBound/atUpperBound` – parametr na granicy limitu;
- `degeneracy` – wartości własne macierzy korelacji, `conditionNumber`, `weakestDirection` (składowe najsłabiej wyznaczonej kombinacji parametrów) i pary o `|corr| > 0.9`;
- `profile.<parametr>` – profil chi2 (`deltaChi2` w funkcji wartości, pozostałe parametry reminimalizowane); włączany w konfiguracji: `"diagnostics": {"profile": {"enabled": true, "points": 11, "nSigma": 3.0}}`.

Gdy względny błąd binu histogramu składnika Signal (`GetBinError/GetBinContent`) przekracza `eventWeights.signalMaxRelativeError` (domyślnie w `config.json` 0.1; 0 wyłącza), tylko czynnik tego histogramu jest zastępowany przez 1; jego składnik nadal wynosi `ratio × 1`, a inne komponenty pozostają bez zmian. JSON podaje liczbę zdarzeń z taką zamianą w `signalWeights.largeRelativeErrorReplaced` oraz liczniki per komponent w `signalWeights.largeRelativeErrorReplacedByComponent`.

`eventWeights.crossChannelMeanWeights.threePi0FromKsPm` oraz `.semileptonicFromKs00` włączają dodatkowe skalowanie odpowiednio 3pi0 średnią wagą `ks_pm` i Semileptonic średnią wagą `ks_00`. Dla każdego zdarzenia danego kanału liczony jest jego własny iloczyn `ratio × correction` z jego `R_ch` i `Δt_tri`, a następnie średnia po zaakceptowanych zdarzeniach tego kanału. Przełączniki domyślnie są wyłączone; wyliczone wartości trafiają do `crossChannelMeanWeights` w JSON.
