# Copilot Instructions for KLOE Physics Analysis Framework

## Project Overview
KLOE physics analysis framework for particle physics data processing using ROOT, RooFit, and custom kinematic fitting algorithms. The project combines:
- **Data Processing**: ROOT TChain analysis with TSelector-based event processing
- **Kinematic Fitting**: Custom constraint-based fitting (Signal, Trilateration, Omega modes)
- **Physics Reconstruction**: Neutral/charged particle reconstruction, PCA correlation analysis
- **Visualization**: ROOT histogram management with statistical analysis and fit result display

## Architecture Quick Reference

### Core Components
1. **KinFitter Framework** (`Include/Codes/src/KinFitter.cpp`)
   - Base class for kinematic fitting with Newton-Raphson optimization
   - Critical: Uses `_constraints` vector of `TF1` objects (ROOT memory management required)
   - Output: `_X` (fitted parameters), `_chiSqr` (optimization metric)
   - Performance bottleneck: ~600-800ms per event (matrix operations, gradient calculations)

2. **Constraint Classes** (`Include/Codes/src/ConstraintsSignal.cpp`, etc.)
   - Derived from `ConstraintsBase` or standalone
   - Key methods: `IntermediateReconstruction(Double_t *p)` - must initialize all member vectors
   - Common crash: `std::vector<float>::operator=` in reconstruction → validate vector sizes

3. **Analysis Entry Points** (`Subanalysis/InitialAnalysis/src/initanalysis_full.cpp`)
   - Sequential execution: `SignalKinFit::Reconstruct()` → `TrilerationKinFit::Reconstruct()` → `OmegaKinFit::Reconstruct()`
   - Critical: `GetResults()` copies fitted parameters to output structures
   - Thread safety: Direct member access preferred over virtual getters in hot paths

4. **Configuration System** (`Subanalysis/Properties/histogram_conf/`)
   - Histogram configs in `histogram1D.conf`, `histogram2D.conf` (semicolon-delimited)
   - Parser: `KLOE::Histograms::LoadHistogramConfigs1D/2D()` in `const.cpp`
   - Parser respects quoted fields with escaped commas - use `line.Tokenize(";")` not `,`

### Data Flow
```
TChain.Process(TSelector) 
  → Begin() [initialize histograms, load configs]
  → Process() [loop events, apply cuts, fill histograms]
  → SlaveBegin/Process/SlaveTerminate [parallel execution if enabled]
  → Terminate() [fits, statistics, CSV output, visualization]
```

## Critical Patterns & Gotchas

### Memory Management (ROOT Objects)
- **TF1 cleanup**: Use `constraint->Delete()` not `delete constraint` (ROOT memory model)
- **Destructor pattern**: Call `_constraints[i]->Delete()` in loop, then clear vector
- **Histogram errors**: After `Scale()`, call `hist->Sumw2(kFALSE)` to disable auto-errors
- **Vector assignment**: Validate `vector.size()` before `operator=` to prevent segfaults

### Kinematic Fitting
- **Direct constraints disabled by default**: `_useDirectConstraints = false` (guard against null pointers)
- **Performance**: Numerical gradient calculation is 10-20× slower than analytical - prioritize `FitFunctionParallel` with internal OpenMP
- **Convergence**: Add timeout protection (~1000ms) to catch ill-conditioned matrices
- **Early termination**: Break loop if `chi2 < target_threshold` after N iterations

### Configuration & Parsing
- **Tokenizer with quotes**: Use custom parser respecting quoted fields:
  ```cpp
  TObjArray *tokens = line.Tokenize(";");
  TString value = TString(((TObjString *)tokens->At(i))->String()).Strip(TString::kBoth);
  value.ReplaceAll("\"", "");  // Remove quotes
  ```
- **Config validation**: Always check `configs.find(name) != configs.end()` before access
- **Guard clauses preferred**: Check histogram entries/validity before drawing:
  ```cpp
  if (!hist || hist->GetEntries() == 0) continue;
  ```

### Histogram & Plotting
- **TLatex vs TPaveText**: Use `TPaveText` for multi-line stats (auto-spacing, background box)
- **Normalized Gaussian**: Use `"gausn"` for amplitude = integral parameter, `"gaus"` for amplitude = peak height
- **Fit initialization**: Always set `SetParameter()` before `Fit()`, use `SetParLimits()` for constraints
- **Draw options**: `"HIST SAME"` for MC (lines), `"PE"` for Data (points+errors)
- **Coordinate systems**: `SetNDC()` for normalized (0-1), without for axis units

### CSV Output & Logging
- **Results structure**: Include Signal/Total counts, efficiency, purity, pull fits, chi2/ndf
- **Cut statistics**: Use `CountEventsAroundCut(histogram, cutValue)` → returns events below/above
- **Overflow handling**: Use `Integral(0, nbins+1)` to include underflow/overflow bins
- **File paths**: Always sanitize option strings with `option.ToLower()` for folder names

## Build & Development

### Build System
```bash
cd /data/ssd/gamrat/KLOE/build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$(nproc)
```
- **OpenMP**: Configured in CMakeLists.txt, enabled via `#pragma omp` or `_OPENMP` checks
- **RooFit**: Requires ROOT compiled with RooFit support (check `root-config --features`)

### Debugging Commands
```bash
# Memory leak detection
valgrind --leak-check=full ./bin/KLSPM00

# AddressSanitizer
cmake -DCMAKE_BUILD_TYPE=Debug ..
export ASAN_OPTIONS="abort_on_error=1:detect_leaks=1"

# Segfault with full stack
gdb --args ./bin/KLSPM00
```

### Common Issues & Fixes
| Issue | Root Cause | Solution |
|-------|-----------|----------|
| `Illegal axis coordinates` | Empty histogram | Check `GetEntries() > 0` before drawing |
| `Missing "}"` in TLatex | Unescaped braces in config | Use `"m^{inv}_{2#gamma}"` with quotes |
| Vector assignment crash | Uninitialized member vector | Call `resize(size)` before `operator=` |
| Fit fails to converge | Poor initial conditions | Set `SetParameter(i, estimatedValue)` |

## Key Files to Know
- **Config**: `Subanalysis/Properties/histogram_conf/histogram*.conf`
- **Utilities**: `Include/Codes/src/const.cpp` (histogram creation, cut loading)
- **Analysis**: `Subanalysis/InitialAnalysis/src/initanalysis_full.cpp` (main workflow)
- **Framework**: `Include/Codes/src/KinFitter.cpp` (core fitting engine)

## Useful ROOT/Physics Commands
```cpp
// Validate histogram
if (!hist || hist->GetEntries() == 0) return;

// Safe vector initialization
std::vector<Float_t> data(size, 0.0f);

// PCA analysis (for 2D correlations)
pca->AddRow(data);
pca->MakePrincipals();
Double_t sigma_transverse = TMath::Sqrt((*pca->GetEigenValues())[1]);

// Efficiency calculation
Double_t efficiency = (Double_t)signal_passing / signal_total;
```

## Testing New Features
1. **Always use guard clauses** - validate input before processing
2. **Add debug output** - use `std::cout << "DEBUG: ..."` with condition checks
3. **Test with empty data** - ensure code handles `GetEntries() == 0`
4. **Check memory with valgrind** - before committing multi-histogram changes
5. **Validate against physics** - check counts match expected signal/background ratios
