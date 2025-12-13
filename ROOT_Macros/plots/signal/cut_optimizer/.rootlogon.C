{
  cout << "Loading shared library..." << endl
       << endl;
  gSystem->Load("/data/ssd/gamrat/KLOE/build/Include/Codes/libLibRec.so");
    // Set of the external libraries
  gSystem->Load("libcurl");

  // Creation of the build dir for selector compilation
  const char *outputDir = "build";
  gSystem->mkdir(outputDir, kTRUE);
  gSystem->SetBuildDir(outputDir, kTRUE);
  gROOT->LoadMacro("src/SignalCutOptimizer.C++");
  gSystem->Load("/data/ssd/gamrat/KLOE/build/Include/Codes/SignalCutOptimizer_C.so");
  cout << "Loaded LibRec.so. Project is finally consistent ;-)" << endl
       << endl;
}
