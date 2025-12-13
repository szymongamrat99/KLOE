#define SignalCutOptimizer_cxx
// The class definition in SignalCutOptimizer.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("SignalCutOptimizer.C")
// root> T->Process("SignalCutOptimizer.C","some options")
// root> T->Process("SignalCutOptimizer.C+")
//

#include "SignalCutOptimizer.h"
#include <TH2.h>
#include <TStyle.h>
#include <const.h>

void SignalCutOptimizer::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();

  fCutsOneSided.clear();

  fCutsOneSided["Chi2SignalKinFit"] = [this](Double_t cutValue) {
    // Example cut: VarX < P
    Double_t varXValue = 0.0;

    if (fVarXName == "Chi2SignalKinFit")
      varXValue = *Chi2SignalKinFit / 10.0; // Example scaling

    return varXValue < cutValue;
  };
}

void SignalCutOptimizer::SlaveBegin(TTree *tree)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();

  fChain = tree;

  TString maxKey;
  Int_t maxValue = 0;

  for (const auto &lumi : channLumi)
  {
    if (lumi.second > maxValue)
    {
      maxValue = lumi.second;
      maxKey = lumi.first;
    }
  }

  for (const auto &lumi : channLumi)
  {
    channFactor[lumi.first] = channLumi.at(maxKey) / lumi.second;

    channEventsCut[lumi.first] = 0;
  }
}

Bool_t SignalCutOptimizer::Process(Long64_t entry)
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // When processing keyed objects with PROOF, the object is already loaded
  // and is available via the fObject pointer.
  //
  // This function should contain the \"body\" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.
  //
  // The processing can be stopped by calling Abort().
  //
  // Use fStatus to set the return value of TTree::Process().
  //
  // The return value is currently not used.

  TString fileNameTmp;

  fReader.SetLocalEntry(entry);

  fmctruth = *mctruth;

  fileNameTmp = (TString)fChain->GetCurrentFile()->GetName();

  if (fmctruth == 0 && *mcflag == 1)
  {
    for (auto const &category : categoryMap)
    {
      if (fileNameTmp.Contains(category.first))
      {
        fmctruth = category.second;
        break;
      }
    }
  }

  if (fmctruth != -1)
  {

  Double_t combinedMassPi0Err = sqrt(pow(pi01Fit[5] - PhysicsConstants::mPi0, 2) + 
                                     pow(pi02Fit[5] - PhysicsConstants::mPi0, 2) );

  fCutsOneSided["Chi2SignalKinFit"](fParamP);

  if (abs(Kchrec[5] - PhysicsConstants::mK0) > 3)
    return kFALSE;

  if (abs(*minv4gam - PhysicsConstants::mK0) > 50)
    return kFALSE;

  if (combinedMassPi0Err > 15)
    return kFALSE;

  if (*Chi2SignalKinFit / 10. > fParamP)
    return kFALSE;

  channEventsCut[KLOE::channName.at(fmctruth)]++;

  }

  return kTRUE;
}

void SignalCutOptimizer::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.
}

void SignalCutOptimizer::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  std::map<TString, Float_t> channEventsCutCorr, channEventsTotalCorr;

  for (const auto &channelType : KLOE::channName)
  {
    channEventsCutCorr[channelType.second] = channEventsCut[channelType.second] * channFactor[channelType.second];

    fNCut += channEventsCutCorr[channelType.second];
  }
}

#ifdef SignalCutOptimizer_cxx
void SignalCutOptimizer::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the reader is initialized.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  fReader.SetTree(tree);
}

Bool_t SignalCutOptimizer::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

#endif // #ifdef SignalCutOptimizer_cxx