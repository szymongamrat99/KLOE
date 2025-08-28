#include "KaonMassReconstructor.h"
#include <iostream>
#include <TChain.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TFile.h>

/**
 * @brief Example of how to use KaonMassReconstructor class
 * 
 * This example shows how to integrate the KaonMassReconstructor with existing
 * KLOE analysis code, replacing the manual implementation in kchrec_Kmass.cpp
 */

void kaonReconstructionExample() {
    // Initialize KLOE utilities
    KLOE::pm00 kloeObj;
    
    // Setup input chain (example)
    TChain* chain = new TChain("T");
    // chain->Add("your_data_file.root");
    
    // Physical constants
    const Float_t mPiCh = 0.13957018; // Charged pion mass (GeV)
    
    // Setup branch addresses (similar to kchrec_Kmass.cpp)
    std::vector<Float_t> 
        *KchboostKL = nullptr,
        *ipKL = nullptr,
        *trkKL0 = nullptr,
        *trkKL1 = nullptr;
    
    chain->SetBranchAddress("KchboostKL", &KchboostKL);
    chain->SetBranchAddress("ipKL", &ipKL);
    chain->SetBranchAddress("trk1KL", &trkKL0);
    chain->SetBranchAddress("trk2KL", &trkKL1);
    
    // Create histograms
    TH1F* hKaonMass = new TH1F("hKaonMass", "Reconstructed Kaon Mass;Mass (GeV);Events", 
                               100, 0.4, 0.6);
    TH1F* hMomDiff = new TH1F("hMomDiff", "Momentum Difference;|p_rec - p_orig| (GeV);Events", 
                              100, 0, 0.5);
    
    Int_t nentries = chain->GetEntries();
    
    for(Int_t i = 0; i < nentries; i++) {
        chain->GetEntry(i);
        
        // Skip events with insufficient data
        if(!KchboostKL || !ipKL || !trkKL0 || !trkKL1) continue;
        if(KchboostKL->size() < 9 || ipKL->size() < 3) continue;
        if(trkKL0->size() < 4 || trkKL1->size() < 4) continue;
        
        // Prepare input data for KaonMassReconstructor
        std::vector<Float_t> kaonBoost = *KchboostKL;  // Copy the vector
        std::vector<Float_t> interactionPoint = *ipKL; // Copy the vector
        
        std::vector<std::vector<Float_t>*> tracks = {trkKL0, trkKL1};
        
        try {
            // Use KaonMassReconstructor instead of manual implementation
            KLOE::KaonReconstructionResult result = 
                KLOE::KaonMassReconstructor::reconstructKaonMass(
                    kaonBoost, interactionPoint, tracks, mPiCh, kloeObj
                );
            
            // Fill histograms with results
            hKaonMass->Fill(result.KaonTwoBody[5]); // Reconstructed mass
            
            // Compare reconstructed vs original momentum (if available)
            Float_t pOrig = sqrt(pow(kaonBoost[0], 2) + 
                                pow(kaonBoost[1], 2) + 
                                pow(kaonBoost[2], 2));
            Float_t pRec = result.KaonTwoBody[4];
            hMomDiff->Fill(abs(pRec - pOrig));
            
            // Print some results for debugging
            if(i < 10) {
                std::cout << "Event " << i << ":" << std::endl;
                std::cout << "  Reconstructed Kaon Mass: " << result.KaonTwoBody[5] << " GeV" << std::endl;
                std::cout << "  Original momentum: " << pOrig << " GeV" << std::endl;
                std::cout << "  Reconstructed momentum: " << pRec << " GeV" << std::endl;
                std::cout << "  Pion 1 momentum: (" 
                          << result.track1TwoBody[0] << ", "
                          << result.track1TwoBody[1] << ", "
                          << result.track1TwoBody[2] << ") GeV" << std::endl;
                std::cout << "  Pion 2 momentum: (" 
                          << result.track2TwoBody[0] << ", "
                          << result.track2TwoBody[1] << ", "
                          << result.track2TwoBody[2] << ") GeV" << std::endl;
                std::cout << std::endl;
            }
            
        } catch(const std::exception& e) {
            std::cerr << "Reconstruction failed for event " << i << ": " << e.what() << std::endl;
        }
    }
    
    // Save results
    TFile* outFile = new TFile("kaon_reconstruction_results.root", "RECREATE");
    hKaonMass->Write();
    hMomDiff->Write();
    
    // Create and save plots
    TCanvas* c1 = new TCanvas("c1", "Kaon Mass", 800, 600);
    hKaonMass->Draw();
    c1->SaveAs("kaon_mass_spectrum.png");
    
    TCanvas* c2 = new TCanvas("c2", "Momentum Difference", 800, 600);
    hMomDiff->Draw();
    c2->SaveAs("momentum_difference.png");
    
    outFile->Close();
    
    std::cout << "Analysis completed. Results saved to kaon_reconstruction_results.root" << std::endl;
}

/**
 * @brief How to replace the manual implementation in kchrec_Kmass.cpp
 * 
 * Instead of the manual calculation in kchrec_Kmass.cpp, you can replace
 * the entire reconstruction section with:
 */
void replaceManualImplementation() {
    std::cout << "\n=== How to replace manual implementation ===" << std::endl;
    std::cout << "In kchrec_Kmass.cpp, replace the section from line ~150 to ~300 with:" << std::endl;
    std::cout << R"(
    // Prepare input vectors
    std::vector<Float_t> kaonBoost = *KchboostKL;
    std::vector<Float_t> interactionPoint = *ipKL;
    std::vector<std::vector<Float_t>*> tracks = {trkKL[0], trkKL[1]};
    
    // Use KaonMassReconstructor
    KLOE::KaonReconstructionResult result = 
        KLOE::KaonMassReconstructor::reconstructKaonMass(
            kaonBoost, interactionPoint, tracks, mPiCh, Obj
        );
    
    // Copy results to existing variables
    baseKin.KchrecKLTwoBody = result.KaonTwoBody;
    trkKLTwoBody1 = result.track1TwoBody;
    trkKLTwoBody2 = result.track2TwoBody;
    )" << std::endl;
}

int main() {
    std::cout << "KaonMassReconstructor Example" << std::endl;
    std::cout << "=============================" << std::endl;
    
    replaceManualImplementation();
    
    // Uncomment to run the actual example:
    // kaonReconstructionExample();
    
    return 0;
}
