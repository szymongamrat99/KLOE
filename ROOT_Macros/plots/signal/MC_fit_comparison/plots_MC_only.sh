#!/bin/bash

root -b <<EOF
TChain *chain = new TChain("h1");

for (Int_t i = 1; i <= 73; i++)\
{\
    chain->Add(Form("../../../../Subanalysis/InitialAnalysis/root_files/2025-11-10/mk0*all_phys_SIGNAL_MIXED_Signal_%d.root",i));\
}
chain->Process("MC_fit_comparison.C", "Chi2SignalReduced+DeltaPhivFit");
.q
EOF

