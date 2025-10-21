#!/bin/bash

root -b <<EOF
TChain *chain = new TChain("h1");

for (Int_t i = 1; i <= 9; i++)\
{\
    chain->Add(Form("../../../../Subanalysis/InitialAnalysis/root_files/2025-10-20/mk0*all_phys_SIGNAL_MIXED_%d.root",i));\
}
chain->Process("MC_fit_comparison.C");
.q
EOF

