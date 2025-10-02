#!/bin/bash

root -b <<EOF
TChain *chain = new TChain("h1");

for (Int_t i = 1; i <= 1; i++)\
{\
    chain->Add(Form("../../../../InitialAnalysis/root_files/2025-10-02/mk0*all_phys_SIGNAL_MIXED_%d.root",i));\
}
chain->Process("MC_fit_comparison.C");
.q
EOF

