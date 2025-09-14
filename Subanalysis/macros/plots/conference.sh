#!/bin/bash

root -b <<EOF
TChain *chain = new TChain("h1");

for (Int_t i = 1; i <= 3; i++)\
{\
    chain->Add(Form("../../InitialAnalysis/root_files/2025-09-14/mk0*all_phys2_SIGNAL_MIXED_%d.root",i));\
    chain->Add(Form("../../InitialAnalysis/root_files/2025-09-14/mk0*all_phys3_SIGNAL_MIXED_%d.root",i));\
}
chain->Process("init_analysis.C");
.q
EOF

