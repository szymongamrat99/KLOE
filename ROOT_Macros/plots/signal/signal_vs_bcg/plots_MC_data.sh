#!/bin/bash

root -b <<EOF
TChain *chain = new TChain("h1");

for (Int_t i = 1; i <= 1; i++)\
{\
    chain->Add(Form("../../../../Subanalysis/InitialAnalysis/root_files/2025-10-20/mk0*all_phys2_SIGNAL_MIXED_%d.root",i));\
}\
for (Int_t i = 1; i <= 1; i++)\
{\
    chain->Add(Form("../../../../Subanalysis/InitialAnalysis/root_files/2025-10-20/dk0*SIGNAL_MIXED_%d.root",i));\
}\
chain->Process("signal_vs_bcg_v2.C");
.q
EOF

