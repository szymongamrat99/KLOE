#!/bin/bash

root -b <<EOF
TChain *chain = new TChain("h1");

for (Int_t i = 1; i <= 7; i++)\
{\
    chain->Add(Form("../../../../InitialAnalysis/root_files/2025-10-16/mk0*all_phys2_SIGNAL_MIXED_%d.root",i));\
}\
for (Int_t i = 1; i <= 5; i++)\
{\
    chain->Add(Form("../../../../InitialAnalysis/root_files/2025-10-16/dk0*SIGNAL_MIXED_%d.root",i));\
}\
chain->Process("MC_fit_comparison.C");
.q
EOF

