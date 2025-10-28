#!/bin/bash

option=$1

root -b <<EOF
TChain *chain = new TChain("h1");
for (Int_t i = 1; i <= 11; i++)\
{\
    chain->Add(Form("../../../../Subanalysis/InitialAnalysis/root_files/2025-10-25/mk0*all_phys2_SIGNAL_MIXED_%d.root",i));\
}\
for (Int_t i = 1; i <= 7; i++)\
{\
    chain->Add(Form("../../../../Subanalysis/InitialAnalysis/root_files/2025-10-25/dk0*SIGNAL_MIXED_%d.root",i));\
}\
chain->Process("signal_vs_bcg_v2.C", "$option");
.q
EOF

