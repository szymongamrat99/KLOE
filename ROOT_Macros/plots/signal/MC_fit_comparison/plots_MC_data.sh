#!/bin/bash

root -b <<EOF
TChain *chain = new TChain("h1");

for (Int_t i = 1; i <= 1; i++)\
{\
    chain->Add(Form("../../../../InitialAnalysis/root_files/2025-10-19/mk0*all_phys3_SIGNAL_MIXED_%d.root",i));\
}\
for (Int_t i = 1; i <= 1; i++)\
{\
    chain->Add(Form("../../../../InitialAnalysis/root_files/2025-10-19/dk0*SIGNAL_MIXED_%d.root",i));\
}\
chain->Process("dupa.C");
.q
EOF

