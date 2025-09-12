#!/bin/bash

root -b <<EOF
TChain *chain = new TChain("h1");
for (Int_t i = 1; i <= 5; i++)\
{\
    chain->Add(Form("../../InitialAnalysis/root_files/2025-09-11/dk0*_SIGNAL_MIXED_%d.root",i));\
}

for (Int_t i = 1; i <= 7; i++)\
{\
    chain->Add(Form("../../InitialAnalysis/root_files/2025-09-11/mk0*_SIGNAL_MIXED_%d.root",i));\
}
chain->Process("init_analysis.C");
.q
EOF

