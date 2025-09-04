#!/bin/bash

root -b <<EOF
TChain *chain = new TChain("INTERF/h1");
chain->Add("/data/4/users/gamrat/old_root_files/MONTE_CARLO/*.root");
chain->Add("/data/4/users/gamrat/old_root_files/DATA/*.root");
chain->AddFriend("h_gen_vars", "/data/4/users/gamrat/scripts/Scripts/Subanalysis/GeneratedVars/root_files/mctruth_1_56.root");
chain->Process("conferenceCuts.C");
.q
EOF

