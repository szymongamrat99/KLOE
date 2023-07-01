#!/bin/bash

# Neutral vertex
cd Neutrec/

root <<EOF
.L six_gamma_selection.cpp
for(Int_t zz = 1; zz <= 99; zz++){ six_gamma_selection(zz,"230623_data","data_stream42_"); cout << "File done: " << zz << endl; };
.q
EOF

root <<EOF
.L four_gamma_selection.cpp
for(Int_t zz = 1; zz <= 99; zz++){ four_gamma_selection(zz,"230623_data","data_stream42_"); cout << "File done: " << zz << endl; };
.q
EOF

# Semileptonic
cd ../Semileptonic/

root <<EOF
.L semi_selection.cpp
for(Int_t zz = 1; zz <= 99; zz++){ semi_selection(zz,"230623_data","data_stream42_"); cout << "File done: " << zz << endl; };
.q
EOF

# Double pipi
cd ../Pipi/

root <<EOF
.L doublepipi_selection.cpp
for(Int_t zz = 1; zz <= 99; zz++){ doublepipi_selection(zz,"230623_data","data_stream42_"); cout << "File done: " << zz << endl; };
.q
EOF

# Generated variables
cd ../Generated_vars/

root <<EOF
.L generated_variables.cpp
for(Int_t zz = 1; zz <= 99; zz++){ generated_variables(zz,"230623_data","data_stream42_"); cout << "File done: " << zz << endl; };
.q
EOF

root <<EOF
.L pipi_find.cpp
for(Int_t zz = 1; zz <= 99; zz++){ pipi_find(zz,"230623_data","data_stream42_"); cout << "File done: " << zz << endl; };
.q
EOF

cd ../
