#!/bin/bash

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:../../build/Include/Codes

g++ back_to_the_future_histograms.cpp -g -o bttf_histo `root-config --cflags --glibs` -Wextra -std=c++14 -L../../build/Include/Codes -lLibRec -lcurl -lboost_filesystem -lboost_system

g++ back_to_the_future_sampling.cpp -g -o bttf_sample `root-config --cflags --glibs` -Wextra -std=c++14 -L../../build/Include/Codes -lLibRec -lcurl -lboost_filesystem -lboost_system

g++ interf_func_bttf_draw.cpp -g -o bttf_draw `root-config --cflags --glibs` -Wextra -std=c++14 -L../../build/Include/Codes -lLibRec -lcurl -lboost_filesystem -lboost_system

g++ interf_func_draw.cpp -g -o func_draw `root-config --cflags --glibs` -Wextra -std=c++14 -L../../build/Include/Codes -lLibRec -lcurl -lboost_filesystem -lboost_system

g++ src/bttf_sample_hist.cpp -g -o hist_sample_draw `root-config --cflags --glibs` -Wextra -std=c++14 -L../../build/Include/Codes -lLibRec -lcurl -lboost_filesystem -lboost_system

g++ src/bttf_fits_to_histos.cpp -g -o hist_fits `root-config --cflags --glibs` -Wextra -std=c++14 -L../../build/Include/Codes -lLibRec -lcurl -lboost_filesystem -lboost_system

g++ src/bttf_fits_to_histos_integral_range.cpp -g -o hist_fits_integral_range `root-config --cflags --glibs` -Wextra -std=c++14 -L../../build/Include/Codes -lLibRec -lcurl -lboost_filesystem -lboost_system