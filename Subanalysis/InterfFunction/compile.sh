#!/bin/bash

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:../../build/Include/Codes

g++ back_to_the_future_histograms.cpp -g -o bttf_histo `root-config --cflags --glibs` -Wall -Wextra -std=c++14 -L../../build/Include/Codes -lLibRec -lcurl -lboost_filesystem -lboost_system

g++ back_to_the_future_sampling.cpp -g -o bttf_sample `root-config --cflags --glibs` -Wall -Wextra -std=c++14 -L../../build/Include/Codes -lLibRec -lcurl -lboost_filesystem -lboost_system

g++ interf_func_bttf_draw.cpp -g -o bttf_draw `root-config --cflags --glibs` -Wall -Wextra -std=c++14 -L../../build/Include/Codes -lLibRec -lcurl -lboost_filesystem -lboost_system

g++ interf_func_draw.cpp -g -o func_draw `root-config --cflags --glibs` -Wall -Wextra -std=c++14 -L../../build/Include/Codes -lLibRec -lcurl -lboost_filesystem -lboost_system

g++ interf_func_compare.cpp -g -o func_compare `root-config --cflags --glibs` -Wall -Wextra -std=c++14 -L../../build/Include/Codes -lLibRec -lcurl -lboost_filesystem -lboost_system