#!/bin/bash

g++ back_to_the_future_histograms.cpp -o bttf_histo `root-config --cflags --glibs` -Wall -Wextra -std=c++14 -L../../build/Include/Codes -lLibRec -lcurl -lboost_filesystem -lboost_system

g++ back_to_the_future_sampling.cpp -o bttf_sample `root-config --cflags --glibs` -Wall -Wextra -std=c++14 -L../../build/Include/Codes -lLibRec -lcurl -lboost_filesystem -lboost_system

g++ interf_func_bttf_draw.cpp -o bttf_draw `root-config --cflags --glibs` -Wall -Wextra -std=c++14 -L../../build/Include/Codes -lLibRec -lcurl -lboost_filesystem -lboost_system

g++ interf_func_draw.cpp -o func_draw `root-config --cflags --glibs` -Wall -Wextra -std=c++14 -L../../build/Include/Codes -lLibRec -lcurl -lboost_filesystem -lboost_system

g++ interf_func_compare.cpp -o func_compare `root-config --cflags --glibs` -Wall -Wextra -std=c++14 -L../../build/Include/Codes -lLibRec -lcurl -lboost_filesystem -lboost_system