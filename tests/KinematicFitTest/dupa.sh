#!/bin/bash

g++ pi0pi0_kinfit.cpp `root-config --libs --cflags` -lLibRec -lboost_filesystem -lboost_system -lcurl -Wl,-rpath,.

./a.out
