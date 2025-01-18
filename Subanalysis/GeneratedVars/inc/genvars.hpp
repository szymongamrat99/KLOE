#ifndef GENVARS_H
#define GENVARS_H

#include <iostream>
#include <fstream>
#include <vector>
#include <omp.h>

#include <TTree.h>
#include <TFile.h>

#include <const.h>
#include <ErrorLogs.h>
#include <MainMenu.h>
#include <fort_common.h>

int genvars(int, int, int);
Int_t split_channels(Int_t, Int_t);

#endif //! GENVARS_H