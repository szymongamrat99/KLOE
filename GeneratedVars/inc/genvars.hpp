#ifndef GENVARS_H
#define GENVARS_H

#include <iostream>
#include <fstream>
#include <vector>

#include <TTree.h>
#include <TFile.h>

#include "../../../Include/const.h"
#include "../../../Include/Codes/ErrorLogs.h"
#include "../../../Include/Codes/MainMenu.h"

int genvars(int, int, int);
Int_t split_channels(Int_t, Int_t);

#endif //! GENVARS_H