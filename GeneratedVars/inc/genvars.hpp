#ifndef GENVARS_H
#define GENVARS_H

#include <iostream>
#include <fstream>
#include <vector>

#include <TTree.h>
#include <TFile.h>

#include "../../../Include/const.h"
#include "../../../Include/Codes/chain_init.cpp"
#include "../../../Include/Codes/ErrorLogs.h"
#include "../../../Include/Codes/MainMenu.h"

int genvars(int, int, int, Controls::DataType);

#endif //! GENVARS_H