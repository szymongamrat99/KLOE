#ifndef OMEGAREC_HPP
#define OMEGAREC_HPP

#include <iostream>
#include <fstream>
#include <vector>

#include <TTree.h>
#include <TFile.h>

#include <const.h>
#include <ErrorLogs.h>
#include <MainMenu.h>
#include <clear_variables.h>

int omegarec(Int_t, Int_t, Short_t, Short_t, Short_t, Controls::DataType);
int plots(int, int, int, int, int, Controls::DataType);

#endif //! OMEGAREC_HPP