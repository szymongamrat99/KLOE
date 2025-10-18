#ifndef COVARIANCEMATRIX_H
#define COVARIANCEMATRIX_H

#include <iostream>
#include <fstream>
#include <vector>

#include <TVectorT.h>
#include <TMatrixT.h>
#include <TChain.h>

#include <const.h>
#include <ErrorLogs.h>
#include <Logs.h>
#include <MainMenu.h>
#include <kloe_class.h>
#include <MeasQualityGraph.h>
#include <MomentumSmearing.h>
#include <CustomGraph.h>
#include <interference.h>
#include <fort_common.h>

/**
* @file covmatrix.hpp
* @brief Header file for Covariance Matrix determination
* @param chain Address to the TChain used to get the data
* @param data_type Parameter, which type of data should be taken into account: 
                  * - 1: MC signal only 
                  * - 2: MC signal + background 
                  * - 3: MC + Data 
                  * - 4: Entire MC signal (before cuts and with all events with errors)
* @param logger Address to the logger for error and info.
* @param Obj Address to the object from KLOE::pm00 class
* @return Integer flag: 
          * - 0: no error
          * - 1: error
*/

int CovarianceMatrixDetermination(TChain &chain, Controls::DataType &data_type, ErrorHandling::ErrorLogs &logger, KLOE::pm00 &Obj);
int CovarianceMatrixDeterminationControlSample(TChain &chain, Controls::DataType &data_type, ErrorHandling::ErrorLogs &logger, KLOE::pm00 &Obj);

#endif //! COVARIANCEMATRIX_H