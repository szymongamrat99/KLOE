#include <fstream>
#include <iostream>

#include <TVector3.h>
#include <TLorentzVector.h>

#include "neutral_mom.h"
#include "ErrorLogs.h"
#include <const.h>

int neu_triangle(Double_t *, Double_t *, Double_t Clu5Vec[4][5], Double_t *ip, Double_t *Phi4Mom, Double_t *Kne4Mom, Double_t *Kne4Vec, Double_t *trc, ErrorHandling::ErrorLogs &logger);
