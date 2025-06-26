#include <fstream>
#include <iostream>

#include <TVector3.h>
#include <TLorentzVector.h>

#include "neutral_mom.h"
#include "ErrorLogs.h"
#include <const.h>

int neu_triangle(Float_t *, Float_t *, Float_t Clu5Vec[4][5], Float_t *ip, Float_t *Phi4Mom, Float_t *Kne4Mom, Float_t *Kne4Vec, Float_t *trc, ErrorHandling::ErrorLogs &logger);
