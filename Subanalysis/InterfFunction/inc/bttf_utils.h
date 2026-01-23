#pragma once

#include <TChain.h>

#include <const.h>

// ======== Configuration struct for Back to the Future analysis ========
struct Config
{
  Double_t  tMinX = 0.0,
            tMinY = 0.0,
            tMaxX = 20.0,
            tMaxY = 300.0,
            resT1 = 1.0,
            resT2 = 1.0;

  Int_t binsX = (tMaxX - tMinX) / resT1,
        binsY = (tMaxY - tMinY) / resT2;

  Double_t
          physicsRe = PhysicsConstants::Re,
          physicsIm = PhysicsConstants::Im_nonCPT;
};

// ======== Function prototypes ========
TChain *LoadPM00Chain();
TChain *LoadPMPMChain();


