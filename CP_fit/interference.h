#ifndef INTERFERENCE_H
#define INTERFERENCE_H

#include "kloe_class.h"
#include "../../Include/const.h" 
#include "TMath.h"
#include "TString.h"
#include "TH1.h"
#include <string>
#include <iostream>   

namespace KLOE
{
    class interference: public pm00
    {
        public:
            UInt_t bin_number = 201;
            Double_t x_min = -100., x_max = 100., y_min = 0., y_max = 1000.;

            UInt_t ev_num[chann_num];

            Double_t *efficiency, *correction;

            interference(TString mode_m)
            {
                mode = mode_m;
            };

            void bin_extraction(UInt_t channel, TH1* histogram);

            Double_t interf_function(const Float_t x, Int_t check = 1, const Double_t *par = 0);

            Double_t interf_chi2_split(const Double_t *xx);
            Double_t interf_chi2_window(const Double_t *xx);
            Double_t interf_chi2_excluded(const Double_t *xx);
            Double_t interf_chi2_all(const Double_t *xx);

            Double_t interf_chi2_mc(const Double_t *xx);
            Double_t interf_chi2_bcg(const Double_t *xx);

            Double_t interf_chi2(const Double_t *xx);

        private:
            TString mode = "split";  // "split", "window", "excluded, "mc", "bcg", "all"

            TH1 *frac[7];
            Double_t *b[7], *e[7]; 

            void histos_init()
            {
                for(Int_t i = 0; i < 7; i++)
                    frac[i] = new TH1F(("histo" + std::to_string(i)).c_str(), "", bin_number, x_min, x_max);
            };

    };
}

#endif