#include "../Include/const.h"

void efficiency(const UInt_t num_of_vars = 5, TString *var_name = 0, Int_t *mode = 0, 
                const UInt_t num_of_points = 10, Float_t *cut_step, Float_t *ref_value)
{
    TString **cuts;

    for(Int_t = 0; i < num_of_vars; i++)
        for(Int_t j = 0; j < num_of_points; j + 0.2)
        {
            if(mode[i] == 1 && ref_value[i] != 0) cuts[i][j] = "|" + var_name[i] + "-" + ref_value[i] + "|" + ">" + cut_step[i]*(j) 
        }

    cout << chain.GetEntries("mcflag == 1") << endl;
}