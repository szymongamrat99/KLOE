// Author: Szymon Gamrat
// Date of last update: 26.05.2024

#include <iostream>

#include <fort_common.h> // Linking of FORTRAN common block interfcommon
#include <fort_func.h>   // Linking of FORTRAN klspm00_lib library

#include <const.h>       // Linking of const.h with standard definitions for klspm00 analysis

using namespace std;

int main()
{
  TChain *chain = new TChain("h1");
  chain->Add(prod2root_path_v26 + "/*"+ Paths::ext_root);

  chain->SetBranchAddress("PxTv",interfcommon_.PxTV);

  Int_t nentries = chain->GetEntries();


  for(Int_t i = 0; i < nentries; i++)
  {
    chain->GetEntry(i);

    std::cout << interfcommon_.PxTV[0] << std::endl;
  }

  clearstruct_();

  double array[2][4] = {{0., 0., 0., 0.},
                        {0., 0., 0., 0.}};
  double array1[4] = {0., 0., 0., 0.};
  float numbers = 0;

  test_array_(*array, &numbers, array1);

  cout << array[0][0] << " " << array[0][1] << " " << array[0][2] << " " << array[0][3] << endl;
  cout << array[1][0] << " " << array[1][1] << " " << array[1][2] << " " << array[1][3] << endl;

  cout << array1[0] << " " << array1[1] << " " << array1[2] << " " << array1[3] << endl;

  cout << numbers << endl;

  cout << "Test interf: " << interfcommon_.a1typ << endl;

  return 0;
}