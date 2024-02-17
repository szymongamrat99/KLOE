#include <iostream>

#include <TMath.h>

#include "inc/trilateration.hpp"

using namespace std;

int main(int argc, char *argv[])
{
  Int_t first = atoi(argv[1]), last = atoi(argv[2]), good_clus = atoi(argv[3]);

  cout << "Choose the reconstruction method: " << endl;
  cout << "1. Bare trilateration." << endl;
  cout << "2. Trilateration with kinematic fit." << endl;
  cout << "3. Trilateration with kinematic fit and T0 correction." << endl;

  Int_t mode;

  cin >> mode;

  switch (mode)
  {

    case 1:
      cout << "Analysis started." << endl;
      tri_neurec(first, last, good_clus);
      break;

    case 2:
      cout << "Analysis started." << endl;
      tri_neurec_kinfit_corr(first, last);
      break;

    default:
      cout << "No option chosen, exiting..." << endl;
      break;
  }

  return 0;
}