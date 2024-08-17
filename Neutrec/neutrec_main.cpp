#include <iostream>
#include <chrono>

#include <TMath.h>

#include "inc/trilateration.hpp"

using namespace std;

namespace Neutrec
{
  int main(int argc, char *argv[])
  {
    using std::chrono::duration;
    using std::chrono::duration_cast;
    using std::chrono::high_resolution_clock;
    using std::chrono::minutes;

    Int_t first = atoi(argv[1]), last = atoi(argv[2]);
    Short_t good_clus = atoi(argv[3]), loopcount = atoi(argv[4]), jmin = atoi(argv[5]), jmax = atoi(argv[6]), ind_data_mc = atoi(argv[7]);

    cout << "Choose the reconstruction method: " << endl;
    cout << "1. Bare trilateration." << endl;
    cout << "2. Trilateration with kinematic fit." << endl;
    cout << "3. Trilateration with kinematic fit and T0 correction." << endl;
    cout << "4. Trilateration with triangle." << endl;

    Int_t mode;

    cin >> mode;

    auto t1 = high_resolution_clock::now();

    switch (mode)
    {

    case 1:
      cout << "Analysis started." << endl;
      tri_neurec(first, last, good_clus);
      break;

    case 2:
      cout << "Analysis started." << endl;
      tri_neurec_kinfit_corr(ind_data_mc, first, last, loopcount, jmin, jmax);
      break;

    case 4:
      cout << "Analysis started." << endl;
      triangle_neurec(first, last, good_clus);
      break;

    default:
      cout << "No option chosen, exiting..." << endl;
      break;
    }

    auto t2 = high_resolution_clock::now();

    /* Getting number of minutes as an integer. */
    auto ms_int = duration_cast<minutes>(t2 - t1);

    /* Getting number of minutes as a double. */
    duration<double, std::milli> ms_double = t2 - t1;

    std::cout << ms_int.count() << " minutes\n";

    return 0;
  }
}