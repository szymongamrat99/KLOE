#include <iostream>
#include <TSystem.h>
#include "kloe_class.h"

using namespace KLOE;

int test()
{
    pm00 first_event;

    first_event.inv_mass_calc(first_event.phi_mom);

    std::cout << first_event.inv_mass << " " << "dupa" << std::endl;

    return 0;
}