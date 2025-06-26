#include "arrays_equality.h"

#include <iostream>
#include<array> 

//! These arrays should have the same length to be equal

Int_t arr_eq(Int_t *arr1, Int_t *arr2)
{
  Int_t bad_clus = 0, bad_clus_all = 0;

  const Int_t length1 = 4;
  const Int_t length2 = 4;

  Int_t sorted_ind1[length1], sorted_ind2[length2];

  if( length1 != length2 )
    return kFALSE;

  std::sort(arr1, arr1 + length1);
  std::sort(arr2, arr2 + length2);

  for(Int_t i = 0; i < length1; i++)
  {
    bad_clus = 0;
    for(Int_t j = 0; j < length2; j++)
    {
      if(arr1[i] != arr2[j] && bad_clus == 0)
      {
        bad_clus++;
      }

      if(arr1[i] == arr2[j])
      {
        bad_clus = 0;
        break;
      }
    }

    bad_clus_all += bad_clus;

  }
  
  return bad_clus_all;
}