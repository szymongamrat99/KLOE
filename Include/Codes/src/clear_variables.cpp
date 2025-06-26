#include <clear_variables.h>

void Clear1DArray(UInt_t M, Float_t *array)
{
  for(UInt_t i = 0; i < M; i++)
    array[i] = 999.;
};

void Clear2DArray(UInt_t M, UInt_t N, Float_t **array)
{
  for(UInt_t i = 0; i < M; i++)
    for(UInt_t j = 0; j < N; j++)
      array[i][j] = 999.;
};