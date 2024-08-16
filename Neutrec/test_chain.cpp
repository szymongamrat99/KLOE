#include "../../Include/Codes/chain_init.h"

UInt_t nentriesret(UInt_t first, UInt_t last)
{
  TChain *chain = new TChain("INTERF/h1");
	chain_init(chain, first, last);

  return (UInt_t)chain->GetEntries();
}