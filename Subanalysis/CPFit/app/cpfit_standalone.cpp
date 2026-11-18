#include <exception>
#include <iostream>

#include <TChain.h>
#include <TROOT.h>

#include "../inc/cpfit.hpp"

int main(int argc, char **argv)
{
  if (argc < 5)
  {
    std::cerr << "Usage: cpfit_standalone CONFIG.json OUTPUT.json TREE FILE.root [FILE.root ...]\n";
    return 2;
  }

  try
  {
    gROOT->SetBatch(kTRUE);
    TChain chain(argv[3]);
    for (int i = 4; i < argc; ++i)
    {
      if (chain.Add(argv[i], 0) == 0)
        throw std::runtime_error(std::string("Cannot add input ROOT file: ") + argv[i]);
    }
    if (chain.GetEntries() <= 0)
      throw std::runtime_error("Input chain has no entries");
    return cp_fit_run(chain, argv[1], argv[2]);
  }
  catch (const std::exception &error)
  {
    std::cerr << "CP fit failed: " << error.what() << '\n';
    return 1;
  }
}
