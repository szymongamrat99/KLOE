#include <iostream>
#include <fstream>

#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

#include <event_data.h>

#include "../inc/initialanalysis.hpp"

int InitialAnalysis_full(TChain &chain, ErrorHandling::ErrorLogs &logger, KLOE::pm00 &Obj)
{
	TTreeReader reader(&chain);
	GeneralEventPropertiesMC eventProps(reader);

	while(reader.Next())
	{
		// Here you would process each entry in the tree.
		// For example, you can read values from the tree and perform calculations.
		// This is a placeholder for your actual analysis logic.

		std::cout << "Processing event with ntmc: " << *eventProps.ntmc << ", nvtxmc: " << *eventProps.nvtxmc << std::endl;
	}

	return 0;
}