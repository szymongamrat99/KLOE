{
    TString fullname = "",
            dirnamemc = "230623_mc", dirnamedata = "230623_data",
            filenamemc = "mc_stream62_mccard2_", filenamedata = "data_stream42_",
            extension = ".root";

    TChain chain("INTERF/h1");

    for(Int_t i = 1; i <= 583; i++)
    {
        if(i != 64 && i != 100)
        {
            fullname = "../ROOT_files/" + dirnamemc + "/backup/" + filenamemc + i + extension;
            chain.Add(fullname);
        }
    }
    /*for(Int_t i = 1; i <= 10; i++)
    {
    	fullname = "../ROOT_files/" + dirnamedata + "/backup/" + filenamedata + i + extension;
    	chain.Add(fullname);
    }*/
}
