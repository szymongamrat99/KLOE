{
    TString fullname = "",
            dirnamemc = "MONTE_CARLO", dirnamedata = "DATA",
            filenamemc = "mc_stream62_mccard2_", filenamedata = "data_stream42_",
            extension = ".root";

    TChain chain("INTERF/h1");

    for(Int_t i = 1; i <= 56; i++)
    {
    	fullname = "/internal/big_one/4/users/gamrat/old_root_files/" + dirnamedata + "/" + filenamedata + i + extension;
    	chain.Add(fullname);
    }
    for(Int_t i = 1; i <= 56; i++)
    {
        fullname = "/internal/big_one/4/users/gamrat/old_root_files/" + dirnamemc + "/" + filenamemc + i + extension;
    	chain.Add(fullname);
    }
}
