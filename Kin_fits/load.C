{
    TString fullname = "",
            dirnamemc = "mc_root_files_vtx_part", dirnamedata = "data_root_files_vtx_part",
            filenamemc = "mc_stream62_mccard2_", filenamedata = "data_stream42_",
            extension = ".root";

    TChain chain("INTERF/h1");

    for(Int_t i = 1; i <= 1; i++)
    {
        if(i != 10)
        {
            fullname = "../../ROOT_files/" + dirnamemc + "/" + filenamemc + i + extension;
            chain.Add(fullname);

            fullname = "../../ROOT_files/" + dirnamedata + "/" + filenamedata + i + extension;
            chain.Add(fullname);
        }
    }
}