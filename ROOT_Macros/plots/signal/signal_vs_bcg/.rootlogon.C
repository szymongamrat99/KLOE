{
  cout << "Loading shared library..." << endl
       << endl;
  gSystem->Load("/data/ssd/gamrat/KLOE/build/Include/Codes/libLibRec.so");
  gSystem->Load("libcurl");
  cout << "Loaded LibRec.so. Project is finally consistent ;-)" << endl
       << endl;
}
