{  
  cout << "Loading shared library..." << endl << endl;
  gSystem->Load("/internal/big_one/4/users/gamrat/scripts/Include/librec.so");

  //Global style of histograms, pads, etc.

  gStyle->SetOptStat("iouMn");
  gStyle->SetStatFormat("6.2f");

  gStyle->SetCanvasDefH(750);
  gStyle->SetCanvasDefW(750);

  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  gStyle->SetOptLogz(1);
  gStyle->SetPalette(1);

  gStyle->SetLabelSize(0.04, "X");
  gStyle->SetLabelSize(0.04, "Y");
  gStyle->SetLabelSize(0.04, "Z");
  gStyle->SetLabelOffset(0.02, "X");
  gStyle->SetLabelOffset(0.02, "Y");
  gStyle->SetLabelOffset(0.02, "Z");
  gStyle->SetLabelFont(62, "X");
  gStyle->SetLabelFont(62, "Y");
  gStyle->SetLabelFont(62, "Z");

  gStyle->SetTitleSize(0.05, "X");
  gStyle->SetTitleSize(0.05, "Y");
  gStyle->SetTitleSize(0.05, "Z");
  gStyle->SetTitleOffset(1.2, "X");
  gStyle->SetTitleOffset(1.5, "Y");
  gStyle->SetTitleOffset(1.2, "Z");
  gStyle->SetTitleFont(62, "X");
  gStyle->SetTitleFont(62, "Y");
  gStyle->SetTitleFont(62, "Z");

  gStyle->SetHistLineWidth(3);
  gStyle->SetTitle("");

  cout << "Units: cm, ns, tau_{S}, MeV/c^{2}, MeV/c, MeV" << endl;
  TString units[6] = {" [cm]", " [ns]", " [#tau_{S}]", " [MeV/c^{2}]", " [MeV/c]", " [MeV]"};
}
