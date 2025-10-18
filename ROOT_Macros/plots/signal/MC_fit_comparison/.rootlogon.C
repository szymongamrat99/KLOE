{
  cout << "Loading shared library..." << endl
       << endl;
  gSystem->Load("/data/ssd/gamrat/KLOE/build/Include/Codes/libLibRec.so");

  // Global style of histograms, pads, etc.

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

  gStyle->SetHistLineWidth(2);
  gStyle->SetTitle("");

  cout << "Units: cm, ns, tau_{S}, MeV/c^{2}, MeV/c, MeV" << endl;
  TString units[6] = {" [cm]", " [ns]", " [#tau_{S}]", " [MeV/c^{2}]", " [MeV/c]", " [MeV]"};

  const unsigned int KLOE::channNum = 7;
  const Color_t channColor[KLOE::channNum] = {kRed, kGreen, kViolet, kCyan, kBlue, kGreen - 1, kYellow};
  const std::vector<TString> channNames = {"K_{S}K_{L}#rightarrow#pi^{+}#pi^{-}#pi^{0}#pi^{0}", "Regeneration",
                                           "#omega#pi^{0}#rightarrow#pi^{+}#pi^{-}#pi^{0}#pi^{0}",
                                           "K_{S}K_{L}#rightarrow#pi^{+}#pi^{-}3#pi^{0}",
                                           "K_{S}K_{L}#rightarrow#pi^{#pm}l^{#mp}#nu#pi^{0}#pi^{0}",
                                           "Other bcg", "K_{S}K_{L}#rightarrow#pi^{+}#pi^{-}#pi^{+}#pi^{-}"};
  const Color_t dataColor = kBlack;

  const Float_t PhysicsConstants::mK0 = 497.611;   // MeV/c^{2}
  const Float_t PhysicsConstants::mPhi = 1019.455; // MeV/c^{2}
  const Float_t PhysicsConstants::mPi0 = 134.9768; // MeV/c^{2}

  const Float_t tau_S = 0.087;
  const Float_t PhysicsConstants::cVel = 30; // cm/ns
}
