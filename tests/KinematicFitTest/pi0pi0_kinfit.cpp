#include <iostream>

#include <TRandom.h>
#include <TMath.h>
#include <TCanvas.h>

#include <ErrorLogs.h>
#include <KinFitter.h>
#include <kloe_class.h>
#include <const.h>

int main()
{
    gErrorIgnoreLevel = 6001;

    Int_t N_free = 2,
          N_const = 1,
          M = 1,
          loopcount = 30;

    TVectorD
        _X_min,
        _X_init_min;

    TMatrixD
        _V_min,
        _V_init;

    Double_t
        _CHISQRMIN;

    _V_min.ResizeTo(N_free + N_const, N_free + N_const);
    _V_init.ResizeTo(N_free + N_const, N_free + N_const);

    _X_min.ResizeTo(N_free + N_const);
    _X_init_min.ResizeTo(N_free + N_const);

    ErrorHandling::ErrorLogs
        logger("dupsko.log");
    KLOE::KinFitter kinFit("Test", N_free, N_const, M, 0, loopcount, 0.01, logger);

    kinFit.ConstraintSet({"MinvConsvNeutralKaon"});

    std::map<TString, TCanvas *> canvasList;
    std::map<TString, TH1 *> histInitList, histAfterList;
    std::vector<TString> histName = {"gamma1_Energy",
                                     "gamma2_Energy",
                                     "gamma1_Energy_relative",
                                     "gamma2_Energy_relative",
                                     "pi0_InvariantMass",
                                     "pi0_momentum",
                                     "chi2",
                                     "prob",
                                     "pull_Egamma1",
                                     "pull_Egamma2"};

    std::map<TString, std::vector<Double_t>> histRange = {
        {"gamma1_Energy", {0., 1000.}},
        {"gamma2_Energy", {0., 1000.}},
        {"gamma1_Energy_relative", {-0.5, 0.5}},
        {"gamma2_Energy_relative", {-0.5, 0.5}},
        {"pi0_InvariantMass", {100., 170.}},
        {"pi0_momentum", {300, 700.}},
        {"chi2", {0., 20.}},
        {"prob", {0., 1.}},
        {"pull_Egamma1", {-5., 5.}},
        {"pull_Egamma2", {-5., 5.}}};

    for (TString name : histName)
    {
        canvasList[name] = new TCanvas(name, name, 750, 750);
        histInitList[name] = new TH1F(name, name, 100, histRange[name][0], histRange[name][1]);
        histAfterList[name] = new TH1F(name + "_after", name + "_after", 100, histRange[name][0], histRange[name][1]);
    }

    TRandom3 randGen;

    KLOE::neutralParticle
        pi0Gen,
        gamma1,
        gamma2,
        gamma1Lab,
        gamma2Lab,
        gamma1LabSmear,
        gamma2LabSmear;

    randGen.SetSeed(0);

    TVector3 zAxis(0., 0., 1.);

    // Generacja zgodna z rozpadem dwuciałowym
    // Określam własności pi0

    pi0Gen.fourMom[0] = 500.;
    pi0Gen.fourMom[1] = 0.;
    pi0Gen.fourMom[2] = 0.; // MeV/c
    pi0Gen.fourMom[3] = sqrt(pow(pi0Gen.fourMom[0], 2) +
                             pow(pi0Gen.fourMom[1], 2) +
                             pow(pi0Gen.fourMom[2], 2) +
                             pow(PhysicsConstants::mPi0, 2));

    Double_t beta = pi0Gen.fourMom[0] / pi0Gen.fourMom[3],
             gamma = 1. / sqrt(1. - beta * beta);

    pi0Gen.SetTotalVectorPhoton();

    Int_t maxPoints = 100000, count = 0;

    for (Int_t i = 0; i < maxPoints; i++)
    {
        // Generacja w układzie środka masy pi0
        // Rozpad pi0 na dwa fotony
        gamma1.fourMom[3] = PhysicsConstants::mPi0 / 2.;
        gamma2.fourMom[3] = PhysicsConstants::mPi0 / 2.;

        Double_t theta1 = randGen.Uniform(0., TMath::TwoPi()),
                 theta2 = TMath::Pi() + theta1;

        gamma1.fourMom[0] = gamma1.fourMom[3] * cos(theta1);
        gamma1.fourMom[1] = gamma1.fourMom[3] * sin(theta1);

        gamma2.fourMom[0] = gamma2.fourMom[3] * cos(theta2);
        gamma2.fourMom[1] = gamma2.fourMom[3] * sin(theta2);

        // Transformacja Lorentza do układu laboratorium

        gamma1Lab.fourMom[3] = gamma * (gamma1.fourMom[3] + beta * gamma1.fourMom[0]);
        gamma1Lab.fourMom[0] = gamma * (gamma1.fourMom[0] + beta * gamma1.fourMom[3]);
        gamma1Lab.fourMom[1] = gamma1.fourMom[1];

        gamma2Lab.fourMom[3] = gamma * (gamma2.fourMom[3] + beta * gamma2.fourMom[0]);
        gamma2Lab.fourMom[0] = gamma * (gamma2.fourMom[0] + beta * gamma2.fourMom[3]);
        gamma2Lab.fourMom[1] = gamma2.fourMom[1];

        Double_t mom1 = sqrt(pow(gamma1Lab.fourMom[0], 2) +
                             pow(gamma1Lab.fourMom[1], 2) +
                             pow(gamma1Lab.fourMom[2], 2));
        Double_t mom2 = sqrt(pow(gamma2Lab.fourMom[0], 2) +
                             pow(gamma2Lab.fourMom[1], 2) +
                             pow(gamma2Lab.fourMom[2], 2));

        Double_t numerator = gamma1Lab.fourMom[0] * gamma2Lab.fourMom[0] +
                             gamma1Lab.fourMom[1] * gamma2Lab.fourMom[1] +
                             gamma1Lab.fourMom[2] * gamma2Lab.fourMom[2];
        Double_t cosTheta = numerator / (mom1 * mom2);

        Double_t angle = acos(cosTheta);

        Double_t relErrEnergy = 0.05,
                 relErrEnergySmear = 0.05; // 5% energy resolution

        // Rozmycie energii zmierzonej
        gamma1LabSmear.fourMom[3] = gamma1Lab.fourMom[3] + randGen.Gaus(0.0, relErrEnergy * gamma1Lab.fourMom[3]);
        gamma2LabSmear.fourMom[3] = gamma2Lab.fourMom[3] + randGen.Gaus(0.0, relErrEnergy * gamma2Lab.fourMom[3]);

        Double_t invMass = 2 * gamma1LabSmear.fourMom[3] * gamma2LabSmear.fourMom[3] * (1. - cos(angle));

        Double_t pi0_momentum = sqrt(pow(gamma1LabSmear.fourMom[0] + gamma2LabSmear.fourMom[0], 2) +
                                     pow(gamma1LabSmear.fourMom[1] + gamma2LabSmear.fourMom[1], 2) +
                                     pow(gamma1LabSmear.fourMom[2] + gamma2LabSmear.fourMom[2], 2));

        invMass = sqrt(invMass);

        Double_t Param[3] = {gamma1LabSmear.fourMom[3], gamma2LabSmear.fourMom[3], angle},
                 Errors[3] = {relErrEnergySmear * gamma1Lab.fourMom[3], relErrEnergySmear * gamma2Lab.fourMom[3], 0.0};

        kinFit.ParameterInitialization(Param, Errors);

        _CHISQRMIN = kinFit.FitFunction();

        if (TMath::Prob(_CHISQRMIN, 1) > 0.04)
        {
            kinFit.GetResults(_X_min, _V_min, _X_init_min, _V_init);

            histInitList["gamma1_Energy"]->Fill(gamma1LabSmear.fourMom[3]);
            histInitList["gamma2_Energy"]->Fill(gamma2LabSmear.fourMom[3]);
            histInitList["pi0_InvariantMass"]->Fill(invMass);
            histInitList["gamma1_Energy_relative"]->Fill((gamma1LabSmear.fourMom[3] - gamma1Lab.fourMom[3]) / (gamma1Lab.fourMom[3]));
            histInitList["gamma2_Energy_relative"]->Fill((gamma2LabSmear.fourMom[3] - gamma2Lab.fourMom[3]) / (gamma2Lab.fourMom[3]));
            histInitList["pi0_momentum"]->Fill(pi0_momentum);

            Double_t invMassCorr = 2 * _X_min[0] * _X_min[1] * (1. - cos(angle));
            invMassCorr = sqrt(invMassCorr);

            Double_t pull1 = (_X_init_min[0] - _X_min[0]) / sqrt(_V_init[0][0] - _V_min[0][0]);
            Double_t pull2 = (_X_init_min[1] - _X_min[1]) / sqrt(_V_init[1][1] - _V_min[1][1]);

            histAfterList["gamma1_Energy"]->Fill(_X_min[0]);
            histAfterList["gamma2_Energy"]->Fill(_X_min[1]);
            histAfterList["pi0_InvariantMass"]->Fill(invMassCorr);
            histAfterList["gamma1_Energy_relative"]->Fill((_X_min[0] - gamma1Lab.fourMom[3]) / (gamma1Lab.fourMom[3]));
            histAfterList["gamma2_Energy_relative"]->Fill((_X_min[1] - gamma2Lab.fourMom[3]) / (gamma2Lab.fourMom[3]));
            histAfterList["chi2"]->Fill(_CHISQRMIN);
            histAfterList["prob"]->Fill(TMath::Prob(_CHISQRMIN, 1));
            histAfterList["pull_Egamma1"]->Fill(pull1);
            histAfterList["pull_Egamma2"]->Fill(pull2);
        }

        if (TMath::Prob(_CHISQRMIN, 1) == 0)
            count++;

        gamma1.fourMom[3] = 0.;
        gamma2.fourMom[3] = 0.;
    }

    std::cout << "Liczba punktów z prob=0: " << count << " z " << maxPoints << " (" << 100. * count / maxPoints << "%)" << std::endl;

    // Wyświetlenie histogramów
    for (const auto &canvasPair : canvasList)
    {
        if (canvasPair.first == "pull_Egamma1" || canvasPair.first == "pull_Egamma2")
            histAfterList[canvasPair.first]->Fit("gaus");

        canvasPair.second->cd();
        histAfterList[canvasPair.first]->SetLineColor(kRed);
        histAfterList[canvasPair.first]->SetLineWidth(2);
        histInitList[canvasPair.first]->SetLineColor(kBlue);
        histInitList[canvasPair.first]->SetLineWidth(2);

        histAfterList[canvasPair.first]->GetYaxis()->SetRangeUser(0., TMath::Max(histAfterList[canvasPair.first]->GetMaximum(), histInitList[canvasPair.first]->GetMaximum()) * 1.2);

        histAfterList[canvasPair.first]->Draw();
        histInitList[canvasPair.first]->Draw("SAME");

        canvasPair.second->Print(canvasPair.first + ".png");
    }

    return 0;
}