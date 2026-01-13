#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TPad.h>
#include <TLine.h>

#include <const.h>

#include <vector>
#include <string>
#include <iostream>

int main()
{
    KLOE::setGlobalStyle();

    TString plot_dir = "theoretical_plots/";
    TString file_path = Paths::img_dir + plot_dir + "interf_func_parameters_theory.root";

    TFile *file = TFile::Open(file_path);
    if (!file || file->IsZombie())
    {
        std::cerr << "Error: Cannot open file " << file_path << std::endl;
        return 1;
    }

    // Przechowuj informacje o folderach
    std::vector<TString> folders_T300;  // Wszystkie T = 300
    std::vector<TString> folders_other; // Inne kombinacje T

    // Przejdź przez wszystkie klucze w pliku
    TIter next(file->GetListOfKeys());
    TKey *key;
    while ((key = (TKey *)next()))
    {
        if (strcmp(key->GetClassName(), "TDirectoryFile") == 0)
        {
            TString dirName = key->GetName();

            // Sprawdź czy wszystkie T = 300
            if (dirName.Contains("T00pm_300") && dirName.Contains("Tpm00_300") &&
                dirName.Contains("Tpmpm1_300") && dirName.Contains("Tpmpm2_300"))
            {
                folders_T300.push_back(dirName);
            }
            else
            {
                folders_other.push_back(dirName);
            }
        }
    }

    std::cout << "Folders with all T=300: " << folders_T300.size() << std::endl;
    std::cout << "Folders with other T: " << folders_other.size() << std::endl;

    // Funkcja pomocnicza do rysowania zestawu funkcji
    auto drawComparison = [&](std::vector<TString> &folders, TString suffix)
    {
        if (folders.size() == 0)
        {
            std::cout << "No folders for " << suffix << std::endl;
            return;
        }

        TString funcNames[] = {"func1D_RA", "func1D_RB", "func1D_RC"};
        TString titles[] = {"R_{A}(t_{1})", "R_{B}(t_{1})", "R_{C}(t_{1})"};

        // Podziel foldery na grupy po 4
        Int_t foldersPerCanvas = 4;
        Int_t nCanvases = (folders.size() + foldersPerCanvas - 1) / foldersPerCanvas;

        for (Int_t iCanvas = 0; iCanvas < nCanvases; iCanvas++)
        {
            Int_t startIdx = iCanvas * foldersPerCanvas;
            Int_t endIdx = TMath::Min((Int_t)folders.size(), (iCanvas + 1) * foldersPerCanvas);
            Int_t nFolders = endIdx - startIdx;

            TCanvas *canvas = new TCanvas(Form("canvas_%s_%d", suffix.Data(), iCanvas),
                                          Form("Comparison %s - part %d", suffix.Data(), iCanvas + 1),
                                          1800, 400 * nFolders);
            canvas->Divide(3, nFolders);

            // Dla każdej funkcji (RA, RB, RC) i każdego folderu
            for (Int_t iRow = 0; iRow < nFolders; iRow++)
            {
                Int_t iDir = startIdx + iRow;

                for (Int_t iFunc = 0; iFunc < 3; iFunc++)
                {
                    Int_t padNum = iRow * 3 + iFunc + 1;
                    canvas->cd(padNum);
                    gPad->SetMargin(0.12, 0.05, 0.12, 0.08);
                    gPad->SetGrid(1, 1);

                    TString path = folders[iDir] + "/" + funcNames[iFunc];
                    TF1 *func = (TF1 *)file->Get(path);

                    if (func)
                    {
                        func->SetLineColor(kBlue);
                        func->SetLineWidth(2);
                        func->SetNpx(1000);

                        // Wyciągnij parametry z nazwy folderu
                        TString label = folders[iDir];
                        label.ReplaceAll("Re_", "Re=");
                        label.ReplaceAll("_Im_", " Im=");
                        label.ReplaceAll("_T00pm_", " T_{00#pm}=");
                        label.ReplaceAll("_Tpm00_", " T_{#pm00}=");
                        label.ReplaceAll("_Tpmpm1_", " T_{#pm#pm}^{1}=");
                        label.ReplaceAll("_Tpmpm2_", " T_{#pm#pm}^{2}=");

                        func->SetTitle(titles[iFunc] + " - " + label);
                        func->GetXaxis()->SetTitle("t_{1} [#tau_{S}]");
                        func->GetYaxis()->SetTitle(titles[iFunc] + " [-]");
                        func->Draw();

                        // Dodaj linię referencyjną y=1 dla RC
                        if (iFunc == 2)
                        {
                            TLine *line = new TLine(0, 1, 20, 1);
                            line->SetLineStyle(2);
                            line->SetLineColor(kBlack);
                            line->SetLineWidth(2);
                            line->Draw("SAME");
                        }
                    }
                }
            }

            canvas->SaveAs(Paths::img_dir + plot_dir + Form("comparison_%s_part%d", suffix.Data(), iCanvas + 1) + ".svg");
            delete canvas;
        }
    };

    // Funkcja filtrująca funkcje według max/min
    auto drawComparisonFiltered = [&](std::vector<TString> &folders, TString suffix, TString filterType)
    {
        if (folders.size() == 0)
        {
            std::cout << "No folders for " << suffix << std::endl;
            return;
        }

        TString funcNames[] = {"func1D_RA", "func1D_RB", "func1D_RC"};
        TString titles[] = {"R_{A}(t_{1})", "R_{B}(t_{1})", "R_{C}(t_{1})"};

        // Dla każdej funkcji oblicz max/min i odfiltruj outlierów
        std::vector<TString> filtered_folders;

        for (Int_t iFunc = 0; iFunc < 3; iFunc++)
        {
            std::vector<Double_t> max_values, min_values;

            // Oblicz max i min dla każdego folderu
            for (size_t iDir = 0; iDir < folders.size(); iDir++)
            {
                TString path = folders[iDir] + "/" + funcNames[iFunc];
                TF1 *func = (TF1 *)file->Get(path);
                if (func)
                {
                    Double_t func_max = -1e10, func_min = 1e10;
                    for (Double_t x = 0; x < 20; x += 0.1)
                    {
                        Double_t y = func->Eval(x);
                        if (y > func_max)
                            func_max = y;
                        if (y < func_min)
                            func_min = y;
                    }
                    max_values.push_back(func_max);
                    min_values.push_back(func_min);
                }
            }

            // Znajdź index folderu z ekstremalną wartością
            Int_t idx_to_remove = -1;
            
            if (iFunc == 0) // Tylko raz, dla pierwszej funkcji
            {
                if (filterType == "max")
                {
                    // Znajdź folder z najwyższym max
                    Double_t highest_max = -1e10;
                    for (size_t iDir = 0; iDir < folders.size(); iDir++)
                    {
                        if (max_values[iDir] > highest_max)
                        {
                            highest_max = max_values[iDir];
                            idx_to_remove = iDir;
                        }
                    }
                    std::cout << funcNames[iFunc] << " - Removing folder with highest max=" << highest_max << std::endl;
                }
                else if (filterType == "min")
                {
                    // Znajdź folder z najniższym min
                    Double_t lowest_min = 1e10;
                    for (size_t iDir = 0; iDir < folders.size(); iDir++)
                    {
                        if (min_values[iDir] < lowest_min)
                        {
                            lowest_min = min_values[iDir];
                            idx_to_remove = iDir;
                        }
                    }
                    std::cout << funcNames[iFunc] << " - Removing folder with lowest min=" << lowest_min << std::endl;
                }

                // Dodaj wszystkie foldery oprócz tego do usunięcia
                for (size_t iDir = 0; iDir < folders.size(); iDir++)
                {
                    if ((Int_t)iDir != idx_to_remove)
                    {
                        filtered_folders.push_back(folders[iDir]);
                        std::cout << "Keeping: " << folders[iDir] << std::endl;
                    }
                    else
                    {
                        std::cout << "Filtering out: " << folders[iDir] << " (max=" << max_values[iDir]
                                  << ", min=" << min_values[iDir] << ")" << std::endl;
                    }
                }
            }
        }

        std::cout << "Filtered folders: " << filtered_folders.size() << " / " << folders.size() << std::endl;

        if (filtered_folders.size() == 0)
        {
            std::cout << "All folders filtered out!" << std::endl;
            return;
        }

        // Podziel przefiltrowane foldery na grupy po 4
        Int_t foldersPerCanvas = 4;
        Int_t nCanvases = (filtered_folders.size() + foldersPerCanvas - 1) / foldersPerCanvas;

        for (Int_t iCanvas = 0; iCanvas < nCanvases; iCanvas++)
        {
            Int_t startIdx = iCanvas * foldersPerCanvas;
            Int_t endIdx = TMath::Min((Int_t)filtered_folders.size(), (iCanvas + 1) * foldersPerCanvas);
            Int_t nFolders = endIdx - startIdx;

            TCanvas *canvas = new TCanvas(Form("canvas_%s_%s_%d", suffix.Data(), filterType.Data(), iCanvas),
                                          Form("Comparison %s %s - part %d", suffix.Data(), filterType.Data(), iCanvas + 1),
                                          1800, 400 * nFolders);
            canvas->Divide(3, nFolders);

            // Dla każdej funkcji (RA, RB, RC) i każdego folderu
            for (Int_t iRow = 0; iRow < nFolders; iRow++)
            {
                Int_t iDir = startIdx + iRow;

                for (Int_t iFunc = 0; iFunc < 3; iFunc++)
                {
                    Int_t padNum = iRow * 3 + iFunc + 1;
                    canvas->cd(padNum);
                    gPad->SetMargin(0.12, 0.05, 0.12, 0.08);
                    gPad->SetGrid(1, 1);

                    TString path = filtered_folders[iDir] + "/" + funcNames[iFunc];
                    TF1 *func = (TF1 *)file->Get(path);

                    if (func)
                    {
                        func->SetLineColor(kBlue);
                        func->SetLineWidth(2);
                        func->SetNpx(1000);

                        // Wyciągnij parametry z nazwy folderu
                        TString label = filtered_folders[iDir];
                        label.ReplaceAll("Re_", "Re=");
                        label.ReplaceAll("_Im_", " Im=");
                        label.ReplaceAll("_T00pm_", " T_{00#pm}=");
                        label.ReplaceAll("_Tpm00_", " T_{#pm00}=");
                        label.ReplaceAll("_Tpmpm1_", " T_{#pm#pm}^{1}=");
                        label.ReplaceAll("_Tpmpm2_", " T_{#pm#pm}^{2}=");

                        func->SetTitle(titles[iFunc] + " (filtered " + filterType + ") - " + label);
                        func->GetXaxis()->SetTitle("t_{1} [#tau_{S}]");
                        func->GetYaxis()->SetTitle(titles[iFunc] + " [-]");
                        func->Draw();

                        // Dodaj linię referencyjną y=1 dla RC
                        if (iFunc == 2)
                        {
                            TLine *line = new TLine(0, 1, 20, 1);
                            line->SetLineStyle(2);
                            line->SetLineColor(kBlack);
                            line->SetLineWidth(2);
                            line->Draw("SAME");
                        }
                    }
                }
            }

            canvas->SaveAs(Paths::img_dir + plot_dir + Form("comparison_%s_%s_part%d", suffix.Data(), filterType.Data(), iCanvas + 1) + ".svg");
            delete canvas;
        }
    };

    // Rysuj porównania
    drawComparison(folders_T300, "T300");
    drawComparison(folders_other, "T_other");

    // Rysuj przefiltrowane porównania - usuń outlierów z max > 110% średniej
    drawComparisonFiltered(folders_T300, "T300", "max");
    drawComparisonFiltered(folders_other, "T_other", "max");

    // Rysuj przefiltrowane porównania - usuń outlierów z min > 110% średniej
    drawComparisonFiltered(folders_T300, "T300", "min");
    drawComparisonFiltered(folders_other, "T_other", "min");

    file->Close();

    std::cout << "Comparison plots saved!" << std::endl;

    return 0;
}
