/*
 * Przykład aplikacji desktopowej używającej ROOT GUI
 * Zintegrowany z istniejącym kodem analizy KLOE
 * 
 * Author: Szymon Gamrat
 * Date: 2024
 */

#ifdef ROOT_AVAILABLE  // Kompiluje tylko gdy ROOT jest dostępny

#include <TApplication.h>
#include <TGWindow.h>
#include <TGFrame.h>
#include <TGButton.h>
#include <TGMenu.h>
#include <TGMenuBar.h>
#include <TGLabel.h>
#include <TGTextEdit.h>
#include <TGProgressBar.h>
#include <TGStatusBar.h>
#include <TGFileDialog.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TLegend.h>
#include <TStyle.h>
#include <iostream>

// Include z istniejącego projektu KLOE (jeśli dostępne)
// #include "../Include/klspm00.hpp"

/**
 * @brief Główne okno GUI używające ROOT
 */
class KLOERootMainFrame : public TGMainFrame 
{
private:
    TGMenuBar*          fMenuBar;
    TGHorizontalFrame*  fButtonFrame;
    TGVerticalFrame*    fMainFrame;
    TGTextEdit*         fLogText;
    TGLabel*            fStatusLabel;
    TGProgressBar*      fProgressBar;
    
    TGTextButton*       fGenVarsButton;
    TGTextButton*       fKchRecButton;
    TGTextButton*       fNeutRecButton;
    TGTextButton*       fInterfButton;
    TGTextButton*       fPlotsButton;
    
    TCanvas*            fCanvas;

public:
    KLOERootMainFrame(const TGWindow* p, UInt_t w, UInt_t h);
    virtual ~KLOERootMainFrame();
    
    // Obsługa zdarzeń
    void DoGenVars();
    void DoKchRec();
    void DoNeutRec();
    void DoInterf();
    void DoPlots();
    void DoExit();
    void DoOpenFile();
    void DoAbout();
    
    // ROOT ClassDef
    ClassDef(KLOERootMainFrame, 0)
};

/**
 * @brief Konstruktor głównego okna
 */
KLOERootMainFrame::KLOERootMainFrame(const TGWindow* p, UInt_t w, UInt_t h)
    : TGMainFrame(p, w, h)
{
    SetWindowName("KLOE Physics Analysis - ROOT GUI");
    
    // Tworzenie menu
    fMenuBar = new TGMenuBar(this, 1, 1, kHorizontalFrame);
    
    TGPopupMenu* fMenuFile = new TGPopupMenu(gClient->GetRoot());
    fMenuFile->AddEntry("&Otwórz plik ROOT", 1);
    fMenuFile->AddSeparator();
    fMenuFile->AddEntry("&Wyjście", 2);
    
    TGPopupMenu* fMenuAnalysis = new TGPopupMenu(gClient->GetRoot());
    fMenuAnalysis->AddEntry("Zmienne &generowane", 10);
    fMenuAnalysis->AddEntry("Rekonstrukcja K→π+π-", 11);
    fMenuAnalysis->AddEntry("Rekonstrukcja K→π0π0", 12);
    fMenuAnalysis->AddEntry("&Interferometria", 13);
    fMenuAnalysis->AddEntry("&Wykresy", 14);
    
    TGPopupMenu* fMenuHelp = new TGPopupMenu(gClient->GetRoot());
    fMenuHelp->AddEntry("&O programie", 20);
    
    fMenuFile->Connect("Activated(Int_t)", "KLOERootMainFrame", this, "HandleMenu(Int_t)");
    fMenuAnalysis->Connect("Activated(Int_t)", "KLOERootMainFrame", this, "HandleMenu(Int_t)");
    fMenuHelp->Connect("Activated(Int_t)", "KLOERootMainFrame", this, "HandleMenu(Int_t)");
    
    fMenuBar->AddPopup("&Plik", fMenuFile, new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0));
    fMenuBar->AddPopup("&Analiza", fMenuAnalysis, new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0));
    fMenuBar->AddPopup("&Pomoc", fMenuHelp, new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0));
    
    AddFrame(fMenuBar, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 0, 0, 1, 1));
    
    // Główna ramka
    fMainFrame = new TGVerticalFrame(this);
    
    // Nagłówek
    TGLabel* titleLabel = new TGLabel(fMainFrame, "KLOE Physics Analysis");
    titleLabel->SetTextFont("-*-helvetica-bold-*-*-*-16-*-*-*-*-*-*-*");
    fMainFrame->AddFrame(titleLabel, new TGLayoutHints(kLHintsCenterX, 10, 10, 10, 10));
    
    // Ramka przycisków
    fButtonFrame = new TGHorizontalFrame(fMainFrame);
    
    fGenVarsButton = new TGTextButton(fButtonFrame, "Zmienne generowane");
    fGenVarsButton->Connect("Clicked()", "KLOERootMainFrame", this, "DoGenVars()");
    fButtonFrame->AddFrame(fGenVarsButton, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
    
    fKchRecButton = new TGTextButton(fButtonFrame, "K→π+π-");
    fKchRecButton->Connect("Clicked()", "KLOERootMainFrame", this, "DoKchRec()");
    fButtonFrame->AddFrame(fKchRecButton, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
    
    fNeutRecButton = new TGTextButton(fButtonFrame, "K→π0π0");
    fNeutRecButton->Connect("Clicked()", "KLOERootMainFrame", this, "DoNeutRec()");
    fButtonFrame->AddFrame(fNeutRecButton, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
    
    fInterfButton = new TGTextButton(fButtonFrame, "Interferometria");
    fInterfButton->Connect("Clicked()", "KLOERootMainFrame", this, "DoInterf()");
    fButtonFrame->AddFrame(fInterfButton, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
    
    fPlotsButton = new TGTextButton(fButtonFrame, "Wykresy");
    fPlotsButton->Connect("Clicked()", "KLOERootMainFrame", this, "DoPlots()");
    fButtonFrame->AddFrame(fPlotsButton, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
    
    fMainFrame->AddFrame(fButtonFrame, new TGLayoutHints(kLHintsExpandX, 10, 10, 5, 5));
    
    // Obszar logów
    TGLabel* logLabel = new TGLabel(fMainFrame, "Logi wykonania:");
    fMainFrame->AddFrame(logLabel, new TGLayoutHints(kLHintsLeft, 10, 10, 10, 5));
    
    fLogText = new TGTextEdit(fMainFrame, 500, 100);
    fLogText->SetReadOnly(kTRUE);
    fMainFrame->AddFrame(fLogText, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10, 10, 0, 5));
    
    // Pasek postępu
    fProgressBar = new TGProgressBar(fMainFrame, TGProgressBar::kStandard, 200);
    fProgressBar->SetRange(0, 100);
    fProgressBar->SetPosition(0);
    fMainFrame->AddFrame(fProgressBar, new TGLayoutHints(kLHintsExpandX, 10, 10, 5, 5));
    
    // Status
    fStatusLabel = new TGLabel(fMainFrame, "Gotowe");
    fMainFrame->AddFrame(fStatusLabel, new TGLayoutHints(kLHintsLeft, 10, 10, 5, 10));
    
    AddFrame(fMainFrame, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
    
    // Mapowanie okna
    MapSubwindows();
    Resize(GetDefaultSize());
    MapWindow();
}

/**
 * @brief Destruktor
 */
KLOERootMainFrame::~KLOERootMainFrame() 
{
    Cleanup();
}

/**
 * @brief Analiza zmiennych generowanych
 */
void KLOERootMainFrame::DoGenVars() 
{
    fLogText->AddText("Uruchamianie analizy zmiennych generowanych...\n");
    fStatusLabel->SetText("Wykonywanie: Analiza zmiennych generowanych");
    fProgressBar->SetPosition(25);
    
    // Tutaj można dodać wywołanie funkcji z projektu KLOE:
    // GenVars_main(chain, eventAnalysis, dataTypeOpt, physConst);
    
    fLogText->AddText("Analiza zmiennych generowanych zakończona!\n");
    fProgressBar->SetPosition(100);
    fStatusLabel->SetText("Gotowe");
}

/**
 * @brief Rekonstrukcja naładowanych kaonów
 */
void KLOERootMainFrame::DoKchRec() 
{
    fLogText->AddText("Uruchamianie rekonstrukcji K->π+π-...\n");
    fStatusLabel->SetText("Wykonywanie: Rekonstrukcja K->π+π-");
    fProgressBar->SetPosition(50);
    
    // KchRec_main(chain, eventAnalysis, dataTypeOpt, physConst);
    
    fLogText->AddText("Rekonstrukcja K->π+π- zakończona!\n");
    fProgressBar->SetPosition(100);
    fStatusLabel->SetText("Gotowe");
}

/**
 * @brief Rekonstrukcja neutralnych kaonów
 */
void KLOERootMainFrame::DoNeutRec() 
{
    fLogText->AddText("Uruchamianie rekonstrukcji K->π0π0...\n");
    fStatusLabel->SetText("Wykonywanie: Rekonstrukcja K->π0π0");
    fProgressBar->SetPosition(75);
    
    // Neutrec_main(chain, eventAnalysis, dataTypeOpt, physConst);
    
    fLogText->AddText("Rekonstrukcja K->π0π0 zakończona!\n");
    fProgressBar->SetPosition(100);
    fStatusLabel->SetText("Gotowe");
}

/**
 * @brief Analiza interferometrii z istniejącego kodu
 */
void KLOERootMainFrame::DoInterf() 
{
    fLogText->AddText("Uruchamianie analizy interferometrii...\n");
    fStatusLabel->SetText("Wykonywanie: Analiza interferometrii");
    fProgressBar->SetPosition(80);
    
    // Kod oparty na interf_func_draw.cpp z projektu KLOE
    fCanvas = new TCanvas("canvas_interf", "Interference Function", 750, 750);
    fCanvas->SetMargin(0.15, 0.15, 0.15, 0.1);
    
    // Funkcja interferometrii (z interf_function.h)
    TF1 *func_with_im = new TF1("Interference function2", "1 + 2*[0]*cos(x) - 2*[1]*sin(x)", -20.0, 20.0);
    func_with_im->SetParameter(0, 0.005);
    func_with_im->SetParameter(1, 0.05);
    
    TF1 *func_wo_im = new TF1("Interference function1", "1 + 2*[0]*cos(x) - 2*[1]*sin(x)", -20.0, 20.0);
    func_wo_im->SetParameter(0, 0);
    func_wo_im->SetParameter(1, 0);
    
    fCanvas->SetGrid(1, 1);
    gStyle->SetGridStyle(3);
    gStyle->SetGridWidth(2);
    
    func_with_im->SetNpx(1E6);
    func_with_im->SetLineColor(kRed);
    func_with_im->SetLineWidth(3);
    
    func_wo_im->SetNpx(1E6);
    func_wo_im->SetLineColor(kBlack);
    func_wo_im->SetLineWidth(3);
    
    func_with_im->SetTitle("");
    func_with_im->GetYaxis()->SetRangeUser(0, 0.4);
    func_with_im->GetYaxis()->SetTitle("I [-]");
    func_with_im->GetYaxis()->CenterTitle(1);
    func_with_im->GetXaxis()->SetTitle("#Deltat [#tau_{S}]");
    func_with_im->GetXaxis()->CenterTitle(1);
    
    func_with_im->Draw();
    func_wo_im->Draw("SAME");
    
    TLegend *legend = new TLegend(0.6, 0.73, 0.9, 0.9, "");
    legend->AddEntry(func_wo_im, "#varepsilon'/#varepsilon = 0", "l");
    legend->AddEntry(func_with_im, "Re(#varepsilon'/#varepsilon) = 5#times10^{-3}, Im(#varepsilon'/#varepsilon) = 5#times10^{-2}", "l");
    legend->SetTextSize(0.025);
    legend->SetTextFont(42);
    legend->Draw();
    
    fCanvas->Update();
    
    fLogText->AddText("Analiza interferometrii zakończona! Sprawdź okno wykresu.\n");
    fProgressBar->SetPosition(100);
    fStatusLabel->SetText("Gotowe");
}

/**
 * @brief Tworzenie wykresów
 */
void KLOERootMainFrame::DoPlots() 
{
    fLogText->AddText("Tworzenie wykresów i wizualizacji...\n");
    fStatusLabel->SetText("Wykonywanie: Tworzenie wykresów");
    fProgressBar->SetPosition(90);
    
    // Tutaj można zintegrować inne wykresy z projektu KLOE
    
    fLogText->AddText("Wykresy utworzone!\n");
    fProgressBar->SetPosition(100);
    fStatusLabel->SetText("Gotowe");
}

/**
 * @brief Obsługa menu
 */
void KLOERootMainFrame::HandleMenu(Int_t id) 
{
    switch (id) {
        case 1:  DoOpenFile(); break;
        case 2:  DoExit(); break;
        case 10: DoGenVars(); break;
        case 11: DoKchRec(); break;
        case 12: DoNeutRec(); break;
        case 13: DoInterf(); break;
        case 14: DoPlots(); break;
        case 20: DoAbout(); break;
        default: break;
    }
}

void KLOERootMainFrame::DoOpenFile() 
{
    const char *filetypes[] = {"ROOT files", "*.root", "All files", "*", 0, 0};
    TGFileInfo fi;
    fi.fFileTypes = filetypes;
    new TGFileDialog(gClient->GetRoot(), this, kFDOpen, &fi);
    
    if (fi.fFilename) {
        fLogText->AddText(Form("Otwarto plik: %s\n", fi.fFilename));
        fStatusLabel->SetText(Form("Plik załadowany: %s", fi.fFilename));
    }
}

void KLOERootMainFrame::DoExit() 
{
    gApplication->Terminate(0);
}

void KLOERootMainFrame::DoAbout() 
{
    new TGMsgBox(gClient->GetRoot(), this, "O programie",
                 "KLOE Physics Analysis Desktop App\n\n"
                 "Aplikacja do analizy danych z eksperymentu KLOE\n"
                 "Przykład implementacji GUI w ROOT\n\n"
                 "Autor: Szymon Gamrat\n"
                 "Data: 2024",
                 kMBIconAsterisk, kMBOk);
}

// ROOT ClassImp
ClassImp(KLOERootMainFrame)

/**
 * @brief Funkcja główna dla aplikacji ROOT GUI
 */
int main(int argc, char **argv) 
{
    TApplication theApp("KLOEApp", &argc, argv);
    
    KLOERootMainFrame *mainframe = new KLOERootMainFrame(gClient->GetRoot(), 800, 600);
    
    theApp.Run();
    
    return 0;
}

#else

#include <iostream>

int main(int argc, char *argv[])
{
    std::cout << "=== Przykład aplikacji ROOT GUI ===" << std::endl;
    std::cout << "Ten przykład wymaga zainstalowanego ROOT framework." << std::endl;
    std::cout << "Aby skompilować:" << std::endl;
    std::cout << "1. Zainstaluj ROOT (https://root.cern/install/)" << std::endl;
    std::cout << "2. Zdefiniuj ROOT_AVAILABLE przy kompilacji" << std::endl;
    std::cout << "3. Użyj CMake z find_package(ROOT)" << std::endl;
    return 0;
}

#endif