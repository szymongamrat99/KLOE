/*
 * Przykład aplikacji desktopowej Qt dla analizy KLOE
 * 
 * To jest przykład jak stworzyć GUI używając Qt
 * Można go skompilować gdy Qt jest zainstalowane
 * 
 * Author: Szymon Gamrat  
 * Date: 2024
 */

#ifdef QT_AVAILABLE  // Kompiluje tylko gdy Qt jest dostępne

#include <QApplication>
#include <QMainWindow>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QWidget>
#include <QPushButton>
#include <QLabel>
#include <QMenuBar>
#include <QMenu>
#include <QAction>
#include <QTextEdit>
#include <QProgressBar>
#include <QStatusBar>
#include <QMessageBox>
#include <QFileDialog>
#include <QSplitter>

/**
 * @brief Główne okno aplikacji KLOE z Qt
 */
class KLOEMainWindow : public QMainWindow
{
    Q_OBJECT

private:
    QWidget* centralWidget;
    QVBoxLayout* mainLayout;
    QHBoxLayout* buttonLayout;
    QTextEdit* logTextEdit;
    QPushButton* genVarsButton;
    QPushButton* kchRecButton;
    QPushButton* neutRecButton;
    QPushButton* interfButton;
    QPushButton* plotsButton;
    QProgressBar* progressBar;
    QLabel* statusLabel;

public:
    KLOEMainWindow(QWidget *parent = nullptr) : QMainWindow(parent)
    {
        setWindowTitle("KLOE Physics Analysis - Desktop App");
        setMinimumSize(800, 600);
        
        setupUI();
        setupMenuBar();
        setupStatusBar();
        connectSignals();
    }

private slots:
    void onGenVarsClicked() 
    {
        logTextEdit->append("Uruchamianie analizy zmiennych generowanych...");
        statusLabel->setText("Wykonywanie: Analiza zmiennych generowanych");
        progressBar->setValue(20);
        
        // Tutaj można dodać wywołanie funkcji z projektu KLOE:
        // GenVars_main(chain, eventAnalysis, dataTypeOpt, physConst);
        
        QMessageBox::information(this, "Info", "Analiza zmiennych generowanych zakończona!");
        progressBar->setValue(100);
        statusLabel->setText("Gotowe");
    }
    
    void onKchRecClicked()
    {
        logTextEdit->append("Uruchamianie rekonstrukcji K->π+π-...");
        statusLabel->setText("Wykonywanie: Rekonstrukcja K->π+π-");
        progressBar->setValue(40);
        
        // KchRec_main(chain, eventAnalysis, dataTypeOpt, physConst);
        
        QMessageBox::information(this, "Info", "Rekonstrukcja K->π+π- zakończona!");
        progressBar->setValue(100);
        statusLabel->setText("Gotowe");
    }
    
    void onNeutRecClicked()
    {
        logTextEdit->append("Uruchamianie rekonstrukcji K->π0π0...");
        statusLabel->setText("Wykonywanie: Rekonstrukcja K->π0π0");
        progressBar->setValue(60);
        
        // Neutrec_main(chain, eventAnalysis, dataTypeOpt, physConst);
        
        QMessageBox::information(this, "Info", "Rekonstrukcja K->π0π0 zakończona!");
        progressBar->setValue(100);
        statusLabel->setText("Gotowe");
    }
    
    void onInterfClicked()
    {
        logTextEdit->append("Uruchamianie analizy interferometrii...");
        statusLabel->setText("Wykonywanie: Analiza interferometrii");
        progressBar->setValue(80);
        
        // Można dodać kod z interf_func_draw.cpp
        
        QMessageBox::information(this, "Info", "Analiza interferometrii zakończona!");
        progressBar->setValue(100);
        statusLabel->setText("Gotowe");
    }
    
    void onPlotsClicked()
    {
        logTextEdit->append("Tworzenie wykresów i wizualizacji...");
        statusLabel->setText("Wykonywanie: Tworzenie wykresów");
        progressBar->setValue(90);
        
        // Tutaj można zintegrować z ROOT TCanvas lub Qt Charts
        
        QMessageBox::information(this, "Info", "Wykresy utworzone!");
        progressBar->setValue(100);
        statusLabel->setText("Gotowe");
    }
    
    void openFile()
    {
        QString fileName = QFileDialog::getOpenFileName(this,
            "Otwórz plik danych", "", "ROOT Files (*.root);;All Files (*)");
        
        if (!fileName.isEmpty()) {
            logTextEdit->append("Otwarto plik: " + fileName);
            statusLabel->setText("Plik załadowany: " + fileName);
        }
    }
    
    void about()
    {
        QMessageBox::about(this, "O programie",
            "KLOE Physics Analysis Desktop App\n\n"
            "Aplikacja do analizy danych z eksperymentu KLOE\n"
            "Przykład implementacji GUI w Qt\n\n"
            "Autor: Szymon Gamrat\n"
            "Data: 2024");
    }

private:
    void setupUI()
    {
        centralWidget = new QWidget(this);
        setCentralWidget(centralWidget);
        
        mainLayout = new QVBoxLayout(centralWidget);
        
        // Nagłówek
        QLabel* titleLabel = new QLabel("KLOE Physics Analysis");
        titleLabel->setStyleSheet("font-size: 18px; font-weight: bold; margin: 10px;");
        titleLabel->setAlignment(Qt::AlignCenter);
        mainLayout->addWidget(titleLabel);
        
        // Przyciski analizy
        buttonLayout = new QHBoxLayout();
        
        genVarsButton = new QPushButton("Zmienne generowane");
        kchRecButton = new QPushButton("Rekonstrukcja K→π+π-");
        neutRecButton = new QPushButton("Rekonstrukcja K→π0π0");
        interfButton = new QPushButton("Interferometria");
        plotsButton = new QPushButton("Wykresy");
        
        buttonLayout->addWidget(genVarsButton);
        buttonLayout->addWidget(kchRecButton);
        buttonLayout->addWidget(neutRecButton);
        buttonLayout->addWidget(interfButton);
        buttonLayout->addWidget(plotsButton);
        
        mainLayout->addLayout(buttonLayout);
        
        // Obszar logów
        QLabel* logLabel = new QLabel("Logi wykonania:");
        mainLayout->addWidget(logLabel);
        
        logTextEdit = new QTextEdit();
        logTextEdit->setMaximumHeight(200);
        logTextEdit->setReadOnly(true);
        mainLayout->addWidget(logTextEdit);
        
        // Pasek postępu
        progressBar = new QProgressBar();
        progressBar->setRange(0, 100);
        progressBar->setValue(0);
        mainLayout->addWidget(progressBar);
    }
    
    void setupMenuBar()
    {
        QMenuBar* menuBar = this->menuBar();
        
        // Menu Plik
        QMenu* fileMenu = menuBar->addMenu("Plik");
        
        QAction* openAction = fileMenu->addAction("Otwórz");
        openAction->setShortcut(QKeySequence::Open);
        connect(openAction, &QAction::triggered, this, &KLOEMainWindow::openFile);
        
        fileMenu->addSeparator();
        
        QAction* exitAction = fileMenu->addAction("Wyjście");
        exitAction->setShortcut(QKeySequence::Quit);
        connect(exitAction, &QAction::triggered, this, &QWidget::close);
        
        // Menu Analiza
        QMenu* analysisMenu = menuBar->addMenu("Analiza");
        
        QAction* genVarsAction = analysisMenu->addAction("Zmienne generowane");
        connect(genVarsAction, &QAction::triggered, this, &KLOEMainWindow::onGenVarsClicked);
        
        QAction* kchRecAction = analysisMenu->addAction("Rekonstrukcja K→π+π-");
        connect(kchRecAction, &QAction::triggered, this, &KLOEMainWindow::onKchRecClicked);
        
        QAction* neutRecAction = analysisMenu->addAction("Rekonstrukcja K→π0π0");
        connect(neutRecAction, &QAction::triggered, this, &KLOEMainWindow::onNeutRecClicked);
        
        // Menu Pomoc
        QMenu* helpMenu = menuBar->addMenu("Pomoc");
        
        QAction* aboutAction = helpMenu->addAction("O programie");
        connect(aboutAction, &QAction::triggered, this, &KLOEMainWindow::about);
    }
    
    void setupStatusBar()
    {
        statusLabel = new QLabel("Gotowe");
        statusBar()->addWidget(statusLabel);
    }
    
    void connectSignals()
    {
        connect(genVarsButton, &QPushButton::clicked, this, &KLOEMainWindow::onGenVarsClicked);
        connect(kchRecButton, &QPushButton::clicked, this, &KLOEMainWindow::onKchRecClicked);
        connect(neutRecButton, &QPushButton::clicked, this, &KLOEMainWindow::onNeutRecClicked);
        connect(interfButton, &QPushButton::clicked, this, &KLOEMainWindow::onInterfClicked);
        connect(plotsButton, &QPushButton::clicked, this, &KLOEMainWindow::onPlotsClicked);
    }
};

/**
 * @brief Funkcja główna dla aplikacji Qt
 */
int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
    
    KLOEMainWindow window;
    window.show();
    
    return app.exec();
}

#include "desktop_app_qt.moc"  // Wymagane dla Q_OBJECT

#else

#include <iostream>

int main(int argc, char *argv[])
{
    std::cout << "=== Przykład aplikacji Qt ===" << std::endl;
    std::cout << "Ten przykład wymaga zainstalowanego Qt." << std::endl;
    std::cout << "Aby skompilować:" << std::endl;
    std::cout << "1. Zainstaluj Qt (sudo apt install qtbase5-dev)" << std::endl;
    std::cout << "2. Zdefiniuj QT_AVAILABLE przy kompilacji" << std::endl;
    std::cout << "3. Użyj qmake lub CMake z Qt" << std::endl;
    return 0;
}

#endif