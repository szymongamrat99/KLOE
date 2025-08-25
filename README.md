# KLOE Physics Analysis Project

Projekt analizy danych z eksperymentu KLOE (K LOng Experiment) - analiza fizyki wysokich energii z użyciem C++ i ROOT framework.

## 🚀 Nowa funkcjonalność: Aplikacje desktopowe

**Odpowiedź na pytanie: "Jak w C++ mógłbym stworzyć aplikację desktopową?"**

Zobacz [`examples/`](examples/) - zawiera kompletne przykłady tworzenia aplikacji desktopowych w C++:

- **🎮 [`desktop_app_demo.cpp`](examples/desktop_app_demo.cpp)** - Szybkie demo (uruchom najpierw!)
- **📚 [`desktop_app_basic.cpp`](examples/desktop_app_basic.cpp)** - Podstawowa struktura
- **🎨 [`desktop_app_qt.cpp`](examples/desktop_app_qt.cpp)** - Nowoczesne GUI z Qt
- **🔬 [`desktop_app_root.cpp`](examples/desktop_app_root.cpp)** - GUI zoptymalizowane dla fizyki

**Szybki start:**
```bash
cd examples
g++ -std=c++14 -o demo desktop_app_demo.cpp && ./demo
```

Pełny przewodnik: [`docs/DESKTOP_APP_GUIDE.md`](docs/DESKTOP_APP_GUIDE.md)

## 🔬 O projekcie KLOE

Projekt zawiera analizy fizyczne dla:
- Rekonstrukcji naładowanych kaonów (K→π+π-)
- Rekonstrukcji neutralnych kaonów (K→π0π0) 
- Analizy interferometrii
- Analizy naruszenia symetrii CP
- Analizy zmiennych generowanych Monte Carlo

## 📁 Struktura projektu

```
KLOE/
├── app/                    # Główna aplikacja konsolowa
├── examples/              # 🆕 Przykłady aplikacji desktopowych
├── docs/                  # 🆕 Dokumentacja GUI
├── Include/
│   ├── Codes/            # Główne biblioteki analizy
│   └── FortranAnalysis/  # Integracja z kodem FORTRAN
├── Subanalysis/          # Moduły specjalistycznych analiz
│   ├── KchRec/          # Rekonstrukcja naładowanych kaonów
│   ├── Neutrec/         # Rekonstrukcja neutralnych kaonów
│   ├── CPFit/           # Analiza naruszenia CP
│   └── ...
└── tests/               # Testy jednostkowe
```

## 🛠️ Kompilacja

**Podstawowa kompilacja (tylko przykłady GUI):**
```bash
cd examples
mkdir build && cd build
cmake ..
make
```

**Pełna kompilacja projektu (wymaga ROOT, Boost, LAPACK):**
```bash
mkdir build && cd build
cmake ..
make
```

## 📋 Wymagania

**Dla przykładów GUI:**
- C++14 lub nowszy
- CMake 3.13+
- Opcjonalnie: Qt5 (dla GUI)
- Opcjonalnie: ROOT framework (dla wizualizacji fizycznych)

**Dla pełnego projektu:**
- ROOT framework
- Boost libraries  
- LAPACK/BLAS
- FORTRAN compiler (gfortran)

## 🎯 Zastosowania

Ten projekt pokazuje jak:
- ✅ Stworzyć aplikację desktopową w C++ (różne podejścia)
- ✅ Zintegrować GUI z analizą naukową
- ✅ Używać ROOT framework do wizualizacji
- ✅ Łączyć C++ z kodem FORTRAN
- ✅ Organizować duży projekt naukowy

## 📚 Dokumentacja

- [`docs/DESKTOP_APP_GUIDE.md`](docs/DESKTOP_APP_GUIDE.md) - Przewodnik tworzenia aplikacji desktopowych
- [`examples/README.md`](examples/README.md) - Opis przykładów kodu
- `docs/html/` - Dokumentacja Doxygen (generowana automatycznie)

## 👨‍💻 Autor

**Szymon Gamrat**
- GitHub: [@szymongamrat99](https://github.com/szymongamrat99)
- Projekt: Analiza danych eksperymentu KLOE

## 📄 Licencja

Projekt edukacyjny/naukowy - zobacz konkretne pliki dla szczegółów licencji.

---

**🎯 Dla osób uczących się tworzenia aplikacji desktopowych w C++:**
Rozpocznij od katalogu [`examples/`](examples/) - zawiera praktyczne przykłady z pełnymi instrukcjami!