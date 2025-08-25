# Jak stworzyć aplikację desktopową w C++

## Przewodnik po tworzeniu aplikacji desktopowych dla analizy fizycznej KLOE

### Wprowadzenie

W projekcie KLOE obecnie używamy aplikacji konsolowej do analizy danych fizycznych. Ten przewodnik pokazuje, jak można rozszerzyć projekt o interfejs graficzny (GUI), aby stworzyć nowoczesną aplikację desktopową.

### Dostępne opcje dla GUI w C++

1. **Qt** - najpopularniejsza, cross-platform
   - Profesjonalny wygląd
   - Bogate komponenty GUI
   - Dobra integracja z C++
   - Wsparcie dla wykresów i wizualizacji

2. **wxWidgets** - natywny wygląd na każdej platformie
   - Używa natywnych kontrolek systemu
   - Mniejszy rozmiar aplikacji
   - Cross-platform

3. **GTK** - głównie dla Linuxa
   - Używany w środowisku GNOME
   - Dobra integracja z systemami Unix

4. **ROOT GUI** - specjalnie dla aplikacji fizycznych
   - Już używamy ROOT w projekcie
   - Zoptymalizowane dla wizualizacji danych fizycznych
   - TCanvas, TBrowser, TGWindow

5. **Dear ImGui** - immediate mode GUI
   - Doskonałe dla narzędzi deweloperskich
   - Szybkie prototypowanie
   - Lekkie i wydajne

### Rekomendacja dla projektu KLOE

Ze względu na to, że już używamy ROOT framework, najlepszym rozwiązaniem będzie:
1. **ROOT GUI** - dla podstawowych okien i kontrolek
2. **Qt** - dla bardziej zaawansowanego interfejsu

### Przykłady implementacji

W katalogu `examples/` znajdziesz konkretne przykłady kodu dla każdego podejścia:

#### 1. Podstawowa aplikacja desktopowa (`desktop_app_basic.cpp`)

```cpp
// Przykład bazowej struktury aplikacji desktopowej
class KLOEDesktopApp {
    // Podstawowa logika menu i obsługi analiz
};
```

**Kompilacja:**
```bash
cd examples
mkdir build && cd build
cmake ..
make desktop_app_basic
./desktop_app_basic
```

#### 2. Aplikacja Qt (`desktop_app_qt.cpp`)

```cpp
// Przykład z pełnym GUI Qt
class KLOEMainWindow : public QMainWindow {
    Q_OBJECT
    // Przyciski, menu, paski postępu
};
```

**Wymagania:**
- Qt5 lub nowszy
- CMake z find_package(Qt5)

**Instalacja Qt na Ubuntu/Debian:**
```bash
sudo apt update
sudo apt install qtbase5-dev qtchooser qt5-qmake qtbase5-dev-tools
```

**Kompilacja:**
```bash
cd examples/build
cmake .. -DQt5_DIR=/usr/lib/x86_64-linux-gnu/cmake/Qt5
make desktop_app_qt
./desktop_app_qt
```

#### 3. Aplikacja ROOT GUI (`desktop_app_root.cpp`)

```cpp
// Przykład używający ROOT TGWindow
class KLOERootMainFrame : public TGMainFrame {
    // Zintegrowane z ROOT Canvas i analizą fizyczną
};
```

**Wymagania:**
- ROOT framework
- CMake z find_package(ROOT)

**Instalacja ROOT:**
```bash
# Na Ubuntu/Debian
sudo apt install root-system

# Lub kompilacja ze źródeł
wget https://root.cern/download/root_v6.XX.XX.source.tar.gz
# Następnie kompilacja według instrukcji ROOT
```

**Kompilacja:**
```bash
cd examples/build
source /usr/bin/thisroot.sh  # lub ścieżka do ROOT
cmake ..
make desktop_app_root
./desktop_app_root
```

### Integracja z istniejącym projektem KLOE

Aby zintegrować GUI z istniejącym kodem analizy:

1. **Dodaj podkatalog examples do głównego CMakeLists.txt:**
```cmake
add_subdirectory(examples)
```

2. **Połącz z istniejącymi bibliotekami:**
```cpp
// W funkcjach GUI wywołuj istniejące analizy
void RunAnalysis() {
    GenVars_main(chain, eventAnalysis, dataTypeOpt, physConst);
    KchRec_main(chain, eventAnalysis, dataTypeOpt, physConst);
    // itd.
}
```

3. **Użyj istniejących klas:**
```cpp
#include "../Include/klspm00.hpp"

// W konstruktorze GUI
KLOE::pm00 eventAnalysis;
Controls::Menu mainMenu(10);
```

### Zalety każdego podejścia

**Podstawowa aplikacja:**
- ✅ Proste w implementacji
- ✅ Bez dodatkowych zależności
- ❌ Tylko interfejs tekstowy

**Qt:**
- ✅ Profesjonalny wygląd
- ✅ Cross-platform
- ✅ Bogate komponenty GUI
- ✅ Dobra dokumentacja
- ❌ Dodatkowa zależność

**ROOT GUI:**
- ✅ Już używasz ROOT w projekcie
- ✅ Zoptymalizowane dla fizyki
- ✅ Bezpośrednia integracja z TCanvas
- ✅ Wykresy i histogramy out-of-the-box
- ❌ Przestarzały wygląd
- ❌ Mniej uniwersalne

### Następne kroki

1. Wybierz podejście odpowiednie dla twoich potrzeb
2. Zainstaluj wymagane biblioteki
3. Skompiluj przykłady
4. Zmodyfikuj kod pod swoje wymagania
5. Zintegruj z istniejącymi analizami KLOE