# Przykłady aplikacji desktopowych w C++

Ten katalog zawiera praktyczne przykłady różnych podejść do tworzenia aplikacji desktopowych w C++ dla projektu analizy fizycznej KLOE.

## 📁 Dostępne przykłady

### 1. `desktop_app_demo.cpp` 
**🚀 Szybki start - uruchom to najpierw!**
```bash
g++ -std=c++14 -o desktop_app_demo desktop_app_demo.cpp
./desktop_app_demo
```
Automatyczne demo pokazujące wszystkie funkcjonalności bez interakcji użytkownika.

### 2. `desktop_app_basic.cpp`
**📚 Podstawowa struktura aplikacji**
```bash
g++ -std=c++14 -o desktop_app_basic desktop_app_basic.cpp
./desktop_app_basic
```
Interaktywne menu tekstowe pokazujące architekturę aplikacji desktopowej.

### 3. `desktop_app_qt.cpp`
**🎨 Nowoczesne GUI z Qt**
```bash
# Wymaga zainstalowanego Qt5
sudo apt install qtbase5-dev
g++ -std=c++14 -DQT_AVAILABLE desktop_app_qt.cpp -o desktop_app_qt $(pkg-config --cflags --libs Qt5Widgets)
./desktop_app_qt
```
Pełna aplikacja GUI z przyciskami, menu i paskami postępu.

### 4. `desktop_app_root.cpp`
**🔬 GUI zoptymalizowane dla fizyki**
```bash
# Wymaga zainstalowanego ROOT
source $ROOTSYS/bin/thisroot.sh
g++ -std=c++14 -DROOT_AVAILABLE desktop_app_root.cpp -o desktop_app_root $(root-config --cflags --libs --glibs)
./desktop_app_root
```
GUI używające ROOT framework z integracją TCanvas i wykresów fizycznych.

## 🛠️ Kompilacja za pomocą CMake

```bash
mkdir build && cd build
cmake ..
make

# Uruchom przykłady
./desktop_app_demo     # Zawsze działa
./desktop_app_basic    # Zawsze działa  
./desktop_app_qt       # Tylko jeśli Qt dostępne
./desktop_app_root     # Tylko jeśli ROOT dostępny
```

## 📖 Dalsze informacje

Zobacz pełny przewodnik: [`../docs/DESKTOP_APP_GUIDE.md`](../docs/DESKTOP_APP_GUIDE.md)

## 🎯 Odpowiedź na pytanie: "Jak w C++ mógłbym stworzyć aplikację desktopową?"

Te przykłady pokazują **cztery różne sposoby** tworzenia aplikacji desktopowych w C++:

1. **Aplikacja konsolowa z menu** - najprostsza, bez dodatkowych bibliotek
2. **Qt GUI** - najbardziej popularne, profesjonalny wygląd
3. **ROOT GUI** - specjalnie dla aplikacji naukowych/fizycznych
4. **Hybrydowe podejście** - łączenie różnych technologii

Każdy przykład zawiera:
- ✅ Kompletny kod źródłowy
- ✅ Instrukcje kompilacji
- ✅ Komentarze w języku polskim
- ✅ Integrację z istniejącym kodem analizy KLOE

**Polecana ścieżka nauki:**
1. Uruchom `desktop_app_demo.cpp` aby zobaczyć demo
2. Przeanalizuj `desktop_app_basic.cpp` aby zrozumieć strukturę
3. Wybierz Qt lub ROOT w zależności od potrzeb
4. Zmodyfikuj przykłady pod swoje wymagania