#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;

// Poprawna składnia dla ROOT 6 (wymaga C++)
#pragma link C++ namespace KLOE;
#pragma link C++ namespace Paths;
#pragma link C++ namespace Filenames;
#pragma link C++ namespace PhysicsConstants;
#pragma link C++ namespace Utils;
#pragma link C++ namespace ErrorHandling;

// Jeśli masz konkretne klasy, które chcesz wywoływać bezpośrednio, dopisz je w ten sposób:
// #pragma link C++ class KLOE::Histograms::HistConfig1D+;

#endif