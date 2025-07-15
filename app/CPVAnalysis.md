# CPVAnalysis.cpp – Parametry wejściowe i logika uruchamiania

## Opis pliku

Plik `CPVAnalysis.cpp` to główny plik wykonywalny analizy KLOE. Odpowiada za:
- Pobranie i walidację parametrów wejściowych od użytkownika (zakres plików, tryb analizy, typ danych),
- Inicjalizację łańcucha plików ROOT (`TChain`),
- Uruchomienie głównej pętli analizy lub trybu inicjalizacyjnego,
- Obsługę logowania błędów i informacji.

## Parametry wejściowe

### 1. Tryb analizy
- **Measurement (1)** – analiza właściwa (dane sygnałowe)
- **Control Sample (2)** – analiza próbki kontrolnej
- Wybierany przez użytkownika na początku programu.

### 2. Zakres plików
- **firstFile** – numer pierwszego pliku do analizy (liczba całkowita, >= 1)
- **lastFile** – numer ostatniego pliku do analizy (liczba całkowita, >= firstFile)
- Maksymalny numer pliku jest pobierany z pliku konfiguracyjnego (np. `properties["variables"]["rootFiles"]["lastFileMax"]`).
- Program sprawdza poprawność zakresu i wymusza podanie prawidłowych wartości.

### 3. Typ danych (dataTypeOpt)
- Wybierany z menu (np. sygnał, tło, MC, itp.)
- Wartość jest walidowana – musi należeć do dozwolonego zakresu enum `Controls::DataType`.

### 4. Flaga `initialAnalysisExec`
- Jeśli `true`, program uruchamia tryb inicjalizacyjny (np. tylko zaczytuje pliki i wykonuje InitAnalysis_main).
- Jeśli `false`, uruchamiana jest pełna pętla menu analizy.
- Wartość pobierana z pliku konfiguracyjnego `properties.json`.

## Przebieg programu

1. **Pobranie parametrów**
   - Jeśli `initialAnalysisExec == false`, wywoływana jest funkcja `InputParamsHandler::getParams`, która:
     - Pyta użytkownika o tryb (measurement/control), zakres plików i typ danych,
     - Waliduje dane wejściowe,
     - Zapisuje wybrane parametry do pliku konfiguracyjnego.

2. **Inicjalizacja TChain**
   - Na podstawie wybranych parametrów, do łańcucha `TChain` dodawane są odpowiednie pliki.
   - W trybie inicjalizacyjnym (`initialAnalysisExec == true`) zakres runów i plików jest ustalany automatycznie na podstawie zawartości katalogu i regexa.

3. **Uruchomienie analizy**
   - W trybie pełnym: uruchamiana jest pętla menu (`mainMenuHandler.runMenuLoop`).
   - W trybie inicjalizacyjnym: wywoływana jest funkcja `InitAnalysis_main`.

4. **Logowanie**
   - Wszystkie błędy i ważne informacje są logowane przez obiekt `ErrorHandling::ErrorLogs`.

## Przykład uruchomienia

```
$ ./CPVAnalysis
Measurement (1) / Control Sample (2)? 1
Choose the first file: 10
Choose the last file: 20
[menu wyboru typu danych]
...
```

## Pliki konfiguracyjne
- Parametry wejściowe i zakresy plików są zapisywane do pliku `properties.json`.
- Struktura pliku konfiguracyjnego pozwala na łatwe odtworzenie analizy z tymi samymi parametrami.

## Uwagi
- Program wymaga poprawnego pliku konfiguracyjnego oraz obecności plików ROOT w zadanych lokalizacjach.
- Wszelkie błędy wejścia są obsługiwane i logowane, a użytkownik proszony jest o ponowne podanie wartości.

---

**Autor:** automatyczny asystent refaktoryzacji
**Data:** 2025-07-09
