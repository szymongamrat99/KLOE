import os

def parse_histogram_configs(filename):
    """
    Parsuje plik tekstowy dla histogramow 1D i 2D.
    Zwraca slownik: { nazwa_kolumn : (model, lista_kolumn) }
    """
    configs = {}

    if not os.path.exists(filename):
        print("Blad: Plik {} nie istnieje!".format(filename))
        return configs

    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            parts = [p.strip() for p in line.split(';')]
            
            # --- Obsluga Histogramow 2D (Twoj nowy format: ok. 11 pol) ---
            if len(parts) >= 11:
                name = parts[0]
                # Kolumny dla 2D sa zazwyczaj rozdzielone przez '_vs_' w nazwie
                # Zakladamy, ze nazwa w pliku to "kolX_vs_kolY"
                cols = name.split('_vs_')
                
                nbinsX = int(parts[1])
                nbinsY = int(parts[2])
                xmin = float(parts[3])
                xmax = float(parts[4])
                ymin = float(parts[5])
                ymax = float(parts[6])
                
                t_main = parts[7].replace('"', '')
                t_x = parts[8].replace('"', '')
                t_y = parts[9].replace('"', '')
                t_z = parts[10].replace('"', '')

                full_title = "{};{};{};{}".format(t_main, t_x, t_y, t_z)
                
                # Model dla Histo2D: (nazwa, tytul, nX, xmin, xmax, nY, ymin, ymax)
                model = (name, full_title, nbinsX, xmin, xmax, nbinsY, ymin, ymax)
                configs[name] = {"type": "2D", "model": model, "cols": cols}

            # --- Obsluga Histogramow 1D (Poprzedni format) ---
            elif len(parts) >= 7:
                name = parts[0]
                nbins = int(parts[1])
                xmin = float(parts[2])
                xmax = float(parts[3])
                
                t_main = parts[4].replace('"', '')
                t_x = parts[5].replace('"', '')
                t_y = parts[6].replace('"', '')

                full_title = "{};{};{}".format(t_main, t_x, t_y)
                
                model = (name, full_title, nbins, xmin, xmax)
                configs[name] = {"type": "1D", "model": model, "cols": [name]}

    return configs