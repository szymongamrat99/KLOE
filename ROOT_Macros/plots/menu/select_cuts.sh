#!/bin/bash

CUTS_FILE="${1:-cut-limits.json}"
HYPOTHESIS="${2:-SIMONA_ANALYSIS}"

PROPERTIES_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../Subanalysis" && pwd)/Properties"
CUTS_PATH="$PROPERTIES_DIR/$CUTS_FILE"
RESULT_FILE="/tmp/selected_cuts_result.txt"

if [[ ! -f "$CUTS_PATH" ]]; then
    echo "Błąd: Plik $CUTS_PATH nie istnieje" >&2
    exit 1
fi

# Parsuj JSON i wyodrębnij cięcia dla danej hipotezy
python3 - "$CUTS_PATH" "$HYPOTHESIS" << 'PYTHON_EOF'
import json
import sys

cuts_file = sys.argv[1]
hypothesis = sys.argv[2]

with open(cuts_file, 'r') as f:
    data = json.load(f)

cuts = data.get("listOfCuts", {}).get(hypothesis, [])

if not cuts:
    print("Błąd: Brak cięć dla hipotezy: " + hypothesis, file=sys.stderr)
    sys.exit(1)

# Sortuj po order
cuts = sorted(cuts, key=lambda x: x.get("order", 999))

# Wydrukuj cięcia
for i, cut in enumerate(cuts, 1):
    fv_mark = "[FV]" if cut.get("isFiducialVolume", False) else "    "
    cut_id = cut.get("cutId", "N/A")
    desc = cut.get("cutDescription", "")
    print(f"{i:2d}. {fv_mark} {cut_id:20s} - {desc}")

print("---", file=sys.stderr)
PYTHON_EOF

# Teraz menu interaktywne
selected=()
exec 3<&0  # Zapisz oryginalny stdin

while true; do
    clear >&2
    
    # Wydrukuj listę cięć z tickmarkami
    python3 - "$CUTS_PATH" "$HYPOTHESIS" "${selected[@]}" << 'PYTHON_DISPLAY'
import json
import sys

cuts_file = sys.argv[1]
hypothesis = sys.argv[2]
selected_list = sys.argv[3:] if len(sys.argv) > 3 else []

with open(cuts_file, 'r') as f:
    data = json.load(f)

cuts = data.get("listOfCuts", {}).get(hypothesis, [])
cuts = sorted(cuts, key=lambda x: x.get("order", 999))

for i, cut in enumerate(cuts, 1):
    fv_mark = "[FV]" if cut.get("isFiducialVolume", False) else "    "
    cut_id = cut.get("cutId", "N/A")
    desc = cut.get("cutDescription", "")
    tick = "✓" if cut_id in selected_list else " "
    print(f"[{tick}] {i:2d}. {fv_mark} {cut_id:20s} - {desc}")

print("---")
PYTHON_DISPLAY
    
    echo "" >&2
    echo "Wybrane cięcia (${#selected[@]}):" >&2
    if [[ ${#selected[@]} -eq 0 ]]; then
        echo "  (brak)" >&2
    else
        for s in "${selected[@]}"; do
            echo "  ✓ $s" >&2
        done
    fi
    
    echo "" >&2
    echo "Opcje:" >&2
    echo "  [numer] - wybierz/usuń cięcie" >&2
    echo "  [d]     - dokończ" >&2
    echo "  [c]     - wyczyść" >&2
    echo "" >&2
    echo -n "Wybór: " >&2
    
    read -u 3 choice
    
    if [[ "$choice" == "d" ]]; then
        break
    elif [[ "$choice" == "c" ]]; then
        selected=()
    elif [[ "$choice" =~ ^[0-9]+$ ]]; then
        idx=$((choice - 1))
        cut_name=$(python3 - "$CUTS_PATH" "$HYPOTHESIS" "$idx" << 'PYTHON_GET'
import json
import sys

cuts_file = sys.argv[1]
hypothesis = sys.argv[2]
idx = int(sys.argv[3])

with open(cuts_file, 'r') as f:
    data = json.load(f)

cuts = data.get("listOfCuts", {}).get(hypothesis, [])
cuts = sorted(cuts, key=lambda x: x.get("order", 999))

if 0 <= idx < len(cuts):
    print(cuts[idx].get("cutId", ""))
PYTHON_GET
)
        
        if [[ -n "$cut_name" ]]; then
            # Sprawdź czy już jest zaznaczone
            found=0
            for i in "${!selected[@]}"; do
                if [[ "${selected[$i]}" == "$cut_name" ]]; then
                    # Usuń to cięcie
                    unset 'selected[$i]'
                    found=1
                    break
                fi
            done
            
            if [[ $found -eq 0 ]]; then
                # Dodaj nowe cięcie
                selected+=("$cut_name")
            fi
            
            # Przeindeksuj tablicę (usuwa puste slots)
            selected=("${selected[@]}")
        fi
    fi
done

exec 3<&-  # Zamknij zapisany stdin

# Sortuj: FV na początek, potem reszta
fv_cuts=()
other_cuts=()

for cut in "${selected[@]}"; do
    is_fv=$(python3 - "$CUTS_PATH" "$HYPOTHESIS" "$cut" << 'PYTHON_FV'
import json
import sys

cuts_file = sys.argv[1]
hypothesis = sys.argv[2]
cut_id = sys.argv[3]

with open(cuts_file, 'r') as f:
    data = json.load(f)

cuts = data.get("listOfCuts", {}).get(hypothesis, [])
for c in cuts:
    if c.get("cutId") == cut_id:
        print("1" if c.get("isFiducialVolume", False) else "0")
        break
PYTHON_FV
)
    
    if [[ "$is_fv" == "1" ]]; then
        fv_cuts+=("$cut")
    else
        other_cuts+=("$cut")
    fi
done

# Łącz: FV na początek
final_cuts=("${fv_cuts[@]}" "${other_cuts[@]}")
result=$(IFS=+; echo "${final_cuts[*]}")

echo "$result" > "$RESULT_FILE"