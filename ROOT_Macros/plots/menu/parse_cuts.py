#!/usr/bin/env python3
"""Parse cuts from JSON file and output in pipe-delimited format."""

import json
import sys

def main():
    if len(sys.argv) < 3:
        print("Usage: parse_cuts.py <cuts_file> <hypothesis>", file=sys.stderr)
        sys.exit(1)
    
    cuts_path = sys.argv[1]
    hypothesis = sys.argv[2]
    
    try:
        with open(cuts_path, 'r') as f:
            data = json.load(f)
    except Exception as e:
        print(f"Błąd: Nie mogę przeczytać pliku {cuts_path}: {e}", file=sys.stderr)
        sys.exit(1)
    
    if hypothesis not in data.get('listOfCuts', {}):
        print(f"Błąd: Nie znaleziono hipotezy '{hypothesis}' w pliku {cuts_path}", file=sys.stderr)
        print("Dostępne hipotezy:", file=sys.stderr)
        for h in data.get('listOfCuts', {}).keys():
            print(f"  - {h}", file=sys.stderr)
        sys.exit(1)
    
    cuts = data['listOfCuts'][hypothesis]
    
    # Sortuj po order
    cuts = sorted(cuts, key=lambda x: x.get('order', 0))
    
    for i, cut in enumerate(cuts):
        cut_id = cut.get('cutId', '')
        cut_desc = cut.get('cutDescription', '')
        is_fv = cut.get('isFiducialVolume', False)
        cut_order = cut.get('order', 0)
        
        # Wypisz dane w formacie do pobrania przez bash
        print(f"{i}|{cut_id}|{cut_desc}|{str(is_fv).lower()}|{cut_order}")

if __name__ == '__main__':
    main()
