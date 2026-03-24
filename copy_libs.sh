#!/bin/bash
cp bin/KLSPM00 /opt/exp_software/kloe/users/gamrat/KLSPM00/bin/KLSPM00.exe
find . -name "*.so" -exec cp {} /opt/exp_software/kloe/users/gamrat/KLSPM00/lib/. \;
