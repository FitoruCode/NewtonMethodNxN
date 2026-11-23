#!/usr/bin/env bash
set -e
mkdir -p output
gcc PROJEKT.c libs/src/nrutil.c -Ilibs/include -I/usr/local/dislin -L/usr/local/dislin -Wl,-rpath,/usr/local/dislin -ldislin -lm -o projekt
./projekt
