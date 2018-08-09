#!/bin/bash
# Changes vmax values around and investigates the effect on fluxes of all
# reactions in the model. Outputs to csv files

# This gets around the memory leak

echo "0" > 'jmin.txt'

while true; do
    python fluxde_min.py
    a=$?
    if [ $a == 1 ]; then
      break
    fi
done
