#!/bin/bash
# Changes vmax values around and investigates the effect on fluxes of all
# reactions in the model. Outputs to csv files

# This gets around the memory leak

echo "0" > 'fluxi.txt'

while true; do
    python flux_noloop.py $1 $2
    a=$?
    if [ $a == 1 ]; then
      break
    fi
done

python combine.py
python getminmax.py
