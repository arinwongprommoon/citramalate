#!/bin/bash
# Runs computations for couples of reactions from a specified list, writes data
# to numpy array files. Run couples_plot.py to generate heatmaps after.

# This gets around the memory leak

echo "0" > 'couplesi.txt'

while true; do
    python couples_compute.py
    a=$?
    if [ $a == 1 ]; then
	break
    fi
done
