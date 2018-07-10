#!/bin/bash
# Runs computations for couples of reactions from a specified list, writes data
# to numpy array files, then plots heatmaps.

# This gets around the memory leak

echo "0" > 'couplesi.txt'

while true; do
      python couples_compute.py
done
