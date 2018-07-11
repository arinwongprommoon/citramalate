#!/bin/bash
# Runs differential evolution for tuples of reactions from a specified list,
# writes data to a txt file.

# This gets around the memory leak

echo "0" > 'dei.txt'

while true; do
    python de_noloop.py
    a=$?
    if [ $a == 1 ]; then
      break
    fi
done
