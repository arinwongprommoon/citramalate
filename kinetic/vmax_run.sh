#!/bin/bash
# Runs computations for reactions within a specified list, writes data to TXT
# file, then plots graphs of productivity against multiple of Vmax and saves
# them as image files

# This gets around the memory leak

echo "0" > 'vmaxi.txt'

while true; do
    python vmax_noloop.py $1
    a=$?
    if [ $a == 1 ]; then
      break
    fi
done
