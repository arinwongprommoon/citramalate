#!/usr/bin/env python

# Replaces elements in all numpy arrays in a specific directory lower than some
# threshold value with -1e-4

import glob, os
import numpy as np
import matplotlib.pyplot as plt
import heatmap

for npz in glob.glob("*.npz"):
    # Read from file into numpy arrays
    container = np.load(npz)
    data = [container[key] for key in container]

    # remove the temp bit when I've proven that it works well
    temp = data[2]
    temp[temp > 1e-8] = -1e-4
    data[2] = temp

    newdata = [data[0], data[1], data[2]]
    np.savez(npz, *newdata)
