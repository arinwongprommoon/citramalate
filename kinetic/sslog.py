#!/usr/bin/env python

# Replaces element x in STEADY*.npz with -log(x)

import glob, os
import numpy as np
import matplotlib.pyplot as plt

for npz in glob.glob("STEADY*.npz"):
    # Read from file into numpy arrays
    container = np.load(npz)
    data = [container[key] for key in container]
    temp = -np.log(data[2])
    newdata = [data[0], data[1], temp]
    np.savez(npz, *newdata)
