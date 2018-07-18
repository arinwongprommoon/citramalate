#!/usr/bin/env python

# Replaces elements in all numpy arrays in a specific directory lower than some
# threshold value with -1e-4

import glob, os
import numpy as np
import matplotlib.pyplot as plt
import heatmap

for cnpz in glob.glob("COUPLES*.npz"):
    # Read from file into numpy arrays

    # Reads the COUPLES*.npz file
    c_container = np.load(cnpz)

    # Reads the corresponding STEADY*.npz file
    snpz = cnpz.replace("COUPLES", "STEADY")
    s_container = np.load(snpz)
    cdata = [c_container[key] for key in c_container]
    sdata = [s_container[key] for key in s_container]

    # if element in STEADY exceeds epsilon, replace corresponding element in
    # COUPLES with -1e-4
    ctemp = cdata[2]
    stemp = sdata[2]
    ctemp[stemp > 1e-6] = -1e-4
    cdata[2] = ctemp

    newcdata = [cdata[0], cdata[1], cdata[2]]
    np.savez(cnpz, *newcdata)
