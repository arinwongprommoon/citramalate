#!/usr/bin/env python

# Runs computations for couples of reactions from a specified list, writes data
# to numpy array files, then plots heatmaps.

import glob, os
import numpy as np
import matplotlib.pyplot as plt
import couples
import heatmap

#Run simulations, writes data both into TXT and NPZ
run = couples.Coupl(xstart=0.5, xend=2.0, ystart=0.5, yend=2.0, points=10) # define vmax ranges here
run.setVmax('CITRA_SYN', 4.0)
run.time0 = 0
run.timef = 2*3600 # final simulation time in seconds
run.npoints = 100 # Number of points to be computed in the simulation
run.listofreactions = ['CITRA_SYN', 'GLT', 'LPD', 'GDH', 'ATP_syn'] # define list of reactions
#run.listofreactions = ['CITRA_SYN', 'GLT', 'LPD', 'ATP_MAINTENANCE', 'GDH', 'ATP_syn', 'ACEA', 'ZWF']

run.writeFileHeader()
run.computeCouples()

for npz in glob.glob("*.npz"):
    # Read from file into numpy arrays
    container = np.load(npz)
    data = [container[key] for key in container]

    # Creates the heatmaps
    XRxns = data[0]
    YRxns = data[1]
    productivity_raw = data[2]
    productivity = 10000*productivity_raw
    # Transpose here because I don't want to mess with couples.py, which
    # generates the data
    productivity = np.transpose(productivity)

    XLabs = ['%.2f' % i for i in XRxns]
    YLabs = ['%.2f' % j for j in YRxns]

    fig, ax = plt.subplots()
    im, cbar = heatmap.heatmap(productivity, XLabs, YLabs, ax=ax,
                       cmap="YlGn", cbarlabel="citramalate productivity [10,000 * raw]")
    texts = heatmap.annotate_heatmap(im, valfmt="{x:.1f}")

    fig.tight_layout()

    filename = npz[0:-4] + '.png'
    plt.savefig(filename, bbox_inches='tight')
    print(filename + ' saved')
    plt.gcf().clear()
    plt.close()
