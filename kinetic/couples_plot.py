#!/usr/bin/env python

# Plots heatmaps from all numpy arrays (*.npz) saved in current directory

import glob, os
import numpy as np
import matplotlib.pyplot as plt
import heatmap

for npz in glob.glob("*.npz"):
    # Read from file into numpy arrays
    container = np.load(npz)
    data = [container[key] for key in container]

    # Creates the heatmaps
    
    XRxns = data[0]
    YRxns = data[1]
    productivity_raw = data[2]
    # multiplies by 10,000 to get numbers that can be easily visualised by eye
    productivity = 10000*productivity_raw
    productivity = np.transpose(productivity)

    XLabs = ['%.2f' % i for i in XRxns]
    YLabs = ['%.2f' % j for j in YRxns]

    fig, ax = plt.subplots()
    im, cbar = heatmap.heatmap(productivity, XLabs, YLabs, ax=ax,
                       cmap="YlGn", cbarlabel="citramalate productivity [10,000 * raw]")
    texts = heatmap.annotate_heatmap(im, valfmt="{x:.1f}", threshold=0, textcolors=["red", "black"])

    fig.tight_layout()

    filename = npz[0:-4] + '.png'
    plt.savefig(filename, bbox_inches='tight')
    print(filename + ' saved')
    plt.gcf().clear()
    plt.close()
