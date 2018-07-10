#!/usr/bin/env python

# Runs computations for couples of reactions from a specified list, writes data
# to numpy array files, then plots heatmaps.

import glob, os
import numpy as np
import matplotlib.pyplot as plt
import couples_noloop
import heatmap

#Run simulations, writes data both into TXT and NPZ
run = couples.Coupl(xstart=0.1, xend=1.0, ystart=0.1, yend=1.0, points=10) # define vmax ranges here
run.setVmax('CITRA_SYN', 4.0)
run.time0 = 0
run.timef = 2*3600 # final simulation time in seconds
run.npoints = 100 # Number of points to be computed in the simulation
run.listofreactions = ['CITRA_SYN', 'GLT', 'LPD', 'GDH','ATP_syn'] # define list of reactions
#run.listofreactions = ['CITRA_SYN', 'GLT', 'LPD', 'ATP_MAINTENANCE', 'GDH', 'ATP_syn', 'ACEA', 'ZWF']

run.writeFileHeader()

# with open('couplesi.txt', 'w') as fobj:
#     fobj.write('0')

run.computeCouples()
