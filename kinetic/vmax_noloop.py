#!/usr/bin/env python

# Runs computations for reactions within a specified list, writes data to TXT
# file, then plots graphs of productivity against multiple of Vmax and saves
# them as image files

# This program stops after each reaction and continues on to
# the next reaction when executed again. To run through all reaction, execute
# ./vmax_run.sh . In a nutshell, this script is not very useful on its own

from __future__ import division, print_function
from ecolicitra import ecolicit, mmCITRA, mmGLC
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import numpy as np
import sys
import csv

# Parameters of kinetic model
include_CITRA = True
ecit = ecolicit(include_CITRA = include_CITRA)
ecit.setVmax('CITRA_SYN', 4.0)
ecit.time0 = 0
ecit.timef = 2*3600 # final simulation time in seconds
ecit.npoints = 100 # Number of points to be computed in the simulation
ecit.getVmaxes()
listofreactions = ecit.reacVmaxes

start = 0.1 # START VALUE IN MULTIPLES OF WT VMAX
end = 5.0 # END VALUE IN MULTIPLES OF WT VMAX
points = 50 # NUMBER OF DATA POINTS

yaxisfixed = False # IF YOU WANT TO FIX THE Y-AXIS
y = (0.0070, 0.23) # DEFINE Y-AXIS LIMITS HERE
yl = 'flux' # DEFINE Y-AXIS LABEL HERE

xFormatter = FormatStrFormatter('%.2f')

# read the index of the reaction from file
with open('vmaxi.txt', 'r') as fobj:
    mystr = fobj.readline()
    idx = int(mystr.strip('\n'))

# stops if all reactions gone through
if idx >= len(listofreactions):
    print('DONE')
    sys.exit(1)

else:
    reaction = listofreactions[idx]

    print(reaction)
    initVmax = ecit.getVmax(reaction)

    # Plot graph and saves image
    X = np.linspace(start*initVmax, end*initVmax, points, endpoint=True) # sets X values from 0.5*default Vmax to 2.0*default Vmax
    P = []
    for i in range(points):
        ecit.setVmax(reaction, X[i])
        P.append(ecit.comflux())

    ecit.setVmax(reaction, initVmax)

    # write the index of the reaction to file
    idx += 1
    with open('vmaxi.txt', 'w') as fobj:
        fobj.write(str(idx))

    # Write data to file
    data = [str(reaction), str(X.tolist()), str(P)]
    with open("VMAXDATA.csv", 'a') as csvfile:
        datawriter = csv.writer(csvfile)
        datawriter.writerow([str(reaction), str(min(P)), str(max(P))])

    sys.exit(2)
