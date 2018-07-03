#!/usr/bin/env python

from __future__ import division, print_function
from ecolicitra_copy import ecolicit, mmCITRA, mmGLC
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import numpy as np
import roadrunner
import libsbml

# Parameters of kinetic model
include_CITRA = True
ecit = ecolicit(include_CITRA = include_CITRA)
ecit.setVmax('CITRA_SYN', 4.0)

ecit.time0 = 0
ecit.timef = 2*3600 # final simulation time in seconds
ecit.npoints = 100 # Number of points to be computed in the simulation

wildprod = ecit.comproducti()
print("Final simulation time: ", ecit.timef)
print("Number of points: ", ecit.npoints)

xFormatter = FormatStrFormatter('%.1f')

# Sets X values from start*default Vmax to 5.0*default Vmax
# Plot uses 'points' number of points e.g. 20
start = 0.5
end = 2.0
points = 25

# Clears data in file, writes header
with open("VMAXDATA.txt", 'w') as fobj:
    fobj.write('time0 ' + str(ecit.time0) + '\n')
    fobj.write('timef ' + str(ecit.timef) + '\n')
    fobj.write('Vmax ' + str(start) + 'to' + str(end) + '\n\n')

ecit.getVmaxes()
listofreactions = ecit.reacVmaxes

for reaction in listofreactions:
    print(reaction)
    initVmax = ecit.getVmax(reaction)

    # Plot graph and saves image
    X = np.linspace(start*initVmax, end*initVmax, points, endpoint=True) # sets X values from 0.5*default Vmax to 2.0*default Vmax
    P = []
    for i in range(points):
        ecit.setVmax(reaction, X[i])
        P.append(ecit.comproducti())
    # redefines X so that I can get a plot of productivity against multiples of
    # Vmax because all the methods I used to get rid of floating point
    # representation errors while plotting do not work. X still retained for
    # writing data into file
    XX = np.linspace(start, end, points, endpoint=True)

    ax = plt.subplot(111)
    plt.plot(XX,P)
    #plt.ylim(0.0007, 0.0017)
    plt.ylim(0.0002, 0.0051)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    ax.xaxis.set_major_formatter(xFormatter)

    filename = 'ONE_' + reaction + '.png'
    plt.savefig(filename, bbox_inches='tight')
    plt.gcf().clear()

    ecit.setVmax(reaction, initVmax)

    # Write data to file
    data = [str(reaction), str(X.tolist()), str(P)]
    with open("VMAXDATA.txt", 'a') as fobj:
        fobj.write('\n'.join(data))
        fobj.write('\n')

# To do: make this faster (takes about 4 minutes)
