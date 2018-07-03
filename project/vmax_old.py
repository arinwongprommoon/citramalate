#!/usr/bin/env python

from __future__ import division, print_function
from ecolicitra_copy import ecolicit, mmCITRA, mmGLC
import matplotlib.pyplot as plt
import numpy as np
import roadrunner
import libsbml

# Parameters of kinetic model
include_CITRA = True
ecit = ecolicit(include_CITRA = include_CITRA)
ecit.setVmax('CITRA_SYN', 4.0)

ecit.time0 = 0
ecit.timef = 2*3600 # final simulation time in seconds
ecit.npoints = 2 # Number of points to be computed in the simulation

wildprod = ecit.comproducti()
print("Final simulation time: ", ecit.timef)
print("Number of points: ", ecit.npoints)

ecit.getVmaxes()
listofreactions = ecit.reacVmaxes

# Clears data in file
with open("20180627TEST.txt", 'w') as fobj:
    pass

pltrange = np.arange(0.5, 2.0, 0.5)
ticks = list(pltrange)

for reaction in listofreactions:
    print(reaction)
    initVmax = ecit.getVmax(reaction)

    # Plot graph and saves image
    X = np.linspace(0.5*initVmax, 2.0*initVmax, 20, endpoint=True) # sets X values from 0.5*default Vmax to 2.0*default Vmax
    P = []
    for i in range(20):
        ecit.setVmax(reaction, X[i])
        P.append(ecit.comproducti())
    plt.plot(X,P)
    plt.ylim(0.0007, 0.0017)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.xticks(list(pltrange*initVmax), ticks)

    filename = '20180627TEST_' + reaction + '.png'
    plt.savefig(filename, bbox_inches='tight')
    plt.gcf().clear()

    ecit.setVmax(reaction, initVmax)

    # Write data to file
    data = [str(reaction), str(X), str(P)]
    with open("20180627TEST.txt", 'a') as fobj:
        fobj.write('\n'.join(data))
        fobj.write('\n')

# To do: make this faster (takes about 4 minutes)
