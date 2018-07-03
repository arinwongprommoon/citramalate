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
start = 0.1
end = 1.0
points = 50

# Clears data in file, writes header
with open("VMAXDATA.txt", 'w') as fobj:
    fobj.write('time0 ' + str(ecit.time0) + '\n')
    fobj.write('timef ' + str(ecit.timef) + '\n')
    fobj.write('Vmax ' + str(start) + 'to' + str(end) + '\n\n')

ecit.getVmaxes()
listofreactions = ['MQO', 'ATP_MAINTENANCE', 'CITRA_SYN', 'CYTBO', 'EDA', 'GDH', 'PFK']

for reaction in listofreactions:
    print(reaction)
    initVmax = ecit.getVmax(reaction)

    # Makes the data
    X = np.linspace(start*initVmax, end*initVmax, points, endpoint=True) # sets X values from 0.5*default Vmax to 2.0*default Vmax
    P = []
    for i in range(points):
        ecit.setVmax(reaction, X[i])
        P.append(ecit.comproducti())

    ecit.setVmax(reaction, initVmax)

    # Write data to file
    data = [str(reaction), str(X.tolist()), str(P)]
    with open("VMAXDATA.txt", 'a') as fobj:
        fobj.write('\n'.join(data))
        fobj.write('\n')
