#!/usr/bin/env python

from __future__ import division, print_function
from ecolicitra_copy import ecolicit, mmCITRA, mmGLC
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

xstart = 0.5
xend = 2.0
ystart = 0.5
yend = 2.0
points = 10

# Clears data in file, writes header
with open("COUPLESDATA.txt", 'w') as fobj:
    fobj.write('time0 ' + str(ecit.time0) + '\n')
    fobj.write('timef ' + str(ecit.timef) + '\n')
    fobj.write('npoints ' + str(ecit.npoints) + '\n')
    fobj.write('XVmax ' + str(xstart) + ' to ' + str(xend) + '\n')
    fobj.write('YVmax ' + str(ystart) + ' to ' + str(yend) + '\n')
    fobj.write('datapoints ' + str(points) + '\n\n')

ecit.getVmaxes()
# listofreactions = ecit.reacVmaxes
listofreactions = ['CITRA_SYN', 'GLT', 'LPD']

for XRxn in listofreactions:
    XVmaxI = ecit.getVmax(XRxn)
    X = np.linspace(xstart*XVmaxI, xend*XVmaxI, points, endpoint=True)
    for YRxn in [r for r in listofreactions if listofreactions.index(r)>listofreactions.index(XRxn)]:
        print(XRxn + ' vs ' + YRxn)
        YVmaxI = ecit.getVmax(YRxn)
        Y = np.linspace(ystart*YVmaxI, yend*YVmaxI, points, endpoint=True)
        P = np.zeros((points, points))
        for ii in range(points):
            ecit.setVmax(XRxn, X[ii])
            for jj in range(points):
                ecit.setVmax(YRxn, Y[jj])
                P[ii, jj] = ecit.comproducti()
        ecit.setVmax(YRxn, YVmaxI)
        #print('X ' + str(X))
        #print('Y ' + str(Y))
        #print(P)
        with open("COUPLESDATA.txt", 'a') as fobj:
            fobj.write(str(XRxn) + ' vs ' + str(YRxn) + '\n')
            fobj.write('X ' + str(X) + '\n')
            fobj.write('Y ' + str(Y) + '\n')
            fobj.write(str(P) + '\n\n')
    ecit.setVmax(XRxn, XVmaxI)
