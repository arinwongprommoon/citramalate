#!/usr/bin/env python

from __future__ import division, print_function
from ecolicitra_copy import ecolicit, mmCITRA, mmGLC
import numpy as np
import roadrunner
import libsbml

include_CITRA = True
ecit = ecolicit(include_CITRA = include_CITRA)
ecit.setVmax('CITRA_SYN', 4.0)

ecit.time0 = 0
ecit.timef = 2*3600 # final simulation time in seconds
ecit.npoints = 100 # Number of points to be computed in the simulation

class Coupl(object):
    """Calculate couples and a bunch of other stuff (will write more descriptive thing later)"""
    def __init__(self, xstart=0.5, xend=2.0, ystart=0.5, yend=2.0, points=10):
        ecit.getVmaxes()
        self.listofreactions = ecit.reacVmaxes
        self.xstart = xstart
        self.xend = xend
        self.ystart = ystart
        self.yend = yend
        self.points = points

        self.listofreactions = ['CITRA_SYN', 'GLT', 'LPD']
        #self.listofreactions = ['CITRA_SYN', 'GLT', 'LPD']

    def writeFileHeader(self, filename = "COUPLESDATA.txt"):
        with open(filename, 'w') as fobj:
            fobj.write('time0 ' + str(ecit.time0) + '\n')
            fobj.write('timef ' + str(ecit.timef) + '\n')
            fobj.write('npoints ' + str(ecit.npoints) + '\n')
            fobj.write('XVmax ' + str(self.xstart) + ' to ' + str(self.xend) + '\n')
            fobj.write('YVmax ' + str(self.ystart) + ' to ' + str(self.yend) + '\n')
            fobj.write('datapoints ' + str(self.points) + '\n\n')

    def writeToText(self, names, data, filename = "COUPLESDATA.txt"):
        with open(filename, 'a') as fobj:
            fobj.write(str(names[0]) + ' vs ' + str(names[1]) + '\n')
            fobj.write('X ' + str(data[0]) + '\n')
            fobj.write('Y ' + str(data[1]) + '\n')
            fobj.write(str(data[2]) + '\n\n')

    def writeToNpz(self, data, filename = "COUPLESDATA.npz"):
        np.savez(filename, *data)

    def computeCouples(self, writemethod = 'writeToBoth'):
        """
        writemethod
            writeToText: TXT only
            writeToNpz: NPZ only
            writeToBoth: both
        """
        for XRxn in self.listofreactions:
            XVmaxI = ecit.getVmax(XRxn)
            X = np.linspace(self.xstart*XVmaxI, self.xend*XVmaxI, self.points, endpoint=True)
            for YRxn in [r for r in self.listofreactions if self.listofreactions.index(r)>self.listofreactions.index(XRxn)]:
                print(XRxn + ' vs ' + YRxn)
                YVmaxI = ecit.getVmax(YRxn)
                Y = np.linspace(self.ystart*YVmaxI, self.yend*YVmaxI, self.points, endpoint=True)
                P = np.zeros((self.points, self.points))
                for ii in range(self.points):
                    ecit.setVmax(XRxn, X[ii])
                    for jj in range(self.points):
                        ecit.setVmax(YRxn, Y[jj])
                        P[ii, jj] = ecit.comproducti()
                ecit.setVmax(YRxn, YVmaxI)
                if writemethod == 'writeToText':
                    names = [XRxn, YRxn]
                    data = [X, Y, P]
                    self.writeToText(names, data)
                elif writemethod == 'writeToNpz':
                    data = [X, Y, P]
                    filename = "COUPLESDATA_" + str(XRxn) + "_" + str(YRxn) + ".npz"
                    self.writeToNpz(data, filename)
                data = [X, Y, P]
                names = [XRxn, YRxn]
                if writemethod == 'writeToText':
                    self.writeToText(names, data)
                elif writemethod == 'writeToNpz':
                    filename = "COUPLESDATA_" + str(XRxn) + "_" + str(YRxn) + ".npz"
                    self.writeToNpz(data, filename)
                elif writemethod == 'writeToBoth':
                    self.writeToText(names, data)
                    filename = "COUPLESDATA_" + str(XRxn) + "_" + str(YRxn) + ".npz"
                    self.writeToNpz(data, filename)
            ecit.setVmax(XRxn, XVmaxI)
