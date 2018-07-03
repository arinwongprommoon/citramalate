#!/usr/bin/env python

from __future__ import division, print_function
from ecolicitra_copy import ecolicit, mmCITRA, mmGLC
import numpy as np
import roadrunner
import libsbml

import time

class Coupl(ecolicit):
    """
        Class that contains functions to vary Vmax values for pairs of enzymes
        within a specified set of enzymes (default is the entire list of enzymes
        in the Millard model plus citramalate synthesis) and simulate
        citramalate productivity accordingly. Outputs data into a text (TXT)
        file or as numpy arrays (NPZ).

        Arguments:
            xstart = lower limit of Vmax values for 1st enzyme
            xend = upper limit of Vmax values for 1st enzyme
            ystart = lower limit of Vmax values for 2nd enzyme
            yend = upper limit of Vmax values ofr 2nd enzyme
            points = number of data points to compute (same for each axis)
    """
    def __init__(self, xstart=0.5, xend=2.0, ystart=0.5, yend=2.0, points=10):
        ecolicit.__init__(self, sbmlfile = "E_coli_Millard2016.xml", Vmax = 4.0, Km = 0.495, include_CITRA = True, initial_CITRA = 0.0)

        self.time0 = 0
        self.timef = 2*3600
        self.npoints = 100

        self.getVmaxes()
        self.listofreactions = self.reacVmaxes
        self.xstart = xstart
        self.xend = xend
        self.ystart = ystart
        self.yend = yend
        self.points = points

    def writeFileHeader(self, filename = "COUPLESDATA.txt"):
        """
            Writes file header for output TXT file, containing info relevant to
            reaction.
        """
        with open(filename, 'w') as fobj:
            fobj.write('time0 ' + str(self.time0) + '\n')
            fobj.write('timef ' + str(self.timef) + '\n')
            fobj.write('npoints ' + str(self.npoints) + '\n')
            fobj.write('XVmax ' + str(self.xstart) + ' to ' + str(self.xend) + '\n')
            fobj.write('YVmax ' + str(self.ystart) + ' to ' + str(self.yend) + '\n')
            fobj.write('datapoints ' + str(self.points) + '\n\n')

    def writeToText(self, names, data, filename = "COUPLESDATA.txt"):
        """
            Writes the actual data for output TXT file
        """
        with open(filename, 'a') as fobj:
            fobj.write(str(names[0]) + ' vs ' + str(names[1]) + '\n')
            fobj.write('X ' + str(data[0]) + '\n')
            fobj.write('Y ' + str(data[1]) + '\n')
            fobj.write(str(data[2]) + '\n\n')

    def writeToNpz(self, data, filename = "COUPLESDATA.npz"):
        """Writes data as numpy array NPZ files"""
        np.savez(filename, *data)

    def computeCouples(self, writemethod = 'writeToBoth'):
        """
            Computes citramalate productivity for pairs of enzymes after
            varying Vmax values.

            Argument:
            writemethod =
                writeToText: TXT only
                writeToNpz: NPZ only
                writeToBoth: both (default)
        """
        # TIME TO EXECUTE INCREASES LINEARLY WITH EACH LOOP - fix this
        for XRxn in self.listofreactions:
            # Computation
            XVmaxI = self.getVmax(XRxn)
            X = np.linspace(self.xstart*XVmaxI, self.xend*XVmaxI, self.points, endpoint=True)
            for YRxn in [r for r in self.listofreactions if self.listofreactions.index(r)>self.listofreactions.index(XRxn)]:
                t0 = time.clock() # track time elapsed
                print(XRxn + ' vs ' + YRxn) # track reaction
                YVmaxI = self.getVmax(YRxn)
                Y = np.linspace(self.ystart*YVmaxI, self.yend*YVmaxI, self.points, endpoint=True)
                P = np.zeros((self.points, self.points)) # initialise
                for ii in range(self.points):
                    self.setVmax(XRxn, X[ii])
                    for jj in range(self.points):
                        self.setVmax(YRxn, Y[jj])
                        P[ii, jj] = self.comproducti()
                self.setVmax(YRxn, YVmaxI)

                # Writing data
                data = [X, Y, P]
                names = [XRxn, YRxn]
                if writemethod == 'writeToText':
                    self.writeToText(names, data)
                elif writemethod == 'writeToNpz':
                    filename = "COUPLESDATA-" + str(XRxn) + "-" + str(YRxn) + ".npz"
                    self.writeToNpz(data, filename)
                elif writemethod == 'writeToBoth':
                    self.writeToText(names, data)
                    filename = "COUPLESDATA-" + str(XRxn) + "-" + str(YRxn) + ".npz"
                    self.writeToNpz(data, filename)

                print(time.clock()) # track time elapsed

            self.setVmax(XRxn, XVmaxI)
