#!/usr/bin/env python

from __future__ import division, print_function
from ecolicitra_copy import ecolicit, mmCITRA, mmGLC
import numpy as np
import itertools
import roadrunner
import libsbml

import time
import gc

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
            yend = upper limit of Vmax values for 2nd enzyme
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

    def pairs(self, mylist):
        return list(itertools.combinations(mylist, 2)) # generator

    def computeCouples(self, writemethod='writeToBoth'):
        """
            Computes citramalate productivity for pairs of enzymes after
            varying Vmax values.

            Argument:
            writemethod =
                writeToText: TXT only
                writeToNpz: NPZ only
                writeToBoth: both (default)
        """
        # Generate list of pairs
        pairslist = self.pairs(self.listofreactions)

        # Initialise to first element of pairslist (useful when calling getVmax only when needed)
        XRxn = pairslist[0][0]
        XVmaxI = self.getVmax(XRxn)
        X = np.linspace(self.xstart*XVmaxI, self.xend*XVmaxI, self.points, endpoint=True)

        YRxn = pairslist[0][1]
        YVmaxI = self.getVmax(YRxn)
        Y = np.linspace(self.ystart*YVmaxI, self.yend*YVmaxI, self.points, endpoint=True)

        ## BIG FOR LOOP ORIGINALLY HERE
        # read the index of the pair from file
        with open('couplesi.txt', 'r') as fobj:
            mystr = fobj.readline()
            idx = int(mystr.strip('\n'))

        # stops if all elements gone through
        if idx >= len(pairslist):
            return 0

        pair = pairslist[idx]

        # Call getVmax() and generate arrays only when needed
        if pair[0] != XRxn:
            XRxn = pair[0]
            XVmaxI = self.getVmax(XRxn)
            X = np.linspace(self.xstart*XVmaxI, self.xend*XVmaxI, self.points, endpoint=True)
        if pair[1] != YRxn:
            YRxn = pair[1]
            YVmaxI = self.getVmax(YRxn)
            Y = np.linspace(self.ystart*YVmaxI, self.yend*YVmaxI, self.points, endpoint=True)

        print(XRxn + ' vs ' + YRxn) # track reaction

        start_time = time.time() # time tracking

        # The real bit
        ij = 0
        #Ptemp = []
        Ptemp = np.empty(self.points**2) # initialise
        for (xi, yi) in itertools.product(X, Y):
            self.setVmax(XRxn, xi)
            self.setVmax(YRxn, yi)
            Ptemp[ij] = self.comproducti()
            ij += 1
            gc.collect()

        elapsed_time = time.time() - start_time
        print(elapsed_time) # time tracking

        self.setVmax(XRxn, XVmaxI)
        self.setVmax(YRxn, YVmaxI)

        P = Ptemp.reshape((self.points, self.points))

        # write the index of the pair to file
        idx += 1
        with open('couplesi.txt', 'w') as fobj:
            fobj.write(str(idx))

        # Writing data
        data = [X, Y, P]
        names = [XRxn, YRxn]
        if writemethod == 'writeToText':
            filename = "COUPLESDATA-" + str(XRxn) + "-" + str(YRxn) + ".txt"
            self.writeToText(names, data, filename)
        elif writemethod == 'writeToNpz':
            filename = "COUPLESDATA-" + str(XRxn) + "-" + str(YRxn) + ".npz"
            self.writeToNpz(data, filename)
        elif writemethod == 'writeToBoth':
            filenamet = "COUPLESDATA-" + str(XRxn) + "-" + str(YRxn) + ".txt"
            self.writeToText(names, data, filenamet)
            filenamen = "COUPLESDATA-" + str(XRxn) + "-" + str(YRxn) + ".npz"
            self.writeToNpz(data, filenamen)
        gc.collect()
