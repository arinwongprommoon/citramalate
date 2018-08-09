#!/usr/bin/env python
# Changes Vmax values of one enzyme at a time in the kinetic model and gets
# the maximum and minimum fluxes of each reaction.
# Used to generate 1D boundaries
# Doesn't loop to get around memory leak
from __future__ import division, print_function
import libsbml
import roadrunner
import numpy as np
import sys
import time
import csv

# Reads wild-type model
reader = libsbml.SBMLReader()
document = reader.readSBMLFromFile("../kinetic/E_coli_Millard2016_CITRA.xml")
model = document.getModel()

# start, end, data points
start = 0.1 # START VMAX
end = 10.0 # END VMAX
points = 200 # NUMBER OF DATA POINTS TO COMPUTE
# NUMBER OF REACTIONS IN THE KINETIC MODEL, USUALLY 69
#noofreactions = 69
noofreactions = len(model.getListOfReactions())

# SIMULATION PARAMETERS
timestart = 0
timeend = 7200
rrpoints = 100

# Stealing useful functions from ecolicita
def setVmax(reacId, value):
    # Units: mM/s
    reac = model.getReaction(reacId)
    lo = reac.getListOfAllElements()
    for el in lo:
        if type(el) is libsbml.Parameter:
            if el.getId() == 'Vmax':
                el.setValue(value)

def getVmax(reacId):
    # get Vmax of reaction reacId. See setVmax for info about units
    reac = model.getReaction(reacId)
    lo = reac.getListOfAllElements()
    for el in lo:
        if type(el) is libsbml.Parameter:
            if el.getId() == 'Vmax':
                return el.getValue()

# Get Vmax values of all reactions that have Vmaxes and stores them in the
# dictionary wtVmaxes
Vmaxes = {}
for reac in model.getListOfReactions():
    vm = getVmax(reac.id)
    if vm:
        Vmaxes[reac.id] = vm
reacVmaxes = sorted(Vmaxes) # ids of reactions sorted alphabetically that have Vmax
iniVmaxes = [Vmaxes[r] for r in reacVmaxes] # initial values of Vmax (as in the kinetic model)
wtVmaxes = dict(zip(reacVmaxes, iniVmaxes))

# Modify Vmax of specified reaction
fluxdata = np.empty(shape=(noofreactions,points))
# read the index of the reaction from file
with open('fluxi.txt', 'r') as fobj:
    mystr = fobj.readline()
    idx = int(mystr.strip('\n'))

# stops if all elements gone through
if idx >= len(reacVmaxes):
    print('DONE')
    sys.exit(1)

else:
    reaction = reacVmaxes[idx]

    # main
    start_time = time.time() # time tracking
    print(idx+1, reaction, "varied ---")
    V = wtVmaxes[reaction]
    X = np.linspace(start*V, end*V, points, endpoint=True)
    i = 0
    for xx in X:
        setVmax(reaction, xx)

        # Simulate
        rr = roadrunner.RoadRunner(libsbml.writeSBMLToString(document))
        result = rr.simulate(timestart, timeend, rrpoints)

        # Put reaction rates into array
        for noreac, reac in enumerate(rr.model.getReactionIds()):
            fluxdata[noreac][i] = rr.model.getReactionRates()[noreac]

        i += 1

    # write the index of reaction to file
    idx += 1
    with open('fluxi.txt', 'w') as fobj:
        fobj.write(str(idx))

    # reassign vmax
    setVmax(reaction, V)

    # outputs to both stdout and a csv file
    # I did separate csv files for each reaction modified so that I can catch
    # any irregularities
    filename = reaction + '.csv'
    with open(filename, 'wb') as csvfile:
        fluxwriter = csv.writer(csvfile)
        for noreac, reac in enumerate(rr.model.getReactionIds()):
            print(reac, ": min ", min(fluxdata[noreac]), " max ", max(fluxdata[noreac]))
            fluxwriter.writerow([reac, min(fluxdata[noreac]), max(fluxdata[noreac])])

    elapsed_time = time.time() - start_time
    print("time taken ", elapsed_time) # time tracking
    print("\n")
    sys.exit(2)
