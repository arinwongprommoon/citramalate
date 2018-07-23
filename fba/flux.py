#!/usr/bin/env python
from __future__ import division, print_function
import libsbml
import roadrunner
import numpy as np
import time

# Reads wild-type model
reader = libsbml.SBMLReader()
document = reader.readSBMLFromFile("../kinetic/E_coli_Millard2016.xml")
model = document.getModel()

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

# start, end, data points
start = 0.1
end = 10.0
points = 200
fluxdata = np.empty(shape=(68,points))

j = 0
for reaction in ['GLT']:
    start_time = time.time()
    print(j+1, reaction, "varied ---")
    V = wtVmaxes[reaction]
    X = np.linspace(start*V, end*V, points, endpoint=True)
    i = 0
    for xx in X:
        setVmax(reaction, xx)

        # Simulate
        rr = roadrunner.RoadRunner(libsbml.writeSBMLToString(document))
        result = rr.simulate(0, 7200, 100)

        # Reaction rates
        for noreac, reac in enumerate(rr.model.getReactionIds()):
            fluxdata[noreac][i] = rr.model.getReactionRates()[noreac]
        
        i += 1

    j += 1

    setVmax(reaction, V)

    for noreac, reac in enumerate(rr.model.getReactionIds()):
        print(reac, ": min ", min(fluxdata[noreac]), " max ", max(fluxdata[noreac]))
    elapsed_time = time.time() - start_time
    print("time taken ", elapsed_time)
    print("\n")
