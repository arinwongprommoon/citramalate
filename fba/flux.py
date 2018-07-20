#!/usr/bin/env python
from __future__ import division, print_function
import libsbml
import roadrunner
import numpy as np

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

# Stores list of ALL reactions
listofreactions = model.getListOfReactions()

# Fuck around with stuff
# (note to self: remove the word 'fuck' from everything later)
reactiontofuckwith = 'GLT'
start = 0.1
end = 10.0
points = 20
V = wtVmaxes[reactiontofuckwith]
X = np.linspace(start*V, end*V, points, endpoint=True)

fluxdata = np.empty(shape=(68,points))

i = 0
# Loops that actually fucks around with stuff
for xx in X:
    setVmax(reactiontofuckwith, xx)

    # Simulate
    rr = roadrunner.RoadRunner(libsbml.writeSBMLToString(document))
    result = rr.simulate(0, 7200, 100)

    # Reaction rates
    # print("loop " + str(i+1))
    for noreac, reac in enumerate(rr.model.getReactionIds()):
        fluxdata[noreac][i] = rr.model.getReactionRates()[noreac]

    i += 1

for noreac, reac in enumerate(rr.model.getReactionIds()):
    print(reac, ": min ", min(fluxdata[noreac]), " max ", max(fluxdata[noreac]))
