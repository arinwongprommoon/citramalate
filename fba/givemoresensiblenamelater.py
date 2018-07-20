#!/usr/bin/env python
from __future__ import division, print_function
import libsbml
import roadrunner

# Reads wild-type model
reader = libsbml.SBMLReader()
document = reader.readSBMLFromFile("../kinetic/E_coli_Millard2016.xml")
model = document.getModel()

# Stealing useful functions from ecolicita
# class modify:
#     def __init__(self, sbmlfile = "../kinetic/E_coli_Millard2016.xml"):
#         reader = libsbml.SBMLReader()
#         self.document = reader.readSBMLFromFile(sbmlfile)
#         self.model = self.document.getModel()
#
#     def setVmax(self, reacId, value):
#         # Units: mM/s
#         reac = self.model.getReaction(reacId)
#         lo = reac.getListOfAllElements()
#         for el in lo:
#             if type(el) is libsbml.Parameter:
#                 if el.getId() == 'Vmax':
#                     el.setValue(value)
#
#     def getVmax(self, reacId):
#         # get Vmax of reaction reacId. See setVmax for info about units
#         reac = self.model.getReaction(reacId)
#         lo = reac.getListOfAllElements()
#         for el in lo:
#             if type(el) is libsbml.Parameter:
#                 if el.getId() == 'Vmax':
#                     return el.getValue()
#
#     def getVmaxes(self):
#         # get Vmaxes of reactions. See setVmax for info about units
#         Vmaxes = {}
#         for reac in self.model.getListOfReactions():
#             vm = self.getVmax(reac.id)
#             if vm:
#                 Vmaxes[reac.id] = vm
#         self.reacVmaxes = sorted(Vmaxes) # ids of reactions sorted alphabetically that have Vmax
#         self.iniVmaxes = [Vmaxes[r] for r in self.reacVmaxes] # initial values of Vmax (as in the kinetic model)

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

def getVmaxes():
    # get Vmaxes of reactions. See setVmax for info about units
    Vmaxes = {}
    for reac in model.getListOfReactions():
        vm = getVmax(reac.id)
        if vm:
            Vmaxes[reac.id] = vm
    reacVmaxes = sorted(Vmaxes) # ids of reactions sorted alphabetically that have Vmax
    iniVmaxes = [Vmaxes[r] for r in reacVmaxes] # initial values of Vmax (as in the kinetic model)

listofreactions = model.getListOfReactions()

# Simulate
rr = roadrunner.RoadRunner(libsbml.writeSBMLToString(document))
time0 = 0  # initial time
timef = 200*3600  # final time
npoints = 1000 # number of points
result = rr.simulate(time0, timef, npoints)

# Reaction rates
print("Reaction rates")
for noreac, reac in enumerate(rr.model.getReactionIds()):
    print(reac, ":", rr.model.getReactionRates()[noreac])
