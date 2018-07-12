#! /usr/bin/env python
from __future__ import division, print_function
import libsbml
import roadrunner

"""Read and print main model features"""
reader = libsbml.SBMLReader()
document = reader.readSBMLFromFile("E_coli_Millard2016.xml")
model = document.getModel()

#model.getParameter('FEED').setValue(0.5) # Uncomment this line if you want to change the feed (units are mmol/s)
#reac = model.getReaction("FBA") # Uncomment this and the lines below if you want to change the Vmax of reaction FBA
#lo = reac.getListOfAllElements()
#for el in lo:
#    if type(el) is libsbml.Parameter:
#        if el.getId() == 'Vmax':
#            el.setValue(1.23456)

print("Number of species:", model.getNumSpecies())
print("Number of reactions:", model.getNumReactions())
print("Glucose feeding rate:", model.getParameter('FEED').getValue(), "mmol/s")

"""Simulate"""
rr = roadrunner.RoadRunner(libsbml.writeSBMLToString(document))
time0 = 0  # initial time
timef = 200*3600  # final time
npoints = 1000 # number of points
result = rr.simulate(time0, timef, npoints)

"""Print results"""
print("Maximum derivative at final time:", max(abs(rr.model.getFloatingSpeciesConcentrationRates())))
growthidx = rr.model.getReactionIds().index("GROWTH")
print("Growth rate at final time:", rr.model.getReactionRates()[growthidx]*3600, "h-1") # This should be the steady state growth rate if a steady state is reached, i.e., if the above maximum derivative is low
print("Concentrations")
for nospec, spec in enumerate(rr.model.getFloatingSpeciesIds()):
    print(spec, ":", rr.model.getFloatingSpeciesConcentrations()[nospec])
print("Reaction rates")
for noreac, reac in enumerate(rr.model.getReactionIds()):
    print(reac, ":", rr.model.getReactionRates()[noreac])
