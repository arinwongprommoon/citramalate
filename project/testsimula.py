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
timef = 100  # final time (s)
npoints = 1000 # number of points
result = rr.simulate(time0, timef, npoints)

rr.plot()
