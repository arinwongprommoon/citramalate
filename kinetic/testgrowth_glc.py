#! /usr/bin/env python
from __future__ import division, print_function
import libsbml
import roadrunner
import numpy as np
import csv

"""Read and print main model features"""
reader = libsbml.SBMLReader()
document = reader.readSBMLFromFile("E_coli_Millard2016_CITRA.xml")
model = document.getModel()

def setVmax(reacId, value):
    # Units: mM/s
    reac = model.getReaction(reacId)
    lo = reac.getListOfAllElements()
    for el in lo:
        if type(el) is libsbml.Parameter:
            if el.getId() == 'Vmax':
                el.setValue(value)

def setKm(reacId, value):
    reac = model.getReaction(reacId)
    lo = reac.getListOfAllElements()
    for el in lo:
        if type(el) is libsbml.Parameter:
            if el.getId() == 'Km':
                el.setValue(value)

#setVmax('XCH_GLC', 2000)

X = np.linspace(0, 1, 21, endpoint=True)

mylist = []

for feed in X:

    model.getParameter('FEED').setValue(feed) # Uncomment this line if you want to change the feed (units are mmol/s)
    print("Glucose feeding rate:", model.getParameter('FEED').getValue(), "mmol/s")

    """Simulate"""
    rr = roadrunner.RoadRunner(libsbml.writeSBMLToString(document))
    time0 = 0  # initial time
    timef = 200*3600  # final time
    npoints = 1000 # number of points
    result = rr.simulate(time0, timef, npoints)

    """Print results"""
    print("Maximum derivative at final time:", max(abs(rr.model.getFloatingSpeciesConcentrationRates())))

    specselection = ["GLCx", "GLCp", "G6P", 'iGROWTH']
    specids = [rr.model.getFloatingSpeciesIds().index(s) for s in specselection]   
    print("Concentrations")
    for i in xrange(len(specselection)):
        print(specselection[i], ":", rr.model.getFloatingSpeciesConcentrations()[specids[i]])
        mylist.append(rr.model.getFloatingSpeciesConcentrations()[specids[i]])

    reacselection = ['GROWTH', 'GLC_feed', 'XCH_GLC', 'PTS_0', 'PTS_1', 'PTS_2', 'PTS_3', 'PTS_4']
    reacids = [rr.model.getReactionIds().index(r) for r in reacselection]   
    print("Reaction rates")  
    for i in xrange(len(reacselection)):
        print(reacselection[i], ":", rr.model.getReactionRates()[reacids[i]])
        
    print('\n')
    
myarray = np.reshape(mylist, (21,4))
np.savetxt('testgrowth_glc.csv', myarray, delimiter=',')
