#! /usr/bin/env python
from __future__ import division, print_function
import libsbml
import roadrunner
import numpy as np

"""Read and print main model features"""
reader = libsbml.SBMLReader()
document = reader.readSBMLFromFile("E_coli_Millard2016.xml")
model = document.getModel()

X = np.linspace(0.1, 0.5, 10, endpoint=True)

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
    growthidx = rr.model.getReactionIds().index("GROWTH")


    specselection = ["GLCx", "GLCp"]
    specids = [rr.model.getFloatingSpeciesIds().index(s) for s in specselection]   
    print("Concentrations")
    for i in xrange(len(specselection)):
        print(specselection[i], ":", rr.model.getFloatingSpeciesConcentrations()[specids[i]])

    reacselection = ['GLC_feed', 'XCH_GLC', 'PTS_0', 'PTS_1', 'PTS_2', 'PTS_3', 'PTS_4']
    reacids = [rr.model.getReactionIds().index(r) for r in reacselection]   
    print("Reaction rates")  
    for i in xrange(len(reacselection)):
        print(reacselection[i], ":", rr.model.getReactionRates()[reacids[i]])
        
    print('\n')
