#!/usr/bin/env python
# VARIANT: using scipy - differential_evolution
# Compute optimal absolute Vmax values for n-tuples of enzymes for maximising
# citramalate productivity, using differential evolution

# ROADRUNNER HAS MEMORY LEAK

from __future__ import division, print_function
import ecolicitra
import itertools
import numpy as np
import matplotlib.pyplot as plt
import roadrunner
import libsbml
import time
from scipy.optimize import differential_evolution

# REDEFINE DIFFERENTIAL EVOULTION PARAMETERS HERE

# F ...
#   IF DITHER DEFINE AS A TUPLE (x,y)
#   IF NO DITHER JUST DEFINE AS A FLOAT
mutation = (0.5,1)

crossp = 0.9784 # CR
popsize = 28 # POPULATION SIZE
its = 10 # MAX NUMBER OF GENERATIONS
strategy='best1bin' # DIFFERENTIAL EVOULTION STRATEGY
# REDEFINE NUMBER OF TUPLES (couples, triples...) HERE
n = 41

# REDEFINE Vmax RANGE HERE
boundsrel = [(0.5, 10.0)] * n

# REDEFINE LIST OF REACTIONS HERE
listofreactions = ["ACEA", "ACEB", "ACK", "ACN_1", "ACN_2", "ACS", "ATP_syn", "CITRA_SYN", "CYTBO", "EDA", "EDD", "ENO", "FBA", "FBP", "FUMA", "GDH", "GLT", "GND", "GPM", "LPD", "MAD", "MDH", "MQO", "PCK", "PDH", "PFK", "PGI", "PGK", "PGL", "PIT", "PPC", "PPS", "PTA", "PYK", "RPE", "RPI", "SDH", "SK", "SQR", "TPI", "ZWF"]

# Create kinetic model
include_CITRA = True
ecit = ecolicitra.ecolicit(include_CITRA = include_CITRA)
ecit.setVmax('CITRA_SYN', 4.0)
ecit.time0 = 0
ecit.timef = 2*3600
ecit.npoints = 100

def choose(mylist, n):
    return list(itertools.combinations(mylist, n))

def productivity(r, x):
    # r is a list of n reactions
    # x is a numpy axis with n elements
    for i in range(n):
        ecit.setVmax(r[i], x[i])
    return ecit.comproducti()

boundsrel = np.asarray(boundsrel)
combolist = choose(listofreactions, n)

# overwrites existing file
with open('de.txt', 'w') as f:
    pass

# main
for combo in combolist:
    start_time = time.time() # time tracking

    print(combo)
    VmaxI = [ecit.getVmax(i) for i in combo]
    print(VmaxI)
    bounds = (VmaxI*(boundsrel.T)).T
    print(bounds)

    def fobj(x):
        return -productivity(combo, x)

    # computation
    deresult = differential_evolution(fobj, bounds, strategy=strategy,
                maxiter=its, popsize=popsize, mutation=mutation,
                recombination=crossp, disp=True)
    result = (deresult.x, deresult.fun)

    elapsed_time = time.time() - start_time
    print(elapsed_time) # time tracking

    # reassigns Vmaxes
    for i in range(len(combo)):
        ecit.setVmax(combo[i], VmaxI[i])

    # printing/writing results
    print(result)
    with open('de.txt', 'a') as f:
        f.write(str(combo) + '\n')
        f.write(str(result[0]) + '\n')
        f.write('Fitness: ' + str(result[1]) + '\n')
