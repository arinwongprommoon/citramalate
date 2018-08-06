#!/usr/bin/env python
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

# REDEFINE NUMBER OF TUPLES (couples, triples...) HERE
n = 2

# REDEFINE Vmax RANGE HERE
boundsrel = [(0.1, 10.0)] * n

# REDEFINE LIST OF REACTIONS HERE
listofreactions = ['GLT', 'LPD']

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
                    
def de(fobj, bounds, mut=0.6607, crossp=0.9426, popsize=28, its=10):
    deresult = differential_evolution(fobj, bounds, strategy='rand1bin', maxiter=its, popsize=popsize, mutation=(0.5,1), recombination=crossp, disp=True)
    return deresult.x, deresult.fun

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
    #result = list(de(fobj, bounds))
    result = de(fobj, bounds)

    elapsed_time = time.time() - start_time
    print(elapsed_time) # time tracking

    # reassigns Vmaxes
    for i in range(len(combo)):
        ecit.setVmax(combo[i], VmaxI[i])
    
    # printing/writing results
    #print(result[-1])
    print(result)
    with open('de.txt', 'a') as f:
        f.write(str(combo) + '\n')
        f.write(str(result[0]) + '\n')
        f.write('Fitness: ' + str(result[1]) + '\n')
