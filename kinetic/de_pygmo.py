#!/usr/bin/env python
# VARIANT: using pygmo - sade (self-adaptive differential evolution)
# Compute optimal absolute Vmax values for n-tuples of enzymes for maximising
# citramalate productivity, using differential evolution

# ROADRUNNER HAS MEMORY LEAK

from __future__ import division, print_function
import ecolicitra
import itertools
import numpy as np
import matplotlib.pyplot as plt
import time
import pygmo as pg

# REDEFINE DIFFERENTIAL EVOULTION PARAMETERS HERE
gen = 10 # GENERATIONS
np = 28 # POPULATION SIZE

# REDEFINE NUMBER OF TUPLES (couples, triples...) HERE
n = 3

# REDEFINE Vmax RANGE HERE
boundsrel = [(0.5, 10.0)] * n

# REDEFINE LIST OF REACTIONS HERE
listofreactions = ["ACEA", "ACEB", "ACK"]

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

    # redefine bounds to satisfy pygmo
    bounds_pygmo = ([],[])
    for i in bounds:
        bounds_pygmo[0].append(i[0])
        bounds_pygmo[1].append(i[1])
    print(bounds_pygmo)

    # not the best way, trying to get it to at least do something
    class citraprod:
        def __init__(self, dim):
            self.dim = n
        def fitness(self, x):
            return [-productivity(combo, x)]
        def get_bounds(self):
            return(bounds_pygmo)

    prob = pg.problem(citraprod(n))
    algo = pg.algorithm(pg.sade(gen = gen))
    algo.set_verbosity(1)
    pop = pg.population(prob, np)
    pop = algo.evolve(pop)
    result = pop.champion_f
    print(result)

    print('Vmaxes')
    for i in combo:
        print(ecit.getVmax(i))

    elapsed_time = time.time() - start_time
    print(elapsed_time) # time tracking

    # printing/writing results: enzymes, their final Vmax values, solution
    print(result)
    with open('de.txt', 'a') as f:
        f.write(str(combo) + '\n')
        for i in combo:
            f.write(str(ecit.getVmax(i)) + '\n')
        f.write('Fitness: ' + str(result) + '\n')

    # reassigns Vmaxes
    for i in range(len(combo)):
        ecit.setVmax(combo[i], VmaxI[i])
