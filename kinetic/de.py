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

# REDEFINE NUMBER OF TUPLES (couples, triples...) HERE
n = 5

# REDEFINE Vmax RANGE HERE
boundsrel = [(0.1, 10.0)] * n

# REDEFINE LIST OF REACTIONS HERE
#listofreactions = ['PGI', 'PFK', 'FBA', 'GDH', 'PGK', 'GPM', 'ENO', 'PYK', 'PDH']
listofreactions = ['PGI', 'GDH', 'ENO', 'PYK', 'PDH']

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

generations = []

# DE algorithm adapted from Pablo R Mier
def de(fobj, bounds, mut=0.6301, crossp=0.7122, popsize=17, its=30):
    dimensions = len(bounds)
    # Initialisation
    pop = np.random.rand(popsize, dimensions)
    min_b, max_b = np.asarray(bounds).T
    diff = np.fabs(min_b - max_b)
    pop_denorm = min_b + pop*diff
    fitness = np.asarray([fobj(ind) for ind in pop_denorm])
    best_idx = np.argmin(fitness)
    for i in range(its):
        for j in range(popsize):
            # Mutation
            idxs = [idx for idx in range(popsize) if idx != j]
            a, b, c = pop[np.random.choice(idxs, 3, replace = False)]
            mutant = np.clip(a + mut*(b-c), 0, 1)
            # Recombination
            cross_points = np.random.rand(dimensions) < crossp
            if not np.any(cross_points):
                cross_points[np.random.randint(0, dimensions)] = True
            trial = np.where(cross_points, mutant, pop[j])
            trial_denorm = min_b + trial*diff
            # Selection
            f = fobj(trial_denorm)
            if f < fitness[j]:
                fitness[j] = f
                pop[j] = trial
                if f < fitness[best_idx]:
                    best_idx = j
                    best = trial_denorm
                    generations.append(i)
                    # I don't actually need this print line but
                    # I like seeing that stuff happens while I run code
                    print(str(i) + ' ' + str(j))
                    yield best, fitness[best_idx]

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
    result = list(de(fobj, bounds))

    elapsed_time = time.time() - start_time
    print(elapsed_time) # time tracking

    # reassigns Vmaxes
    for i in range(len(combo)):
        ecit.setVmax(combo[i], VmaxI[i])
    
    # printing/writing results
    print(result[-1])
    with open('de.txt', 'a') as f:
        f.write(str(combo) + '\n')
        f.write(str(result[-1][0]) + '\n')
        f.write('Fitness: ' + str(result[-1][1]) + '\n')

    # plot convergence
    x, f = zip(*result)
    ff = list(f)
    plt.plot(generations, ff)
    plt.show()
