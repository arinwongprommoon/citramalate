#!/usr/bin/env python
# Aim for it to work for couples only for now, but sections of code purposefully
# written for easy change to triples, etc (e.g. changing 2 to 3)
# ROADRUNNER HAS MEMORY LEAK

from __future__ import division, print_function
import ecolicitra # note ecolicitra instead of ecolicitra_copy
import itertools
import gc
import numpy as np
import roadrunner
import libsbml

# REDEFINE NUMBER OF TUPLES (couples, triples...) HERE
n = 2

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

# DE algorithm adapted from Pablo R Mier
def de(fobj, bounds, mut=0.8, crossp=0.7, popsize=10, its=10):
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
    yield best, fitness[best_idx]

boundsrel = [(0.5, 2.0)] * n
boundsrel = np.asarray(boundsrel)

listofreactions = ['CITRA_SYN', 'GLT', 'LPD']
combolist = choose(listofreactions, n)

# overwrites existing file
with open('de.txt', 'w') as f:
    pass

# main
for combo in combolist:
    gc.collect() # NOT USEFUL IN DEALING WITH MEMORY LEAK
    print(combo)
    VmaxI = [ecit.getVmax(i) for i in combo]
    print(VmaxI)
    bounds = (VmaxI*(boundsrel.T)).T
    print(bounds)
    # probably not the best way to do this, but I'll get it to work first
    def fobj(x):
        # minus sign there because we're maximising productivity, but DE
        # minimises the value of a function
        return -productivity(combo, x)

    # computation
    result = list(de(fobj, bounds))

    # printing results
    print(result[-1])
    with open('de.txt', 'a') as f:
        f.write(str(combo) + '\n')
        f.write(str(result[-1][0]) + '\n')
        f.write('Fitness: ' + str(result[-1][1]) + '\n')

    # reassigns Vmaxes
    for i in range(len(combo)):
        ecit.setVmax(combo[i], VmaxI[i])

    # hit return to continue -> NOT USEFUL IN DEALING WITH MEMORY LEAK
    print('Next combination?')
    raw_input()
