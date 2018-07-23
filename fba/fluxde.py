#!/usr/bin/env python

from __future__ import division, print_function
import itertools
import numpy as np
import matplotlib.pyplot as plt
import roadrunner
import libsbml
import time

# Reads wild-type model
reader = libsbml.SBMLReader()
document = reader.readSBMLFromFile("../kinetic/E_coli_Millard2016.xml")
model = document.getModel()

# Stealing useful functions from ecolicita
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

# Get Vmax values of all reactions that have Vmaxes and stores them in the
# dictionary wtVmaxes
Vmaxes = {}
for reac in model.getListOfReactions():
    vm = getVmax(reac.id)
    if vm:
        Vmaxes[reac.id] = vm
reacVmaxes = sorted(Vmaxes) # ids of reactions sorted alphabetically that have Vmax
iniVmaxes = [Vmaxes[r] for r in reacVmaxes] # initial values of Vmax (as in the kinetic model)
wtVmaxes = dict(zip(reacVmaxes, iniVmaxes))

# REDEFINE NUMBER OF TUPLES (couples, triples...) HERE
n = 2

# REDEFINE Vmax RANGE HERE
boundsrel = [(0.1, 10.0)] * n

# REDEFINE LIST OF REACTIONS HERE
#listofreactions = reacVmaxes
listofreactions = ['GLT', 'GDH', 'LPD']

def choose(mylist, n):
    return list(itertools.combinations(mylist, n))

def flux(r, x):
    # r is a list of n reactions
    # x is a numpy axis with n elements
    for i in range(n):
        setVmax(r[i], x[i])

    # Simulate
    rr = roadrunner.RoadRunner(libsbml.writeSBMLToString(document))
    result = rr.simulate(0, 7200, 100)

    noreac = 0
    return rr.model.getReactionRates()[noreac]

generations = []

# DE algorithm adapted from Pablo R Mier
def de(fobj, bounds, mut=0.6607, crossp=0.9426, popsize=28, its=10):
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
with open('fluxde.txt', 'w') as f:
    pass

# main
for combo in combolist:
    start_time = time.time() # time tracking

    print(combo)
    VmaxI = [getVmax(i) for i in combo]
    print(VmaxI)
    bounds = (VmaxI*(boundsrel.T)).T
    print(bounds)

    def fobj(x):
        #return flux(combo, x) # minimum
        return -flux(combo, x) # maximum

    # computation
    result = list(de(fobj, bounds))

    elapsed_time = time.time() - start_time
    print(elapsed_time) # time tracking

    # reassigns Vmaxes
    for i in range(len(combo)):
        setVmax(combo[i], VmaxI[i])

    # printing/writing results
    print(result[-1])
    with open('fluxde.txt', 'a') as f:
        f.write(str(combo) + '\n')
        f.write(str(result[-1][0]) + '\n')
        f.write('Fitness: ' + str(result[-1][1]) + '\n')
