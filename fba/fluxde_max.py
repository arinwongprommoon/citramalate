#!/usr/bin/env python
# Use differential evolution with specified dimensions to vary the Vmax values
# of enzymes in a list in the kinetic model, in order to obtain the minimum
# fluxes through all enzymes in the kinetic model.

from __future__ import division, print_function
import itertools
import numpy as np
import roadrunner
import libsbml
import time
import sys

# COMPREHENSIVE LIST OF ALL REACTIONS
#allreactions = ['PGI', 'PFK', 'FBA', 'TPI', 'GDH', 'PGK', 'GPM', 'ENO', 'PYK', 'ZWF', 'PGL', 'GND', 'RPE', 'RPI', 'X5P_GAP_TKT', 'F6P_E4P_TKT', 'S7P_R5P_TKT', 'F6P_GAP_TAL', 'S7P_E4P_TAL', 'FBP', 'PPC', 'PCK', 'PPS', 'MAD', 'PDH', 'GLT', 'ACN_1', 'ACN_2', 'ICD', 'LPD', 'SK', 'SDH', 'FUMA', 'MQO', 'MDH', 'ACEA', 'ACEB', 'ACEK_1', 'ACEK_2', 'EDD', 'EDA', 'NADH_req', 'PNT_req', 'ADK', 'ATP_syn', 'CYA', 'DOS', 'ACK', 'ACS', 'PTA', 'PTS_0', 'PTS_1', 'PTS_2', 'PTS_3', 'PTS_4', 'GLC_feed', 'CYTBO', 'SQR', 'NDHII', 'GROWTH', 'ATP_MAINTENANCE', 'XCH_GLC', 'PIT', 'XCH_P', 'XCH_ACE1', '_ACE_OUT', 'XCH_ACE2', 'GL6P_HYDROLYSIS']

# REMOVED GLC_FEED 2018-07-25 10:00 BECAUSE IT CAUSES PROBLEMS
# LIST OF REACTIONS TO FIND FLUX BOUNDARIES FOR
allreactions = ['PGI', 'PFK', 'FBA', 'TPI', 'GDH', 'PGK', 'GPM', 'ENO', 'PYK',
'ZWF', 'PGL', 'GND', 'RPE', 'RPI', 'X5P_GAP_TKT', 'F6P_E4P_TKT', 'S7P_R5P_TKT',
'F6P_GAP_TAL', 'S7P_E4P_TAL', 'FBP', 'PPC', 'PCK', 'PPS', 'MAD', 'PDH', 'GLT',
 'ACN_1', 'ACN_2', 'ICD', 'LPD', 'SK', 'SDH', 'FUMA', 'MQO', 'MDH', 'ACEA',
  'ACEB', 'ACEK_1', 'ACEK_2', 'EDD', 'EDA', 'NADH_req', 'PNT_req', 'ADK',
  'ATP_syn', 'CYA', 'DOS', 'ACK', 'ACS', 'PTA', 'PTS_0', 'PTS_1', 'PTS_2',
  'PTS_3', 'PTS_4', 'CYTBO', 'SQR', 'NDHII', 'GROWTH', 'ATP_MAINTENANCE',
  'XCH_GLC', 'PIT', 'XCH_P', 'XCH_ACE1', '_ACE_OUT', 'XCH_ACE2', 'GL6P_HYDROLYSIS']

# REDEFINE LIST OF REACTIONS VARIED IN DE HERE
#listofreactions = reacVmaxes
listofreactions = ["ACEA", "ACEB", "ACK", "ACN_1", "ACN_2", "ACS", "ATP_syn",
"CITRA_SYN", "CYTBO", "EDA", "EDD", "ENO", "FBA", "FBP", "FUMA", "GDH", "GLT",
"GND", "GPM", "LPD", "MAD", "MDH", "MQO", "PCK", "PDH", "PFK", "PGI", "PGK",
 "PGL", "PIT", "PPC", "PPS", "PTA", "PYK", "RPE", "RPI", "SDH", "SK", "SQR",
  "TPI", "ZWF"]

# REDEFINE NUMBER OF TUPLES (couples, triples...) HERE
n = 41

# REDEFINE Vmax RANGE HERE
boundsrel = [(0.5, 10.0)] * n

# DIFFERENTIAL EVOLUTION PARAMETERS
mut=0.6876 # F
crossp=0.9784 # CR
popsize=48 # NP
its=5 # generations

# SIMULATION PARAMETERS
timestart = 0
timeend = 7200
rrpoints = 100

# Reads wild-type model
reader = libsbml.SBMLReader()
document = reader.readSBMLFromFile("../kinetic/E_coli_Millard2016_CITRA.xml")
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


def choose(mylist, n):
    return list(itertools.combinations(mylist, n))

def flux(reacid, r, x):
    # k is the reaction whose flux we care about - as string (2018-07-25)
    # r is a list of n reactions
    # xpy is a numpy axis with n elements
    for i in range(n):
        setVmax(r[i], x[i])

    # Simulate
    rr = roadrunner.RoadRunner(libsbml.writeSBMLToString(document))
    result = rr.simulate(timestart, timeend, rrpoints)

    k = rr.model.getReactionIds().index(reacid)
    return rr.model.getReactionRates()[k]

generations = []

# DE algorithm adapted from Pablo R Mier
def de(fobj, bounds, mut=mut, crossp=crossp, popsize=popsize, its=its):
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

with open('jmax.txt', 'r') as fobj:
    mystr = fobj.readline()
    jdx = int(mystr.strip('\n'))

if jdx >= len(allreactions):
    print('DONE')
    sys.exit(1)

else:
    print(jdx)
    print(allreactions[jdx])
    # main
    for combo in combolist:
        start_time = time.time() # time tracking

        print(combo)
        VmaxI = [wtVmaxes[i] for i in combo]
        print(VmaxI)
        bounds = (VmaxI*(boundsrel.T)).T
        print(bounds)

        reacid = allreactions[jdx]

        def fobj(x):
            return -flux(reacid, combo, x) # MAXIMUM

        # computation
        result = list(de(fobj, bounds))

        elapsed_time = time.time() - start_time
        print(elapsed_time) # time tracking

        # reassigns Vmaxes
        for i in range(len(combo)):
            setVmax(combo[i], VmaxI[i])

        # printing/writing results
        print(result[-1])
        with open('fluxdemax.txt', 'a') as f:
            f.write(str(allreactions[jdx]) + '\n')
            f.write(str(combo) + '\n')
            f.write(str(result[-1][0]) + '\n')
            f.write('Fitness: ' + str(result[-1][1]) + '\n')

        # write the index of reaction to file
        jdx += 1
        with open('jmax.txt', 'w') as fobj:
            fobj.write(str(jdx))


        sys.exit(2)
