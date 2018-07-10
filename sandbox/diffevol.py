#!/usr/bin/env python3
# Based on Pablo R Mier
# https://pablormier.github.io/2017/09/05/a-tutorial-on-differential-evolution-with-python/

import numpy as np

def fobj(x):
  value = 0
  for i in range(len(x)):
      value += x[i]**2
  return value / len(x)

bounds = [(-100, 100)]*32

def de(fobj, bounds, mut=0.8, crossp=0.7, popsize=20, its=3000):
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

result = list(de(fobj, bounds))
print(result[-1])
