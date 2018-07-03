#!/usr/bin/env python

from __future__ import division, print_function
from ecolicitra_copy import ecolicit, mmCITRA, mmGLC
import matplotlib.pyplot as plt
import numpy as np
import roadrunner
import libsbml

# Parameters of kinetic model
include_CITRA = True
ecit = ecolicit(include_CITRA = include_CITRA)
ecit.setVmax('CITRA_SYN', 5.0)

ecit.time0 = 0
ecit.timef = 2*3600 # final simulation time in seconds
ecit.npoints = 10000 # Number of points to be computed in the simulation

wildprod = ecit.comproducti()
ss1 = ecit.compsteadystate()
ss2 = ecit.altcompsteadystate()

print("Final simulation time: ", ecit.timef)
print("Number of points: ", ecit.npoints)
print("Productivity (wild type): ", wildprod, " h-1")
print("Steady state (method 1): ", ss1)
print("Steady state (method 2): ", ss2)
