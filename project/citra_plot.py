#!/usr/bin/env python

from __future__ import division, print_function
from ecolicitra import ecolicit, mmCITRA, mmGLC
import matplotlib.pyplot as plt
import numpy as np
import roadrunner
import libsbml

# Parameters of kinetic model
include_CITRA = True
ecit = ecolicit(include_CITRA = include_CITRA)
ecit.setVmax('CITRA_SYN', 4.0)

ecit.time0 = 0
ecit.timef = 2*3600 # final simulation time in seconds
ecit.npoints = 100 # Number of points to be computed in the simulation

wildprod = ecit.comproducti()

print("Final simulation time: ", ecit.timef)
print("Number of points: ", ecit.npoints)
print("Productivity (wild type): ", wildprod, " h-1")

ecit.plot("CITRA")
