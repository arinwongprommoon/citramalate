#!/usr/bin/env python

from __future__ import division, print_function
from ecolicitra_copy import ecolicit, mmCITRA, mmGLC
import numpy as np
import itertools
import roadrunner
import libsbml

import time
import gc
import random

include_CITRA = True
ecit = ecolicit(include_CITRA = include_CITRA)
ecit.setVmax('CITRA_SYN', 4.0)
ecit.time0 = 0
ecit.timef = 2*3600 # final simulation time in seconds
ecit.npoints = 10 # Number of points to be computed in the simulation

for i in range(0,1000):
    if i % 10 == 0:
        print(i)
    #wildprod = 0
    wildprod = ecit.comproducti()
    gc.collect()
