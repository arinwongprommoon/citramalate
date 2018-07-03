#!/usr/bin/env python
# Adds citramalate reaction to kinetic model
# (this script only used once when E_coli_Millard2016_CITRA.xml didn't exist)

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
ecit.writeToFile()
