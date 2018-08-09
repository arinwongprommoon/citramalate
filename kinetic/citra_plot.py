#!/usr/bin/env python
# Demonstrate and use the plot function in ecolicitra

from __future__ import division, print_function
from ecolicitra import ecolicit, mmCITRA, mmGLC

# Parameters of kinetic model
include_CITRA = True
ecit = ecolicit(include_CITRA = include_CITRA)
ecit.setVmax('CITRA_SYN', 4.0)
ecit.time0 = 0
ecit.timef = 2*3600 # final simulation time in seconds
ecit.npoints = 100 # Number of points to be computed in the simulation
ecit.plot("FBA") # Change ID to change the reaction being plotted
