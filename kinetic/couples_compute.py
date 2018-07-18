#!/usr/bin/env python

# Runs computations for couples of reactions from a specified list, writes data
# to numpy array files.

# Using couples_noloop, this program stops after each pair and continues on to
# the next pair when executed again. To run through all pairs, execute
# ./couples_run (couples_run.sh)

from couples_noloop import Coupl
import sys

#Run simulations, writes data both into TXT and NPZ
run = Coupl(xstart=0.1, xend=10.0, ystart=0.1, yend=10.0, points=10) # define vmax ranges here
run.setVmax('CITRA_SYN', 4.0)
run.time0 = 0
run.timef = 2*3600 # final simulation time in seconds
run.npoints = 100 # Number of points to be computed in the simulation
#run.listofreactions = ['CITRA_SYN', 'GLT', 'LPD'] # define list of reactions
run.listofreactions = ['CITRA_SYN', 'GLT', 'LPD', 'ATP_MAINTENANCE', 'GDH', 'ATP_syn', 'ACEA', 'ZWF']

run.writeFileHeader()
out = run.computeCouples(writemethod='steadystate')
sys.exit(out)
