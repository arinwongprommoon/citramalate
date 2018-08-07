#!/usr/bin/env python

# Runs computations for couples of reactions from a specified list, writes data
# to numpy array files.

# This program stops after each pair and continues on to
# the next pair when executed again. To run through all pairs, execute
# ./couples_run (couples_run.sh). In a nutshell, this script is not very useful
# on its own

from couples import Coupl
import sys

#Run simulations, writes data both into TXT and NPZ
run = Coupl(xstart=0.1, xend=10.0, ystart=0.1, yend=10.0, points=10) # define vmax ranges here
run.setVmax('CITRA_SYN', 4.0)
run.time0 = 0
run.timef = 2*3600 # final simulation time in seconds
run.npoints = 100 # Number of points to be computed in the simulation
run.tolerance = 1e-6
#run.listofreactions = ['CITRA_SYN', 'GLT', 'LPD'] # define list of reactions
run.listofreactions = ['CITRA_SYN', 'GLT', 'LPD', 'GDH', 'ATP_syn', 'ACEA', 'PYK', 'ZWF', 'NDHII', 'MQO']
#run.listofreactions = ['CYTBO', 'MQO', 'MDH', 'ZWF', 'GLT', 'GDH', 'ATP_syn', 'ACK', 'ACEA', 'EDD']

run.writeFileHeader()
out = run.computeCouplesNoLoop(writemethod='writeToText')
sys.exit(out)
