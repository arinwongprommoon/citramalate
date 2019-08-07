#! /usr/bin/env python
#Test kinetic model of E. coli plus citramalate reaction with different glucose feeds
from __future__ import division, print_function
from ecolicitra import ecolicit, mmCITRA, mmGLC
import roadrunner
import libsbml

# Parameters of kinetic model
include_CITRA = True
ecit = ecolicit(include_CITRA = include_CITRA)
ecit.setVmax('CITRA_SYN', 4.0)

ecit.time0 = 0
ecit.timef = 2*3600 # final simulation time in seconds
ecit.npoints = 100 # Number of points to be computed in the simulation

ecit.setFEED(1.74)
print(ecit.getFEED())

rr = roadrunner.RoadRunner(libsbml.writeSBMLToString(ecit.document))
result = rr.simulate(ecit.time0, ecit.timef, ecit.npoints)
print("Concentrations")
for nospec, spec in enumerate(rr.model.getFloatingSpeciesIds()):
    print(spec, ":", rr.model.getFloatingSpeciesConcentrations()[nospec])

wildprod = ecit.comproducti()
print("Productivity: ", wildprod, " h-1")

y = ecit.comyield()
print("Citramalate yield: ", y, " (dimensionless)")
