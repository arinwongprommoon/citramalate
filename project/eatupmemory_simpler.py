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

sbmlfile = "E_coli_Millard2016_CITRA.xml"

for i in range(0,10000):

    if i % 100 == 0:
        print(i)

    reader = libsbml.SBMLReader()
    document = reader.readSBMLFromFile(sbmlfile)
    model = document.getModel()

    out = 0
    selection = ["CITRA", "iGROWTH'"]

    pointer = libsbml.writeSBMLToString(document)
    rr = roadrunner.RoadRunner(pointer)
    rr.timeCourseSelections = selection
    result = rr.simulate(0, 7200, 10)
    Y_PS = (result[-1,selection.index("CITRA")]*mmCITRA)/(0.25*7200*mmGLC)
    mu = result[-1,selection.index("iGROWTH'")]*3600

    gc.collect()
