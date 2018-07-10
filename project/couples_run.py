#!/usr/bin/env python
import os
import couples_noloop
import couples_compute

with open('couplesi.txt', 'w') as fobj:
    fobj.write('0')

while True:
    couples_compute.ccrun()
