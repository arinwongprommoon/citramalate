#! /usr/bin/env python

from __future__ import division, print_function
import roadrunner
import libsbml

import gc

mmCITRA = 146.098 # g/mol # Molecular mass of Citramalate https://pubchem.ncbi.nlm.nih.gov/compound/5460281
mmGLC = 180.156 # g/mol  # Molecular mass of Glucose https://pubchem.ncbi.nlm.nih.gov/compound/79025

class ecolicit:
    """SBML model of E. coli plust reaction for citramalate synthesis"""

    def __init__(self, sbmlfile = "E_coli_Millard2016_CITRA.xml"):
        """
        sbmlfile: SBML file with the model
        Vmax: Vmax of the citramalate synthesis reaction (mM/s)
        Km: Km of the citramalate synthesis reaction (mM)
        include_CITRA: Include citramalate species in the model (notice that it grows monotonically and hence it should not be included for steady state analysis)
        initial_CITRA: Initial concentration of citramalate (mM)
        """
        reader = libsbml.SBMLReader()
        self.document = reader.readSBMLFromFile(sbmlfile)
        self.model = self.document.getModel()
        self.getVmaxes()

    def writeToFile(self, filename = "test.xml"):
        libsbml.writeSBMLToFile(self.document, filename)

    def setVmax(self, reacId, value):
        # Units: mM/s
        reac = self.model.getReaction(reacId)
        lo = reac.getListOfAllElements()
        for el in lo:
            if type(el) is libsbml.Parameter:
                if el.getId() == 'Vmax':
                    el.setValue(value)

    def getVmax(self, reacId):
        # get Vmax of reaction reacId. See setVmax for info about units
        reac = self.model.getReaction(reacId)
        lo = reac.getListOfAllElements()
        for el in lo:
            if type(el) is libsbml.Parameter:
                if el.getId() == 'Vmax':
                    return el.getValue()

    def getVmaxes(self):
        # get Vmaxes of reactions. See setVmax for info about units
        Vmaxes = {}
        for reac in self.model.getListOfReactions():
            vm = self.getVmax(reac.id)
            if vm:
                Vmaxes[reac.id] = vm
        self.reacVmaxes = sorted(Vmaxes) # ids of reactions sorted alphabetically that have Vmax
        self.iniVmaxes = [Vmaxes[r] for r in self.reacVmaxes] # initial values of Vmax (as in the kinetic model)

if __name__ == '__main__':
    ecit = ecolicit()
    for i in range(0,10000):
        if i % 50 == 0:
            print(i)
        if i % 2 == 0:
            ecit.setVmax('ACEA', 0.1234)
        else:
            ecit.setVmax('ACEA', 0.5647)
        #modelstring = libsbml.writeSBMLToString(ecit.document)
        ecit.writeToFile()
