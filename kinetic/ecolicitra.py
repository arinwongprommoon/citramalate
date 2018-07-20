#! /usr/bin/env python
#Class for SBML kinetic model of E. coli plus reaction for the production of citramalate
from __future__ import division, print_function
import roadrunner
import libsbml
import matplotlib.pyplot as plt

mmCITRA = 146.098 # g/mol # Molecular mass of Citramalate https://pubchem.ncbi.nlm.nih.gov/compound/5460281
mmGLC = 180.156 # g/mol  # Molecular mass of Glucose https://pubchem.ncbi.nlm.nih.gov/compound/79025

class ecolicit:
    """SBML model of E. coli plust reaction for citramalate synthesis"""

    def __init__(self, sbmlfile = "E_coli_Millard2016.xml", Vmax = 4.0, Km = 0.495, include_CITRA = True, initial_CITRA = 0.0):
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
#        print("Number of species: ", self.model.getNumSpecies())
#        print("Number of reactions: ", self.model.getNumReactions())

        # Create species that models growth. The growth rate will be equal to
        # the derivative of the species concentration
        igrowth = self.model.createSpecies()
        igrowth.setId('iGROWTH')
        igrowth.setName('iGROWTH')
        igrowth.setCompartment('cell')
        igrowth.setConstant(False)
        igrowth.setInitialConcentration(0.0)
        igrowth.setBoundaryCondition(False)
        igrowth.setHasOnlySubstanceUnits(False)

        growthr = self.model.getReaction('GROWTH')
        spig = growthr.createProduct()
        spig.setSpecies("iGROWTH")

        # Create citramalate species
        if include_CITRA:
            citra = self.model.createSpecies()
            citra.setId('CITRA')
            citra.setName('CITRA')
            citra.setCompartment('cell')
            citra.setConstant(False)
            citra.setInitialConcentration(initial_CITRA)
            citra.setBoundaryCondition(False)
            citra.setHasOnlySubstanceUnits(False)

        # Create reaction for citramalate synthesis
        citrasyn = self.model.createReaction()
        citrasyn.setId("CITRA_SYN")
        citrasyn.setName("CITRA_SYN");

        spr1 = citrasyn.createReactant()
        spr1.setSpecies("ACCOA")
        spr2 = citrasyn.createReactant()
        spr2.setSpecies("PYR")
        spr3 = citrasyn.createReactant()
        spr3.setSpecies("H2O")
        spp1 = citrasyn.createProduct()
        spp1.setSpecies("COA")
        spp2 = citrasyn.createProduct()
        spp2.setSpecies("Hin")
        if include_CITRA:
            spp3 = citrasyn.createProduct()
            spp3.setSpecies("CITRA")

        kl = citrasyn.createKineticLaw()
        para = kl.createParameter()
        para.setId("Vmax")
        para.setValue(Vmax)
        para = kl.createParameter()
        para.setId("Km")
        para.setValue(Km)
        math_ast = libsbml.parseL3Formula('(Vmax * ACCOA) / (ACCOA + Km)')
        kl.setMath(math_ast)

        # Vmaxes of the reactions
        self.getVmaxes()

    def writeToFile(self, filename = "E_coli_Millard2016_CITRA.xml"):
        libsbml.writeSBMLToFile(self.document, filename)

    def setFEED(self, value):
        # value: mmol/s of glucose fed into the extracellular compartment.
        # Notice that the volume of the extracellular compartment is 100l, thus
        # the increase in concentration of glucose produced by the feed in the
        # extracellular glucose is value/100 mM/s
        self.model.getParameter('FEED').setValue(value)

    def getFEED(self):
        # see setFEED for info about units
        return self.model.getParameter('FEED').getValue()

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

    def comproducti(self, tol=99999):
        """
            Computes steady state productivity.
            Argument:
                tol = 'epsilon' value to check the maximum absolute value of
                      floating species concentration among all species in the
                      model against to determine if the system has reached
                      steady state or not

                      If the system has not reached steady state, this function
                      will return -1e-4 instead of the productivity value.
        """
        selection = ["CITRA", "iGROWTH'"]
        rr = roadrunner.RoadRunner(libsbml.writeSBMLToString(self.document))
        rr.timeCourseSelections = selection
        result = rr.simulate(self.time0, self.timef, self.npoints)
        # -2: removes GROWTH and CITRA from the list because they aren't steady
        # state anyway
        st = max(abs(rr.model.getFloatingSpeciesConcentrationRates())[:-2])
        if steadystate == True:
            return st
        else:
            if st < tol:
                Y_PS = (result[-1,selection.index("CITRA")]*mmCITRA)/(self.getFEED()*self.timef*mmGLC)
                mu = result[-1,selection.index("iGROWTH'")]*3600
                return mu*Y_PS
            else:
                return -1e-4

    def compsteady(self):
        """
            Computes maximum absoute value of floating species concentration
            among all species in the model - useful in checking if system has
            reached steady state or not
        """
        selection = ["CITRA", "iGROWTH'"]
        rr = roadrunner.RoadRunner(libsbml.writeSBMLToString(self.document))
        rr.timeCourseSelections = selection
        result = rr.simulate(self.time0, self.timef, self.npoints)
        # -2: removes GROWTH and CITRA from the list because they aren't steady
        # state anyway
        st = max(abs(rr.model.getFloatingSpeciesConcentrationRates())[:-2])
        return st

    def plot(self, species):
        """
            Plots concentration of specified species over time course specified
            Input can also be any other SBML value that can be selected in
            RoadRunner.
        """
        rr = roadrunner.RoadRunner(libsbml.writeSBMLToString(self.document))
        rr.timeCourseSelections = [species]
        result = rr.simulate(self.time0, self.timef, self.npoints)
        plt.plot(result[:,result.colnames.index(species)])
        plt.show()
