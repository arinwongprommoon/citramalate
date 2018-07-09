import roadrunner
import libsbml
import gc

mmCITRA = 146.098 # g/mol # Molecular mass of Citramalate https://pubchem.ncbi.nlm.nih.gov/compound/5460281
mmGLC = 180.156 # g/mol  # Molecular mass of Glucose https://pubchem.ncbi.nlm.nih.gov/compound/79025

reader = libsbml.SBMLReader()
document = reader.readSBMLFromFile("E_coli_Millard2016_CITRA.xml")
model = document.getModel()
sbmlstring = libsbml.writeSBMLToString(document)

def setVmax(reacId, value):
    # Units: mM/s
    reac = model.getReaction(reacId)
    lo = reac.getListOfAllElements()
    for el in lo:
        if type(el) is libsbml.Parameter:
            if el.getId() == 'Vmax':
                el.setValue(value)

for i in range(0,100):
    gc.collect()
    print(i)
    setVmax('GLT', 0.1)
    sbmlstring = libsbml.writeSBMLToString(document)
    rr = roadrunner.RoadRunner(sbmlstring)
    selection = ["CITRA", "iGROWTH'"]
    rr.timeCourseSelections = selection
    result = rr.simulate(0, 7200, 100)
    Y_PS = (result[-1,selection.index("CITRA")]*mmCITRA)/(0.23*7200*mmGLC)
    mu = result[-1,selection.index("iGROWTH'")]*3600
    productivity = mu*Y_PS
    print(productivity)
