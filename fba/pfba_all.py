#!/usr/bin/env python3
# Uses bondaries specified in CSV file to set boundaries for FBA on the
# stoichiometric model. Finds optimal solution (maximise sense) cycling through
# all objective functions possible
from __future__ import division, print_function
import csv
import cobra.test
from cobra import Reaction, Metabolite

# DEFINE BOUNDARY FILE HERE
boundariesfile = '41dBoundaries.csv' # NAME OF FILE
fileloc = 'boundaries_citra/'+boundariesfile # DIRECTORY

modelfile = "MODEL1108160000" # DEFINE MODEL FILE HERE

 # DEFINE OBJECTIVE SENSE HERE 'minimize' or 'maximize'
objective_sense = 'maximize'

def addCimA(model):
    """Add CimA reaction and sink for citramalate to cobra model"""
    reaccima = Reaction('CIMA')
    reaccima.name = '(R)-Citramalate production'
    reaccima.lower_bound = 0.0
    reaccima.upper_bound = 1000.0

    """CIMA reaction"""
    pyr_c = model.metabolites.get_by_id("pyr_c") # Pyruvate
    accoa_c = model.metabolites.get_by_id("accoa_c") # Acetyl-CoA
    h2o_c = model.metabolites.get_by_id("h2o_c") # H2O
    rcitramalate_c = Metabolite(
        'citramalate_c',
        formula='C5H6O5',
        name='(R)-citramalate',
        charge=-2,
        compartment='c')
    coa_c = model.metabolites.get_by_id("coa_c") # CoA

    reaccima.add_metabolites({pyr_c: -1.0,
                              accoa_c: -1.0,
                              h2o_c: -1.0,
                              rcitramalate_c: 1.0,
                              coa_c: 1.0})
    reaccima.gene_reaction_rule = 'CimA37'
#    print(reaccima.reaction)
#    print(reaccima.genes)
    model.add_reaction(reaccima)
    reaccima.objective_coefficient = 0.0

    """Sink for Citramalate"""
    reaccisink = Reaction('CitraSink')
    reaccisink.name = 'Sink needed to allow (R)-Citramalate to leave the system'
    reaccisink.lower_bound = 0.0
    reaccisink.upper_bound = 1000.0

    reaccisink.add_metabolites({rcitramalate_c: -1.0})
#    print(reaccisink.reaction)
#    print(reaccisink.genes)
    model.add_reaction(reaccisink)
    reaccisink.objective_coefficient = 0.0

model = cobra.io.read_sbml_model(modelfile+'.xml')
addCimA(model)

# Loops through all 2,585 reactions - prints objective values to CSV
# (here as a control)
print('Cobra results before change')
with open('Objectives_pFBAdefault.csv', 'w') as fobj:
    writer = csv.writer(fobj)
    for reaction in model.reactions:
        model.objective = reaction.id
        solution = model.optimize(objective_sense=objective_sense)
        pfba_solution = cobra.flux_analysis.pfba(model)
        print(reaction.id, '; Status:', solution.status, '; FBA Solution:', solution.objective_value, '; pFBA Solution:', pfba_solution.objective_value)
        writer.writerow([reaction.id, solution.objective_value, pfba_solution.objective_value])

### Reads CSV file listing reactions and intended lower and upper bounds
with open(fileloc, 'rt') as fobj:
    reader = csv.reader(fobj)
    boundslist = list(reader)
    boundslist = boundslist[1:] # removes header
    for row in boundslist:
        reac = model.reactions.get_by_id(row[0])
        reac.lower_bound = float(row[1])
        reac.upper_bound = float(row[2])

print('Bounds changed')

# Loops through all 2,585 reactions - prints objective values to CSV
filename = 'Objectives_bounds-pFBA' + boundariesfile
with open(filename, 'w') as fobj:
    writer = csv.writer(fobj)
    for reaction in model.reactions:
        model.objective = reaction.id
        solution = model.optimize(objective_sense=objective_sense)
        pfba_solution = cobra.flux_analysis.pfba(model)
        #print(reaction.id, '; Status:', solution.status, '; Solution:', solution.objective_value)
        writer.writerow([reaction.id, solution.objective_value, pfba_solution.objective_value])
