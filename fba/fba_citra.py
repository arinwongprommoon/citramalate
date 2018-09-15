#!/usr/bin/env python3
# Uses bondaries specified in CSV file to set boundaries for FBA on the
# stoichiometric model. Finds optimal solution (maximise sense) using
# citramalate flux as the objective function, and prints out additional info
from __future__ import division, print_function
import csv
import cobra.test
from cobra import Reaction, Metabolite

# DEFINE FILE TO BE LOADED FOR FBA BOUNDS HERE
loc = 'result/obj/Objectives_bounds-41dBoundaries_ALL.csv'

modelfile = "MODEL1108160000" # DEFINE MODEL FILE HERE
objective = 'CitraSink' # DEFINE OBJECTIVE REACTION HERE
#objective = 'Ec_biomass_iJO1366_core_53p95M'

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

print('Cobra results before change')
model.objective = objective
solution = model.optimize(objective_sense='maximize')
print('Status:', solution.status, '; Solution:', solution.objective_value)
productivity = solution.fluxes[7]*solution.objective_value
print('Productivity: ', productivity)

# Reads CSV file listing reactions and intended lower and upper bounds
with open(loc, 'rt') as fobj:
    reader = csv.reader(fobj)
    boundslist = list(reader)
    boundslist = boundslist[1:] # removes header
    for row in boundslist:
        reac = model.reactions.get_by_id(row[0])
        reac.lower_bound = float(row[1])
        reac.upper_bound = float(row[2])
        
print('Bounds changed')

print('Cobra results after change')
solution = model.optimize(objective_sense='maximize')
print('Status:', solution.status, '; Solution:', solution.objective_value)
productivity = solution.fluxes[7]*solution.objective_value
print('Productivity: ', productivity)

f = solution.fluxes
output = f[f != 0]
output.to_csv('FluxesAfterBound.csv')
print('Output to CSV')

print('Model summary.....')
model.summary()
print('Citramalate summary....')
model.metabolites.citramalate_c.summary()

