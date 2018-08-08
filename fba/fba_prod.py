#!/usr/bin/env python3
# Transform cobra model to Flexible Net (FN)
from __future__ import division, print_function
import csv
import pandas
import cobra.test
from cobra import Reaction, Metabolite
from cobra.flux_analysis import (
    single_gene_deletion, single_reaction_deletion, double_gene_deletion,
    double_reaction_deletion)

def addCimA(model):
    """Add CimA reaction and sink for citramalate to cobra model model"""
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

#########  MAIN

modelfile = "MODEL1108160000"
knockouts = []

model = cobra.io.read_sbml_model(modelfile+'.xml')
addCimA(model)
print('Reactions:', len(model.reactions),
      '; Metabolites', len(model.metabolites),
      '; Genes:', len(model.genes))
      
# Re-defining objective function to be the product of flux through citramalate
# and the growth reaction
model.solver = 'scipy' # can't do it - must use the commercial 'cplex'
product_objective = model.problem.Objective(
    model.reactions.get_by_id('Ec_biomass_iJO1366_core_53p95M').flux_expression * 
    model.reactions.get_by_id('CitraSink').flux_expression, direction = 'max')
model.objective = product_objective

print('Cobra results before change')
if knockouts:
    # WARNING: only the first gene in knockouts is deleted
    sol, sta = single_gene_deletion(model,[model.genes.get_by_id(knockouts[0])])
    print('Status:', sta, 'Solution:', sol)
else:
    solution = model.optimize(objective_sense='maximize')
    print('Status:', solution.status, '; Solution:', solution.objective_value)

### Reads CSV file listing reactions and intended lower and upper bounds
with open('boundaries_citra/41dBoundaries.csv', 'rt') as fobj:
    reader = csv.reader(fobj)
    boundslist = list(reader)
    boundslist = boundslist[1:] # removes header
    for row in boundslist:
        reac = model.reactions.get_by_id(row[0])
        reac.lower_bound = float(row[1])
        reac.upper_bound = float(row[2])

print('Bounds changed')

print('Cobra results after change')
if knockouts:
    # WARNING: only the first gene in knockouts is deleted
    sol, sta = single_gene_deletion(model,[model.genes.get_by_id(knockouts[0])])
    print('Status:', sta, 'Solution:', sol)
else:
    solution = model.optimize(objective_sense='maximize')
    print('Status:', solution.status, '; Solution:', solution.objective_value)

f = solution.fluxes
output = f[f != 0]
output.to_csv('FluxesAfterBound.csv')
print('Output to CSV')

print('Model summary.....')
model.summary()
print('Citramalate summary....')
model.metabolites.citramalate_c.summary()
print('Flux through growth....')
print(f[7]) # Ec_biomass_iJO1366_core_53p95M is the 8th reaction
