#!/usr/bin/python
# Transform cobra model to Flexible Net (FN)
from __future__ import division, print_function
import csv
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

test = 'all'
#knockouts = ['b1380', 'b0903'] # list of genes to knock out
#knockouts = ['b0903'] # list of genes to knock out
knockouts = []

if test == 'core':
    modelfile = "ecoli_core_model"
    objective = 'Biomass_Ecoli_core_w_GAM'
elif test == 'all':
    modelfile = "MODEL1108160000"
#        objective = 'Ec_biomass_iJO1366_core_53p95M'
    objective = 'CitraSink'

model = cobra.io.read_sbml_model(modelfile+'.xml')
addCimA(model)
print('Reactions:', len(model.reactions),
      '; Metabolites', len(model.metabolites),
      '; Genes:', len(model.genes))

print('Cobra results before change')
model.objective = objective
if knockouts:
    # WARNING: only the first gene in knockouts is deleted
    sol, sta = single_gene_deletion(model,[model.genes.get_by_id(knockouts[0])])
    print('Status:', sta, 'Solution:', sol)
else:
    solution = model.optimize()
    print('Status:', solution.status, '; Solution:', solution.objective_value)

# Reads in bounds specified in CSV, stores in boundslist
with open('boundaries/OldBoundaries.csv', 'rt') as fobj:
    reader = csv.reader(fobj)
    boundslist = list(reader)
    boundslist = boundslist[1:] # removes header

# Stores default bounds in originalbounds
originalbounds = {}
for row in boundslist:
    reac = model.reactions.get_by_id(row[0])
    original_lower_bound = reac.lower_bound
    original_upper_bound = reac.upper_bound
    originalbounds[reac] = (original_lower_bound, original_upper_bound)

# Sets bounds according to CSV    
for row in boundslist:
    reac = model.reactions.get_by_id(row[0])
    reac.lower_bound = float(row[1])
    reac.upper_bound = float(row[2])

# Reverts bounds of each reaction in turn
for row in boundslist:
    reac = model.reactions.get_by_id(row[0])
    reac.lower_bound = originalbounds[reac][0]
    reac.upper_bound = originalbounds[reac][1]
    print('%s bounds reverted' % row[0])
    print('Cobra results after change')
    model.objective = objective
    if knockouts:
        # WARNING: only the first gene in knockouts is deleted
        sol, sta = single_gene_deletion(model,[model.genes.get_by_id(knockouts[0])])
        print('Status:', sta, 'Solution:', sol)
    else:
        solution = model.optimize()
        print('Status:', solution.status, '; Solution:', solution.objective_value)
    reac.lower_bound = float(row[1])
    reac.upper_bound = float(row[2])
