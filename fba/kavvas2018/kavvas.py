#!/usr/bin/env python3
# Transform cobra model to Flexible Net (FN)
from __future__ import division, print_function
import glob, os
import csv
import pandas
import cobra.test
from cobra import Reaction, Metabolite
from cobra.flux_analysis import (
    single_gene_deletion, single_reaction_deletion, double_gene_deletion,
    double_reaction_deletion)
    
for modelfile in glob.glob("*.json"):
    model = cobra.io.load_json_model(modelfile)
    print(modelfile)
    print('Reactions:', len(model.reactions),
      '; Metabolites', len(model.metabolites),
      '; Genes:', len(model.genes))
    solution = model.optimize(objective_sense='maximize')
    print('Status:', solution.status, '; Solution:', solution.objective_value)
    
    csvname = modelfile.replace("json","csv")
    with open(csvname, 'w') as fobj:
        writer = csv.writer(fobj)
        for reac in model.reactions:
            writer.writerow([reac.id, reac.lower_bound, reac.upper_bound])
