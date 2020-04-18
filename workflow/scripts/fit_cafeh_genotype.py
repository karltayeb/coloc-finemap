import pandas as pd
import pickle
import numpy as np
from coloc.independent_model2 import IndependentFactorSER
from coloc.misc import component_scores, make_variant_report, strip_and_dump

data = pickle.load(open(snakemake.input[0], 'rb'))
model = IndependentFactorSER(**data, K=snakemake.params.k)

fit_args = {
    'max_iter': 300,
    'update_covariate_weights': True,
    'update_weights': True,
    'update_pi': True,
    'ARD_weights': True,
    'update_variance': True,
    'verbose': True
}

print('training model')
for arg in fit_args:
    print('\t{}: {}'.format(arg, fit_args[arg]))
model.fit(**fit_args)

base_path = snakemake.output[0][:-len('.model')]

print('generating scores and variant file')
try:
    component_scores(model).to_json('{}.scores'.format(base_path))
    gene = base_path.split('/')[-2]
    make_variant_report(model, gene).to_csv('{}.variants.bed'.format(base_path), sep='\t')
except Exception:
    print('There was an error generating secondary files')

print('saving model')
strip_and_dump(model, snakemake.output[0])
