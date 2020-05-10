import pandas as pd
import pickle
import numpy as np
from coloc.independent_model_ss import IndependentFactorSER as GSS
from coloc.misc import component_scores, make_variant_report, strip_and_dump

data = pickle.load(open(snakemake.input[0], 'rb'))
model = GSS(**data, K=snakemake.params.k)

fit_args = {
    'max_iter': 300,
    'update_covariate_weights': False,
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


model.prior_activity = np.ones(K) * 0.1
print('fitting full model')
print(model.dims)
fit_args = {
    'max_iter': 300,
    'update_covariate_weights': False,
    'update_active': True,
    'update_weights': True,
    'update_pi': True,
    'ARD_weights': True,
    'update_variance': True,
    'verbose': False
}
model.fit(**fit_args)
print('model fit:\n\titers:{}\n\tELBO:{}\n\trun-time:{}'.format(len(model.elbos), model.elbos[-1], model.run_time))
strip_and_dump(model, snakemake.output.model, False)


base_path = snakemake.output[0][:-len('.model')]
print('generating scores and variant file')
try:
    gene = base_path.split('/')[-2]
    make_variant_report(model, gene).to_csv('{}.variants.bed'.format(base_path), sep='\t')
except Exception:
    print('There was an error generating secondary files')

print('saving model')
strip_and_dump(model, snakemake.output[0])
