import pickle
from coloc.independent_model_ss import IndependentFactorSER as GSS
from coloc.misc import *

data = make_gtex_genotype_data_dict(snakemake.input.expression, snakemake.input.genotype)
model = pickle.load(open(snakemake.input.model, 'rb'))

K = model.dims['K']
gss = GSS(**data, K=model.dims['K'])
gss.prior_activity = np.ones(K) * 0.01

model.__dict__.pop('X', None)
model.__dict__.pop('Y', None)
model.__dict__.pop('Y', None)
model.__dict__.pop('records', None)
model.__dict__.pop('precompute', None)
gss.__dict__.update(model.__dict__)

fit_args = {
    'max_iter': 300,
    'update_covariate_weights': True,
    'update_active': True,
    'update_weights': True,
    'update_pi': True,
    'ARD_weights': True,
    'update_variance': True,
    'verbose': True
}
print('training model')
for arg in fit_args:
    print('\t{}: {}'.format(arg, fit_args[arg]))
gss.fit(**fit_args)

print('saving model')
strip_and_dump(gss, snakemake.output[0])