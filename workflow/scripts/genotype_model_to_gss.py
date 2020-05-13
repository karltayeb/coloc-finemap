import pickle
from coloc.independent_model_ss import IndependentFactorSER as GSS
from coloc.misc import *

print('loading data')
data = make_gtex_genotype_data_dict(snakemake.input.expression, snakemake.input.genotype)

# load and rehydrate model
model = pickle.load(open(snakemake.input.model, 'rb'))
rehydrate_model(model)

# instantiate gss model
K = model.dims['K']
gss = GSS(**data, K=model.dims['K'])
gss.prior_activity = np.ones(K) * snakemake.params.pi

# initialize with genotype model
model.__dict__.pop('X', None)
model.__dict__.pop('Y', None)
model.__dict__.pop('Y', None)
model.__dict__.pop('records', None)
model.__dict__.pop('precompute', None)
gss.__dict__.update(model.__dict__)

# fit with spike and slab
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

# compute recrods
compute_records(gss)

# save variant report
base_path = snakemake.output[0][:-len('.model')]
variant_report_path = '{}.variants.bed'.format(base_path)
print('saving variant report to: {}'.format(variant_report_path))
try:
    gene = base_path.split('/')[-2]
    make_variant_report(model, gene).to_csv(variant_report_path, sep='\t')
except Exception:
    print('There was an error generating secondary files')

#save model
print('saving model')
strip_and_dump(gss, snakemake.output[0])