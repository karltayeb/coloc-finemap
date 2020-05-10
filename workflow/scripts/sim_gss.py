import numpy as np
import pandas as pd
from coloc.independent_model_ss import IndependentFactorSER as GSS
import scipy as sp
import pickle

def compute_records(model):
    """
    save the model with data a weight parameters removed
    add 'mini_weight_measn' and 'mini_weight_vars' to model
    the model can be reconstituted in a small number of iterations
    """
    credible_sets, purity = model.get_credible_sets(0.999)
    active = np.array([purity[k] > 0.1 for k in range(model.dims['K'])])

    try:
        snps = np.unique(np.concatenate([
            credible_sets[k] for k in range(model.dims['K']) if active[k]]))
    except Exception:
        snps = np.unique(np.concatenate([
            credible_sets[k][:5] for k in range(model.dims['K'])]))
    mask = np.isin(model.snp_ids, snps)

    wv = model.weight_vars[:, :, mask]
    wm = model.weight_means[:, :, mask]

    records = {
        'active': active,
        'purity': purity,
        'credible_sets': credible_sets,
        'EXz': model.pi @ model.X,
        'mini_wm': wm,
        'mini_wv': wv,
        'snp_subset': mask
    }
    model.records = records


def strip_and_dump(model, path, save_data=False):
    """
    save the model with data a weight parameters removed
    add 'mini_weight_measn' and 'mini_weight_vars' to model
    the model can be reconstituted in a small number of iterations
    """
    compute_records(model)

    # purge precompute
    for key in model.precompute:
        model.precompute[key] = {}
    model.__dict__.pop('weight_means', None)
    model.__dict__.pop('weight_vars', None)
    if not save_data:
        model.__dict__.pop('X', None)
        model.__dict__.pop('Y', None)
        model.__dict__.pop('covariates', None)
    pickle.dump(model, open(path, 'wb'))

# load data
data = pickle.load(open(snakemake.input.data, 'rb'))
info = pickle.load(open(snakemake.input.info, 'rb'))

##### TRAIN CAFEH GENOTYPE
K = np.max([10, info['causal_snps'].size * 2])
model = GSS(**data, K=K)
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
