import numpy as np
import pandas as pd
from coloc.independent_model2 import IndependentFactorSER as M
from coloc.independent_model_summary import IndependentFactorSER as M2
from coloc.misc import *
import scipy as sp

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


def strip(model, save_data=False):
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

# load data
data = pickle.load(open(snakemake.input.data, 'rb'))

fit_args = {
    'max_iter': 20,
    'update_covariate_weights': True,
    'update_weights': True,
    'update_pi': True,
    'ARD_weights': True,
    'update_variance': True,
    'verbose': False
}

susies = {}

for t in range(10):
    print('training SuSiE tissue: {}'.format(t))
    data_t = {
        'Y': data['Y'][t][None],
        'X': data['X'],
        'snp_ids': data['snp_ids']
    }
    susie_t = M(**data_t, K=10)

    fit_args = {
        'max_iter':100,
        'ARD_weights':False,
        'verbose': False
    }
    susie_t.tolerance = 1e-8
    print(susie_t.dims)
    susie_t.fit(**fit_args)
    print('model fit:\n\titers:{}\n\tELBO:{}\n\trun-time:{}'.format(
        len(susie_t.elbos), susie_t.elbos[-1], susie_t.run_time))

    compute_records(susie_t)
    strip(susie_t)
    susies[t] = susie_t

print('saving susie')
pickle.dump(susies, open(snakemake.output.susie, 'wb'))

