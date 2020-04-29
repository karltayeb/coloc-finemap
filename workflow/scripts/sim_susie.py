import numpy as np
import pandas as pd
from coloc.independent_model2 import IndependentFactorSER as M
from coloc.independent_model_summary import IndependentFactorSER as M2
from coloc.misc import *
import scipy as sp

def make_simulation(genotype, T, pve, sparsity):
    # select snps
    NN = genotype.shape[1]
    start = np.random.choice(NN-1000)
    X = genotype.values.T[start:start+1000]
    snp_ids = genotype.columns.values[start:start+1000]

    N, M = X.shape
    causal_snps = np.random.choice(N, 10)
    true_effects = np.zeros((T, N))
    true_effects[:, causal_snps] = \
        np.random.binomial(1, sparsity, ((T, 10))) \

    if snakemake.params.sample_effects:
        true_effects = true_effects \
            * np.random.normal(size=true_effects.shape)
    tissue_variance = np.array([
        compute_sigma2(X, te, pve) for te in true_effects
    ])

    #simulate expression
    expression = (true_effects @ X) + \
        np.random.normal(size=(T, M)) * np.sqrt(tissue_variance)[:, None]
    expression = expression - expression.mean(1)[:, None]
    pos = np.array([int(x.split('_')[1]) for x in genotype.columns.values])

    # simulate expression 
    data = {
        'X': X,
        'Y': expression,
        'snp_ids': snp_ids,
    }

    sim_info = {
        'causal_snps': causal_snps[np.abs(true_effects[:, causal_snps]).sum(0) > 0],
        'true_effects': true_effects,
        'tissue_variance': tissue_variance,
        'expression': expression
    }
    return data, sim_info


def compute_sigma2(X, true_effect, pve):
    var = np.var(true_effect @ X)
    sigma2_t = var/pve - var
    if sigma2_t == 0:
        # if variance is 0, there were no causal variants-- dont care what the variance is
        sigma2_t = 1.0
    return sigma2_t


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


def strip(model):
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
data = pickle.load(open(snakemake.input.info, 'rb'))

### RUN SUSIE
print('training susie')
susie = M(**data, K=T*5)
a = sp.linalg.block_diag(*[np.ones((5)) * 1e10 for t in range(T)])
a = 1 / (a + 1e-10)
susie.a = a

susie.b = np.ones((susie.dims['T'], susie.dims['K']))
susie.weight_precision_a = susie.a

fit_args = {
    'max_iter': 2,
    'update_covariate_weights': True,
    'update_weights': True,
    'update_pi': True,
    'ARD_weights': True,
    'update_variance': True,
    'verbose': True,
}
susie.fit(**fit_args)

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
    susie_t = M(**data_t, K=5)

    fit_args = {
        'max_iter':100,
        'verbose':True,
        'ARD_weights':False
    }
    susie_t.fit(**fit_args)
    compute_records(susie_t)
    strip(susie_t)
    susies[t] = susie_t

print('saving susie')
pickle.dump(susies, open(snakemake.output.susie, 'wb'))

