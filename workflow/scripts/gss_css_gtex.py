import pandas as pd
import pickle
import numpy as np

from coloc.independent_model_ss import IndependentFactorSER as GSS
from coloc.cafeh_ss import CAFEH as CSS

from coloc.misc import *
import seaborn as sns
import matplotlib.pyplot as plt

ep = snakemake.input.expression
ap = snakemake.input.associations
gp = snakemake.input.genotype

def compute_records_gss(model):
    """
    save the model with data a weight parameters removed
    add 'mini_weight_measn' and 'mini_weight_vars' to model
    the model can be reconstituted in a small number of iterations
    """
    credible_sets, purity = model.get_credible_sets(0.999)
    active = model.active.max(0) > 0.5
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
        'mini_wm': wm,
        'mini_wv': wv,
        'snp_subset': mask
    }
    model.records = records


def compute_records_css(model):
    """
    save the model with data a weight parameters removed
    add 'mini_weight_measn' and 'mini_weight_vars' to model
    the model can be reconstituted in a small number of iterations
    """
    credible_sets, purity = model.get_credible_sets(0.999)
    active = model.active.max(0) > 0.5
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
    # purge precompute
    for key in model.precompute:
        model.precompute[key] = {}
    model.__dict__.pop('weight_means', None)
    model.__dict__.pop('weight_vars', None)
    if not save_data:
        model.__dict__.pop('X', None)
        model.__dict__.pop('Y', None)
        model.__dict__.pop('covariates', None)
        model.__dict__.pop('B', None)
        model.__dict__.pop('S', None)
        model.__dict__.pop('LD', None)
    pickle.dump(model, open(path, 'wb'))

print('loading data...')
data = make_gtex_genotype_data_dict(ep, gp)

print('computing summary stats...')
B, S = compute_summary_stats(data)

print('fetching GTEx summary stats...')
Ba, Sa, n = get_gtex_summary_stats(ap)

common_snps = np.intersect1d(Ba.columns, B.columns)
X = pd.DataFrame(data['X'], index=data['snp_ids']).loc[common_snps].values

Ba = Ba.loc[B.index, common_snps]
Sa = Sa.loc[B.index, common_snps]
B = B.loc[:, common_snps]
S = S.loc[:, common_snps]

# set all nan tests to high noise
mask = ~np.isnan(Sa.values)
Ba.fillna(0, inplace=True)
Sa.fillna(1e10, inplace=True)


# fit gss
gss = GSS(X=X, Y=data['Y'], covariates=data['covariates'], snp_ids=common_snps,
          tissue_ids=data['tissue_ids'], sample_ids=data['sample_ids'], K=20)
gss.prior_activity = np.ones(20) * 0.01
fit_args = {
    'max_iter': 100,
    'update_covariate_weights': True,
    'update_weights': True,
    'update_pi': True,
    'ARD_weights': True,
    'update_variance': True,
    'verbose': True
}

print('fitting gss')
gss.fit(**fit_args, update_active=False)
gss.fit(**fit_args, update_active=True)
compute_records_gss(gss)
strip_and_dump(gss, snakemake.output.gss)

#fit css
LD = np.corrcoef(X)

fit_args = {
    'update_weights': True,
    'update_pi': True,
    'ARD_weights': True,
    'update_variance': False,
    'verbose': True
}
print('fitting css with computed summary stats')
css = CSS(LD, B.values, S.values, K=20)
css.prior_activity = np.ones(20) * 0.01
css.fit(**fit_args, update_active=False, max_iter=50)
css.fit(**fit_args, update_active=True, max_iter=50)
compute_records_css(css)
strip_and_dump(css, snakemake.output.css)

#fit css.gtex
fit_args = {
    'update_weights': True,
    'update_pi': True,
    'ARD_weights': True,
    'update_variance': False,
    'verbose': True
}
print('fitting css with GTEx summary stats')
css = CSS(LD, Ba.values, Sa.values, K=20)
css.prior_activity = np.ones(20) * 0.01
css.fit(**fit_args, update_active=False, max_iter=50)
css.fit(**fit_args, update_active=True, max_iter=50)
compute_records_css(css)
strip_and_dump(css, snakemake.output.gtex_css)