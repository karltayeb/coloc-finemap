import pandas as pd
import pickle
import numpy as np

from coloc.independent_model_ss import IndependentFactorSER as GSS
from coloc.cafeh_ss import CAFEH as CSS

from coloc.misc import *
import seaborn as sns
import matplotlib.pyplot as plt

ep = snakemake.input.expression
gp = snakemake.input.genotype
ap = snakemake.input.associations

def need_to_flip(variant_id):
    _, _, major, minor, _, ref = variant_id.strip().split('_')
    if minor != ref:
        return True
    else:
        return False

flip = lambda x: (x-1)*-1 + 1

def make_gtex_genotype_data_dict(ep, gp):
    gene_expression = pd.read_csv(ep, sep='\t', index_col=0)

    #load genotype
    genotype = pd.read_csv(gp, sep=' ')
    genotype = genotype.set_index('IID').iloc[:, 5:]

    # recode genotypes
    coded_snp_ids = np.array([x.strip() for x in genotype.columns])
    snp_ids = np.array(['_'.join(x.strip().split('_')[:-1]) for x in coded_snp_ids])

    flips = np.array([need_to_flip(vid) for vid in coded_snp_ids])
    genotype.iloc[:, flips] = genotype.iloc[:, flips].applymap(flip)

    # center, mean immpute
    genotype = (genotype - genotype.mean(0))
    genotype = genotype.fillna(0)

    # standardize
    genotype = genotype / genotype.std(0)

    # drop individuals that do not have recorded expression
    gene_expression = gene_expression.loc[:, ~np.all(np.isnan(gene_expression), 0)]

    # filter down to relevant individuals
    genotype = genotype.loc[gene_expression.columns]

    # filter down to common individuals
    individuals = np.intersect1d(genotype.index.values, gene_expression.columns.values)
    genotype = genotype.loc[individuals]
    gene_expression = gene_expression.loc[:, individuals]

    # load covariates
    covariates = pd.read_csv('/work-zfs/abattle4/karl/cosie_analysis/output/GTEx/covariates.csv',
                             sep='\t', index_col=[0, 1])
    covariates = covariates.loc[gene_expression.index]
    covariates = covariates.loc[:, genotype.index.values]

    X = genotype.values.T

    data = {
        'X': X,
        'Y': gene_expression.values,
        'covariates': covariates,
        'snp_ids': snp_ids,
        'sample_ids': genotype.index.values,
        'tissue_ids': gene_expression.index.values
    }
    return data

def compute_summary_stats(data):
    B = {}
    S = {}
    for i, tissue in enumerate(data['tissue_ids']):
        cov = data['covariates'].loc[tissue]
        mask = ~np.isnan(cov.iloc[0])
        cov = cov.values[:, mask]
        y = data['Y'][i, mask]
        X = data['X'][:, mask]

        #H = cov.T @ np.linalg.solve(cov @ cov.T, cov)
        H = (np.linalg.pinv(cov) @ cov)
        yc = y - y @ H
        Xc = X - X @ H

        # prep css data
        B[tissue] = (Xc @ yc) / np.einsum('ij,ij->i', Xc, Xc)
        r = yc - B[tissue][:, None]*Xc
        V = np.einsum('ij,ij->i', r, r) / np.einsum('ij,ij->i', Xc, Xc) / (yc.size)
        S[tissue] = np.sqrt(B[tissue]**2/yc.size + V)

    B = pd.DataFrame(B, index=snp_ids).T
    S = pd.DataFrame(S, index=snp_ids).T
    return B, S

def get_summary_stats_from_gtex(ap):
    associations = pd.read_csv(ap)
    associations.loc[:, 'sample_size'] = (associations.ma_count / associations.maf / 2)
    Ba = associations.pivot('tissue', 'variant_id', 'slope')
    Va = associations.pivot('tissue', 'variant_id', 'slope_se')**2
    n = associations.pivot('tissue', 'variant_id', 'sample_size')
    Sa = np.sqrt(Ba**2/n + Va)
    return Ba, Sa, n

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
Ba, Sa, n = get_summary_stats_from_gtex(ap)

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
strip_and_dump(gss, 'gss')

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
strip_and_dump(css, 'css')
strip_and_dump(css, 'css')