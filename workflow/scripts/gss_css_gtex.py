import pandas as pd
import pickle
import numpy as np
import json

from coloc.independent_model_ss import IndependentFactorSER as GSS
from coloc.cafeh_ss import CAFEH as CSS

from coloc.misc import *
import seaborn as sns
import matplotlib.pyplot as plt

def fit_css(LD, B, S, K, fit_args, path):
    if os.path.isfile(path):
        print('loading model')
        css = pickle.load(open(path, 'rb'))
    else:
        print('fitting model')
        css = CSS(LD, B, S, K=20)
        css.prior_activity = np.ones(20) * 0.01
        css.fit(**fit_args, update_active=False)
        css.fit(**fit_args, update_active=True)
        compute_records_css(css)
        strip_and_dump(css, path)
    rehydrate_model(css)
    css.LD = LD
    css.B = B
    css.S = S
    return css

def fit_gss(X, Y, covariates, snps, tissues, samples, K, fit_args, path):
    if os.path.isfile(path):
        print('loading model')
        gss = pickle.load(open(path, 'rb'))
    else:
        print('fitting model')
        gss = GSS(
            X=X, Y=Y, covariates=covariates,
            snp_ids=snps, tissue_ids=tissues, sample_ids=samples,
            K=20)
        gss.prior_activity = np.ones(20) * 0.01
        gss.fit(**fit_args, update_active=False)
        gss.fit(**fit_args, update_active=True)
        compute_records_gss(gss)
        strip_and_dump(gss, path)
    
    rehydrate_model(gss)
    gss.X = X
    gss.Y = Y
    gss.covariates = covariates
    return gss

annotations = pd.read_csv(
    '/work-zfs/abattle4/lab_data/GTEx_v8/sample_annotations/'
    'GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt', sep='\t', index_col=0)

ep = snakemake.input.expression
gp = snakemake.input.genotype_gtex
gp1kG = snakemake.input.genotype_1kG
ap = snakemake.input.assocations

v2r = json.load(open(snakemake.input.snp2rsid, 'r'))

# Load GTEx and 1kG genotype
print('loading genotypes...')
genotype, ref = load_genotype(gp)
genotype1kG, ref1kG = load_genotype(gp1kG)
genotype.rename(columns=v2r, inplace=True)

# flip 1kG genotypes to agrese with GTEx encoding
flip_1kG = {}
for snp in ref.keys():
    if snp in v2r and v2r[snp] in ref1kG:
        rsid = v2r[snp]
        flip_1kG[rsid] = ref[snp] == ref1kG[rsid]

flip_1kG = pd.Series(flip_1kG)
flip_1kG = flip_1kG[~flip_1kG].index.values
genotype1kG.loc[:, flip_1kG] = genotype1kG.loc[:, flip_1kG].applymap(flip)

# load data
print('loading data...')
data = make_gtex_genotype_data_dict(ep, gp)

# load GTEx summary stats
print('loading associations...')
B, S, V, n = get_gtex_summary_stats(ap)
[x.rename(columns=v2r, inplace=True) for x in [B, S, V, n]];

# filter down to list of snps present in GTEx and 1kG
print('filtering down to common snps')
common_snps = np.intersect1d(
    np.intersect1d(genotype1kG.columns, genotype.columns), B.columns)


X = np.nan_to_num(
    (genotype - genotype.mean(0)).loc[:, common_snps].values)
L = X / np.sqrt(np.nansum(X**2, 0))

L1kG = np.nan_to_num(
    (genotype1kG - genotype1kG.mean(0)).loc[:, common_snps].values)
L1kG = L1kG / np.sqrt(np.nansum(L1kG**2, 0))

K = 20
# fit gss
gss_fit_args = {
    'max_iter': 100,
    'update_covariate_weights': True,
    'update_weights': True,
    'update_pi': True,
    'ARD_weights': True,
    'update_variance': True,
    'verbose': True
}
gss = fit_gss(
    X.T, data['Y'], data['covariates'],
    common_snps, data['tissue_ids'], data['sample_ids'],
    K, gss_fit_args, snakemake.output.gss)

css_fit_args = {
    'update_weights': True,
    'update_pi': True,
    'ARD_weights': True,
    'update_variance': False,
    'verbose': True,
    'max_iter': 20
}
# fit css with LD estimate from GTEx
css = fit_css(L.T @ L, B.values, S.values, K, css_fit_args, snakemake.output.css_gtex)

# fit css with LD estimate from 1kG
css = fit_css(L1kG.T @ L1kG, B.values, S.values, K, css_fit_args, snakemake.output.css_1kG)

# fit css with corrected LD estimate from 1kG
alpha = 0.9
LD = alpha * L1kG.T @ L1kG + (1-alpha) * np.corrcoef(B.T)
css = fit_css(LD, B.values, S.values, K, css_fit_args, snakemake.output.css_1kG_corrected)

