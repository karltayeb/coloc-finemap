import pandas as pd
import pickle
import numpy as np
import json
import os
import sys

from coloc.independent_model_ss import IndependentFactorSER as GSS
from coloc.cafeh_ss import CAFEH as CSS

from coloc.misc import *
import seaborn as sns
import matplotlib.pyplot as plt


def fit_css(init_args, fit_args, path):
    print('fitting model')
    css = CSS(**init_args)
    css.prior_activity = np.ones(20) * 0.01
    css.fit(**fit_args, update_active=False)
    css.fit(**fit_args, update_active=True)
    compute_records_css(css)
    strip_and_dump(css, path)
    rehydrate_model(css)
    return css

def fit_gss(init_args, fit_args, path):
    print('fitting model')
    gss = GSS(**init_args)
    gss.prior_activity = np.ones(20) * 0.01
    gss.fit(**fit_args, update_active=False)
    gss.fit(**fit_args, update_active=True)
    compute_records_gss(gss)
    strip_and_dump(gss, path)
    rehydrate_model(gss)
    return gss

def make_snp_format_table(gp, gp1kG, v2rp):
    with open(gp, 'r') as f:
        snps = f.readline().strip().split()[6:]

    with open(gp1kG, 'r') as f:
        rsids = f.readline().strip().split()[6:]

    v2r = json.load(open(v2rp, 'r'))

    vid_codes = {'_'.join(x.split('_')[:-1]): x.split('_')[-1] for x in snps}
    rsid_codes = {x.split('_')[0]: x.split('_')[1] for x in rsids}
    table = []
    for vid in vid_codes:
        ref = vid.split('_')[-2]
        rsid = v2r.get(vid, '-')
        table.append({
            'variant_id': vid,
            'rsid': v2r.get(vid, '-'),
            'ref': ref,
            'flip_gtex': ref != vid_codes.get(vid, '-'),
            'flip_1kG': ref != rsid_codes.get(rsid, '-')
        })
    return pd.DataFrame(table)

annotations = pd.read_csv(
    '/work-zfs/abattle4/lab_data/GTEx_v8/sample_annotations/'
    'GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt', sep='\t', index_col=0)

ep = snakemake.input.expression
gp = snakemake.input.genotype_gtex
gp1kG = snakemake.input.genotype_1kG
ap = snakemake.input.associations
v2rp = snakemake.input.snp2rsid
v2r = json.load(open(v2rp, 'r'))

table = make_snp_format_table(gp, gp1kG, v2rp)

# Load GTEx and 1kG genotype
# flip genotype encoding to be consistent with GTEx associations
print('loading genotypes...')
genotype, ref = load_genotype(gp)
flip_gtex = table[table.flip_gtex].variant_id.values
genotype.loc[:, flip_gtex] = genotype.loc[:, flip_gtex].applymap(flip)
genotype.rename(columns=v2r, inplace=True)

genotype1kG, ref1kG = load_genotype(gp1kG)
flip_1kG = table[table.flip_1kG & (table.rsid != '-')].rsid.values
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

B = B.loc[data['tissue_ids'], common_snps]
S = S.loc[data['tissue_ids'], common_snps]
V = V.loc[data['tissue_ids'], common_snps]

# replace missnig values with uninformative values
mask = ~np.isnan(S.values)
B.fillna(0, inplace=True)
S.fillna(1e2, inplace=True)
V.fillna(1e2, inplace=True)

X = np.nan_to_num(
    (genotype - genotype.mean(0)).loc[:, common_snps].values)
L = X / np.sqrt(np.nansum(X**2, 0))

L1kG = np.nan_to_num(
    (genotype1kG - genotype1kG.mean(0)).loc[:, common_snps].values)
L1kG = L1kG / np.sqrt(np.nansum(L1kG**2, 0))

K = 10

gss_init_args = {
    'X': X.T,
    'Y': data['Y'],
    'covariates': data['covariates'],
    'snp_ids': common_snps,
    'tissue_ids': data['tissue_ids'],
    'sample_ids': data['sample_ids'],
    'K': 20
}
# fit gss
gss_fit_args = {
    'max_iter': 50,
    'update_covariate_weights': True,
    'update_weights': True,
    'update_pi': True,
    'ARD_weights': True,
    'update_variance': True,
    'verbose': True
}
gss = fit_gss(gss_init_args, gss_fit_args, snakemake.output.gss)


css_fit_args = {
    'update_weights': True,
    'update_pi': True,
    'ARD_weights': True,
    'update_variance': False,
    'verbose': True,
    'max_iter': 20
}
# fit css with LD estimate from GTEx
css_init_args = {
    'LD': L.T @ L,
    'B': B.values,
    'S': S.values,
    'K': K,
    'snp_ids': common_snps,
    'tissue_ids': B.index.values
}
css = fit_css(css_init_args, css_fit_args, snakemake.output.css_gtex)

css_init_args = {
    'LD': L1kG.T @ L1kG,
    'B': B.values,
    'S': S.values,
    'K': K,
    'snp_ids': common_snps,
    'tissue_ids': B.index.values
}
# fit css with LD estimate from 1kG
css = fit_css(css_init_args, css_fit_args, snakemake.output.css_1kG)

# fit css with corrected LD estimate from 1kG
alpha = 0.9
LD = alpha * L1kG.T @ L1kG + (1-alpha) * np.corrcoef(B.T)
css_init_args = {
    'LD': LD,
    'B': B.values,
    'S': S.values,
    'K': K,
    'snp_ids': common_snps,
    'tissue_ids': B.index.values
}
css = fit_css(css_init_args, css_fit_args, snakemake.output.css_1kG_corrected)
