import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib import cm

from collections import defaultdict
import json
from collections import namedtuple
import ast
import sys
sys.path.append('/work-zfs/abattle4/karl/cosie_analysis/utils/')
from misc import *

from cafeh.cafeh_ss import CAFEH as CSS
from cafeh.fitting import weight_ard_active_fit_procedure, fit_all


import pysam
import copy

from collections import namedtuple
import ast

sample_ld = lambda g: np.corrcoef(center_mean_impute(g), rowvar=False)

def cast(s):
    try:
        return ast.literal_eval(s)
    except Exception:
        return s

def _load_ukbb_chromosome(phenotype, chromosome):
    ukbb_columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'ac', 'af', 'num_cases',
       'num_controls', 'beta', 'sebeta', 'Tstat', 'pval',
       'pval_SAIGE_NoSPA', 'Is_Converged', 'varT', 'varTstar']
    path = ('/work-zfs/abattle4/marios/HF_GWAS/'
            'summary_stats_for_genetic_correlation/Phecode4_{}.sumstats.txt'.format(phenotype))
    chr2startidx = json.load(open('output/UKBB/{}/chr2startidx'.format(phenotype)))
    start = chr2startidx.get(str(chromosome))
    nrows = chr2startidx.get(str(int(chromosome) + 1)) - start + 1
    df = pd.read_csv(path, sep='\t', skiprows=start, nrows=nrows, header=None)
    df.columns = ukbb_columns
    return df

def load_ukbb_gwas(gene, phenotype, variants=None):
    tss = get_tss(gene)
    chrom = get_chr(gene)[3:]

    df = _load_ukbb_chromosome(phenotype, chrom)
    df = df[(df.POS > tss-1e6) & (df.POS < tss+1e6) & (df.iloc[:, 0] == int(chrom))]
    df.rename(columns={
        'REF': 'ref',
        'ALT': 'alt',
        'beta': 'slope',
        'sebeta': 'slope_se',
        'pval': 'pval_nominal',
        'ID': 'rsid'}, inplace=True)

    df['ref'] = df['ref'].str.upper()
    df['alt'] = df['alt'].str.upper()
    df.loc[:, 'S'] = df['slope_se']  # large sample approximation

    df.loc[:, 'study'] = phenotype

    df=df.rename(columns={
        'study': 'tissue',
        'P-value': 'pval_nominal'
    })

    df.loc[:, 'sample_size'] = df.num_cases + df.num_controls
    return df


def load_grasp_gwas(gene, phenotype):
    df = pd.read_csv('output/GRASP/{}.txt'.format(phenotype), sep='\t')
    df.columns = ['chr', 'pos', 'variant_id', 'rsid', 'slope', 'slope_se', 'pval_nominal', 'sample_size']

    tss = get_tss(gene)
    chrom = int(get_chr(gene)[3:])

    df.loc[:, 'pos'] = df.variant_id.apply(lambda x: int(x.split('_')[1]))
    df.loc[:, 'ref'] = df.variant_id.apply(lambda x: x.split('_')[2])
    df.loc[:, 'alt'] = df.variant_id.apply(lambda x: x.split('_')[3])
    df = df[(df.chr==chrom) & (df.pos > tss-1e6) & (df.pos < tss+1e6)]

    df['ref'] = df['ref'].str.upper()
    df['alt'] = df['alt'].str.upper()
    df.loc[:, 'S'] = df['slope_se']  # large sample approximation

    df.loc[:, 'study'] = phenotype

    df=df.rename(columns={
        'study': 'tissue',
        'P-value': 'pval_nominal'
    })
    return df


def flip(associations, gwas):
    """
    flip gwas asoociations to match GTEx
    """
    A = associations.drop_duplicates('rsid').set_index('rsid').loc[:, ['alt', 'ref']]

    B = gwas.set_index('rsid').loc[:, ['alt', 'ref']]
    B = B.loc[B.index.drop_duplicates()]
    B = B.applymap(lambda x: np.str.upper(str(x)))

    agree = pd.DataFrame(A).merge(B, left_index=True, right_index=True)
    agree.loc[:, 'flip'] = (agree.alt_x != agree.alt_y)

    gwas = gwas.set_index('rsid')
    gwas.loc[agree[agree.flip == True].index, 'slope'] \
        = gwas.loc[agree[agree.flip == True].index, 'slope']*-1
    gwas = gwas.reset_index()
    return gwas

phenotype = snakemake.wildcards.phenotype

GENE = snakemake.wildcards.gene
ASSOCIATION_PATH = snakemake.input.associations


# load gwas and associations
associations = load_gtex_associations(GENE)

if snakemake.wildcards.study == 'UKBB':
    print('loading UKBB gwas')
    gwas = flip(associations, load_ukbb_gwas(GENE, phenotype))
if snakemake.wildcards.study == 'GRASP':
    print('loading GRASP gwas')
    # gwas = flip(associations, load_grasp_gwas(GENE, phenotype))
    gwas = load_grasp_gwas(phenotype)

gtex_genotype = load_gtex_genotype(GENE, use_rsid=True)
common_columns = np.intersect1d(associations.columns, gwas.columns)
all_associations = pd.concat([associations.loc[:, common_columns], gwas.loc[:, common_columns]])

rsid2pos = associations.set_index('rsid').loc[:, 'pos'].to_dict()
all_associations.loc[:, 'pos'] = all_associations.rsid.apply(lambda x: rsid2pos.get(x, np.nan))

# filter down to common variants
common_variants = np.intersect1d(
    associations.rsid.unique(), gwas.rsid.unique())
common_variants = np.intersect1d(
    common_variants, gtex_genotype.columns)
gtex_genotype = gtex_genotype.loc[:, common_variants]
gtex_genotype = gtex_genotype.loc[:, ~gtex_genotype.columns.duplicated()]
all_associations = all_associations[all_associations.rsid.isin(common_variants)]
all_associations = all_associations.drop_duplicates(['rsid', 'tissue'])
all_associations.loc[:, 'z'] = all_associations.slope / all_associations.slope_se
all_associations.loc[:, 'zS'] = np.sqrt((all_associations.z**2 / all_associations.sample_size) + 1)

# make summary stat tables
z = all_associations[~associations.duplicated(['tissue', 'rsid'])].pivot('tissue', 'rsid', 'z')
zS = all_associations[~associations.duplicated(['tissue', 'rsid'])].pivot('tissue', 'rsid', 'zS')

gwas_variants = z.columns[~z.loc[phenotype].isna()].values
fully_observed_idx = (~np.any(z.isna(), 0)).values
fully_observed_variants = z.columns[fully_observed_idx].values

print('{} variants in gwas'.format(gwas_variants.size))
print('{} variant fully observed in GTEx'.format(fully_observed_variants.size))

#############################
# fit fully observed model  #
#############################
"""
print('{} variants in gwas'.format(gwas_variants.size))
print('{} variant fully observed in GTEx'.format(fully_observed_variants.size))

variants = fully_observed_variants
studies = z.index.values
B = z.loc[:, variants].values
S = zS.loc[:, variants].values
K = 20

init_args = {
    'LD': sample_ld(gtex_genotype.loc[:, variants]),
    'B': B,
    'S': S,
    'K': K,
    'snp_ids': variants,
    'study_ids': studies,
    'tolerance': 1e-8
}
css = CSS(**init_args)
css.prior_activity = np.ones(K) * 0.1
print('fit model with fully observed variants')
weight_ard_active_fit_procedure(css, verbose=False, max_iter=50)

print('saving model to {}'.format(snakemake.output[0]))
css.save(snakemake.output[0])
"""
#############################
#  impute zscores for gtex  #
#############################

# compute LD

# A = (X.T @ X + eI)^-1
X = center_mean_impute(gtex_genotype.loc[:, fully_observed_variants]).values
X = X / (X.std(0) * np.sqrt(X.shape[0]))
m, n = X.shape
e = 0.1
Finv = np.eye(m) + 1/e * X @ X.T
R = sp.linalg.solve_triangular(np.linalg.cholesky(Finv).T, np.eye(m))
V = X.T @ R / e
A = np.eye(n)/e - V @ V.T

LD = X.T @ X
LD_oo = LD[fully_observed_idx][:, fully_observed_idx]
LD_uo = LD[~fully_observed_idx][:, fully_observed_idx]

"""
LD = np.corrcoef(np.nan_to_num(
    gtex_genotype.loc[:, common_variants].values
    - gtex_genotype.loc[:, common_variants].mean(0).values[None]), rowvar=False)
LD_oo = LD[fully_observed_idx][:, fully_observed_idx]
LD_uo = LD[~fully_observed_idx][:, fully_observed_idx]

A = np.linalg.solve(LD_oo + np.eye(fully_observed_idx.sum()) * 0.01, LD_uo.T)
"""

# impute z-scores within cis-window
_z_imp = pd.DataFrame(
    z.loc[:, fully_observed_idx].values @ A,
    index=z.index, columns=z.columns[~fully_observed_idx])
z_imp = z.copy()
z_imp.loc[:, _z_imp.columns] = _z_imp

#######################
#  fit imputed model  #
#######################
variants = gwas_variants
studies = z_imp.index.values
B = z_imp.loc[:, variants].values
S = zS.loc[:, variants].fillna(1.0).values
K = 20

init_args = {
    'LD': sample_ld(gtex_genotype.loc[:, variants]),
    'B': B,
    'S': S,
    'K': K,
    'snp_ids': variants,
    'study_ids': studies,
    'tolerance': 1e-8
}
css = CSS(**init_args)
css.prior_activity = np.ones(K) * 0.1
css.weight_precision_b = np.ones_like(css.weight_precision_b) * 10

print('fit model with imputed z-score')
weight_ard_active_fit_procedure(css, max_iter=10, verbose=True)
fit_all(css, max_iter=50, verbose=True)

print('saving model to {}'.format(snakemake.output[0]))
css.save(snakemake.output[0])
