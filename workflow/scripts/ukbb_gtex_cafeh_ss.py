import numpy as np
import pandas as pd
import scipy as sp

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


def filter_and_flip(gtex, gwas, variants):
    """
    filter down to common variants with unique coding
    flip gwas to match gtex
    """
    common_variants = np.intersect1d(gwas.rsid, gtex.rsid)
    common_variants = np.intersect1d(common_variants, variants)

    a = gtex.loc[gtex.rsid.isin(common_variants), ['rsid', 'ref', 'alt']]
    a = a.drop_duplicates().drop_duplicates('rsid', keep=False).set_index('rsid')

    b = gwas.loc[gwas.rsid.isin(common_variants), ['rsid', 'ref', 'alt']]
    b = b.drop_duplicates().drop_duplicates('rsid', keep=False).set_index('rsid')

    common_variants = np.intersect1d(b.index, a.index)

    a = a.loc[common_variants]
    b = b.loc[common_variants]
    flip = common_variants[(a.ref != b.ref)]
    
    gtex = gtex[gtex.rsid.isin(common_variants)]
    gwas = gwas[gwas.rsid.isin(common_variants)]
    
    gwas.loc[gwas.rsid.isin(flip), 'slope'] = -1 * gwas.loc[gwas.rsid.isin(flip), 'slope']
    
    ref = gwas.loc[gwas.rsid.isin(flip), 'ref']
    gwas.loc[gwas.rsid.isin(flip), 'ref'] = gwas.loc[gwas.rsid.isin(flip), 'alt']
    gwas.loc[gwas.rsid.isin(flip), 'alt'] = ref
    return gtex, gwas, flip


gene = snakemake.wildcards.gene
study = snakemake.wildcards.study
phenotype = snakemake.wildcards.phenotype

# load gtex and gtex genotype
gtex_genotype = load_gtex_genotype(gene, use_rsid=True)
gtex = load_gtex_associations(gene)

print(study)
# load gwas
if 'UKBB' in study:
    gwas = load_ukbb_gwas(gene, phenotype)
elif 'GRASP' in study:
    gwas = load_grasp_gwas(gene, phenotype)
else:
    gwas = load_cad_gwas(gene)
# flip signs in gwas, filter down to variants with genotype, eqtl, and gwas
gtex, gwas, flip = filter_and_flip(gtex, gwas, gtex_genotype.columns)

common_columns = np.intersect1d(gtex.columns, gwas.columns)
all_associations = pd.concat([gtex.loc[:, common_columns], gwas.loc[:, common_columns]])

rsid2pos = associations.set_index('rsid').loc[:, 'pos'].to_dict()
all_associations.loc[:, 'pos'] = all_associations.rsid.apply(lambda x: rsid2pos.get(x, np.nan))

all_associations.loc[:, 'z'] = all_associations.slope / all_associations.slope_se
all_associations.loc[:, 'zS'] = np.sqrt((all_associations.z**2 / all_associations.sample_size) + 1)

z = all_associations[~all_associations.duplicated(['tissue', 'rsid'])].pivot('tissue', 'rsid', 'z')
zS = all_associations[~all_associations.duplicated(['tissue', 'rsid'])].pivot('tissue', 'rsid', 'zS')

variants = all_associations.rsid.unique()
fully_observed_idx = (~np.any(z.isna(), 0)).values
fully_observed_variants = z.columns[fully_observed_idx].values

print('{} variants in gwas'.format(variants.size))
print('{} variant fully observed in GTEx'.format(fully_observed_variants.size))

#############################
#  impute zscores for gtex  #
#############################

if snakemake.params.impute:
    # A = (X.T @ X + eI)^-1
    print('imputing missing genotype')
    X = center_mean_impute(gtex_genotype.loc[:, variants]).values
    X = X / (X.std(0) * np.sqrt(X.shape[0]))
    LD = X.T @ X

    idx = np.isin(variants, fully_observed_variants)
    X_o = X[:, idx]

    m, n = X_o.shape
    e = 0.01
    Finv = np.eye(m) + 1/e * X_o @ X_o.T
    R = sp.linalg.solve_triangular(np.linalg.cholesky(Finv).T, np.eye(m))
    V = X_o.T @ R / e
    A = (np.eye(n)/e - V @ V.T) @ X_o.T @ X[:, ~idx]

    # impute z-scores within cis-window
    _z_imp = pd.DataFrame(
        z.loc[:, fully_observed_idx].values @ A,
        index=z.index, columns=variants[~idx])

    z_imp = z.copy()
    z_imp.loc[:, _z_imp.columns] = _z_imp

    # data for model
    studies = z_imp.index.values
    B = z_imp.loc[:, variants].values
    S = zS.loc[:, variants].fillna(1.0).values

else:
    variants = fully_observed_variants
    X = center_mean_impute(gtex_genotype.loc[:, variants]).values
    X = X / (X.std(0) * np.sqrt(X.shape[0]))
    LD = X.T @ X

    variants = fully_observed_variants
    studies = z.index.values
    B = z.loc[:, variants].values
    S = zS.loc[:, variants].values

K = snakemake.params.K

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
