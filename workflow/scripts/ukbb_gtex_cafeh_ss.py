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
from cafeh.fitting import weight_ard_active_fit_procedure


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


def load_ukbb_gwas(gene, phenotype, variants=None):
    gwas = pysam.TabixFile('output/UKBB/sumstats/Phecode4_{}.sumstats.txt.gz'.format(phenotype))
    tss = get_tss(gene)
    chrom = int(get_chr(gene)[3:])
    import pdb; pdb.set_trace()

    df = pd.DataFrame(
        list(map(cast, x.strip().split('\t')) for x in
             gwas.fetch(chrom, np.clip(tss-1e6, 0, None), tss+1e6)),
        columns=gwas.header[0][1:].strip().split('\t')
    )
    if variants is not None:
        df = df[df.oldID.isin(variants)]

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
gwas = flip(associations, load_ukbb_gwas(GENE, phenotype))

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


#############################
# fit fully observed model  #
#############################
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

#######################
#  fit imputed model  #
#######################

# compute LD
LD = np.corrcoef(np.nan_to_num(
    gtex_genotype.loc[:, common_variants].values
    - gtex_genotype.loc[:, common_variants].mean(0).values[None]), rowvar=False)
LD_oo = LD[fully_observed_idx][:, fully_observed_idx]
LD_uo = LD[~fully_observed_idx][:, fully_observed_idx]

A = np.linalg.solve(LD_oo + np.eye(fully_observed_idx.sum()) * 0.01, LD_uo.T)

# impute z-scores within cis-window
_z_imp = pd.DataFrame(
    z.loc[:, fully_observed_idx].values @ A,
    index=z.index, columns=z.columns[~fully_observed_idx])
z_imp = z.copy()
z_imp.loc[:, _z_imp.columns] = _z_imp

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
print('fit model with imputed z-score')
weight_ard_active_fit_procedure(css, verbose=False, max_iter=50)
print('saving model to {}'.format(snakemake.output[1]))
css.save(snakemake.output[1])
