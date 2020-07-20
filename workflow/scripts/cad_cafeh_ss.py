import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib import cm

from collections import defaultdict
import json
from collections import namedtuple
import ast
from coloc.misc import *
from coloc.cafeh_ss import CAFEH as CSS

import pysam
import copy

from collections import namedtuple
import ast


gc = pd.read_csv('../../output/GTEx/gencode.tss.bed', sep='\t', header=None)
gc.loc[:, 'g'] = gc.iloc[:, 3].apply(lambda x: x.split('.')[0])

sample_ld = lambda g: np.corrcoef(center_mean_impute(g), rowvar=False)

def cast(s):
    try:
        return ast.literal_eval(s)
    except Exception:
        return s


def load_cad_gwas(gene, variants=None):
    gwas = pysam.TabixFile('../../output/CAD/CAD_META.sorted.txt.gz')
    tss = gc[gc.iloc[:, 3]==gene].iloc[0][1]
    chrom = int(get_chromosome(gene)[3:])
    df = pd.DataFrame(
        list(map(autocast, x.strip().split('\t')) for x in
             gwas.fetch(chrom, np.clip(tss-1e6, 0, None), tss+1e6)),
        columns=gwas.header[0][1:].strip().split('\t')
    )
    if variants is not None:
        df = df[df.oldID.isin(variants)]
        
    df.rename(columns={
        'Allele1': 'ref',
        'Allele2': 'alt',
        'Effect': 'slope',
        'StdErr': 'slope_se',
        'P-val': 'pval_nominal',
        'oldID': 'rsid'}, inplace=True)

    df['ref'] = df['ref'].str.upper()
    df['alt'] = df['alt'].str.upper()
    df.loc[:, 'S'] = df['slope_se']  # large sample approximation
    df.loc[:, 'study'] = 'CAD'

    df=df.rename(columns={
        'study': 'tissue',
        'P-value': 'pval_nominal'
    })
    return df


def load_gtex_associations(gene):
    ap = '/work-zfs/abattle4/karl/cosie_analysis/output/GTEx/{}/{}/{}.associations'.format(
        get_chromosome(gene), gene, gene)
    v2rp = '/work-zfs/abattle4/karl/cosie_analysis/output/GTEx/{}/{}/{}.snp2rsid.json'.format(
        get_chromosome(gene), gene, gene)
    v2r = json.load(open(v2rp, 'r'))

    associations = pd.read_csv(ap, index_col=0)
    associations.loc[:, 'rsid'] = associations.variant_id.apply(lambda x: v2r.get(x, '-'))
    associations.loc[:, 'pos'] = associations.variant_id.apply(lambda x: int(x.split('_')[1]))
    associations.loc[:, 'ref'] = associations.variant_id.apply(lambda x: x.split('_')[2])
    associations.loc[:, 'alt'] = associations.variant_id.apply(lambda x: x.split('_')[3])
    associations.loc[:, 'sample_size'] = (associations.ma_count / associations.maf / 2)
    associations.loc[:, 'S'] = np.sqrt(
        associations.slope**2/associations.sample_size + associations.slope_se**2)
    return associations


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


def association2summary(associations):
    return SimpleNamespace(**{
        'B': associations.pivot('tissue', 'rsid', 'slope'),
        'S': associations.pivot('tissue', 'rsid', 'S')
    })


def init_css(summary_stats, genotype, K=10, pi0=1.0, dispersion=1.0, name='css', **kwargs):
    # TODO epsilon--- smooth expression?
    init_args = {
        'LD': sample_ld(genotype),
        'B': summary_stats.B.fillna(0).values,
        'S': summary_stats.S.fillna(100).values,
        'K': K,
        'snp_ids': summary_stats.B.columns.values,
        'tissue_ids': summary_stats.B.index.values,
        'tolerance': 1e-8
    }
    fit_args = {
        'update_weights': True,
        'update_pi': True,
        'ARD_weights': True,
        'update_variance': False,
        'verbose': True,
        'max_iter': 50
    }
    print('initializing summary stat model')
    css = CSS(**init_args)
    css.prior_activity = np.ones(K) * pi0
    css.tissue_precision_b = np.ones(css.dims['T']) * dispersion
    css.name = name
    return css, fit_args


def fit_css(summary_stats, genotype, save_path=None, **kwargs):
    css, fit_args = init_css(summary_stats, genotype, **kwargs)
    print('fitting model')
    css.fit(**fit_args, update_active=False)
    css.fit(**fit_args, update_active=True)
    compute_records_css(css)
    if save_path is not None:
        print('saving model to {}'.format(save_path))
        strip_and_dump(css, save_path)
        rehydrate_model(css)
    return css

# load gwas and associations
gene = snakemake.wildcards.gene
save_path = snakemake.output[0]

gwas = load_cad_gwas(gene)
associations = load_gtex_associations(gene)
gtex_genotype = load_gtex_genotype(gene)

common_columns = np.intersect1d(associations.columns, gwas.columns)
all_associations= pd.concat([associations.loc[:, common_columns], gwas.loc[:, common_columns]])

rsid2pos = associations.set_index('rsid').loc[:, 'pos'].to_dict()
all_associations.loc[:, 'pos'] = all_associations.rsid.apply(lambda x: rsid2pos.get(x, np.nan))
common_variants = np.intersect1d(
    associations.rsid.unique(), gwas.rsid.unique())
summary_stats = association2summary(all_associations[all_associations.rsid.isin(common_variants)])

css = fit_css(
    summary_stats,
    gtex_genotype.loc[:, common_variants], save_path=save_path, pi0=0.01)
