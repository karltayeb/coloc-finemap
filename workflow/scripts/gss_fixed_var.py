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


def fit_gss(init_args, fit_args, path):
    print('fitting model')
    gss = GSS(**init_args)
    gss.tissue_precision_b = np.nanvar(gss.Y, 1)

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
flip_gtex = np.intersect1d(flip_gtex, genotype.columns)
genotype.loc[:, flip_gtex] = genotype.loc[:, flip_gtex].applymap(flip)
genotype.rename(columns=v2r, inplace=True)

# load data
print('loading data...')
data = make_gtex_genotype_data_dict(ep, gp)

# load GTEx summary stats
print('loading associations...')
B, S, V, n = get_gtex_summary_stats(ap)
[x.rename(columns=v2r, inplace=True) for x in [B, S, V, n]];

# filter down to list of snps present in GTEx and 1kG
print('filtering down to common snps')
common_snps = np.intersect1d(B.columns, table.rsid)
common_snps = np.intersect1d(common_snps, genotype1kG.columns)

X = np.nan_to_num(
    (genotype - genotype.mean(0)).loc[:, common_snps].values)

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
    'update_variance': False,
    'verbose': True
}
gss = fit_gss(gss_init_args, gss_fit_args, snakemake.output.gss)