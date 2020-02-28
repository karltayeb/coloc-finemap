import pickle
import pandas as pd
import numpy as np
import os

# load data
ld = pd.read_csv(snakemake.input.ld, sep='\t', header=None).values
associations = pd.read_csv(snakemake.input.associations, sep='\t', index_col=0)
snplist = np.squeeze(pd.read_csv(snakemake.input.snps, header=None).values)

# filter out nans
associations = associations.loc[:, ~np.any(np.isnan(associations), axis=0)]
mask = ~np.any(np.isnan(ld), axis=1)
ld = ld[mask][:, mask]
snplist = snplist[mask]

# filter to intersection
intersect = np.intersect1d(snplist, associations.columns.values)
associations = associations.loc[:, intersect]
mask = np.isin(snplist, intersect)
ld = ld[mask][:, mask]
snplist = snplist[mask]

data = {
    'LD': ld,
    'zscores': associations.values,
    'tissue_ids': associations.index.values,
    'variant_ids': associations.columns.values
}

pickle.dump(data, open(snakemake.output[0], 'wb'))

"""
gene = snakemake.wildcards.gene
print('Training model for {}'.format(gene))

#################
# load genotype #
#################
print('loading genotype matrix...')
#genotype_path = "/work-zfs/abattle4/lab_data/GTEx_v8_trans_eqtl_data_processed_by_brian/processed_genotypes/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_chr22_dosage_MAF_05.txt"
#chr22_genotype = pd.read_csv(genotype_path, sep='\t', index_col=0, na_values='-')

cis_variants = pd.read_csv(snakemake.input[0], index_col=0)

#####################
# load associations #
#####################
print('loading associations...')
associations = pd.read_csv('/work-zfs/abattle4/karl/gp_fine_mapping/eQTL_sign/associations/chr22/{}.associations.txt'.format(gene), sep='\t', header=None)
column_names = ['tissue', 'gene_id', 'variant_id', 'tss_distance','ma_samples', 'ma_count', 'maf', 'pval_nominal', 'slope', 'slope_se']
associations.columns = column_names
tissues = np.unique(associations.tissue)
associations.loc[:, 'z_score'] = associations.slope / associations.slope_se

###################
# load covariates #
###################
print('loading covariates...')
covariates = {tissue: pd.read_csv(
    '/work-zfs/abattle4/lab_data/GTEx_v8/ciseQTL/GTEx_Analysis_v8_eQTL_covariates/{}.v8.covariates.txt'.format(tissue),
    index_col=0, sep='\t') for tissue in tissues
}
covariates = pd.concat(covariates, sort=True)
samples = covariates.columns

covariates = {tissue: pd.read_csv(
    '/work-zfs/abattle4/lab_data/GTEx_v8/ciseQTL/GTEx_Analysis_v8_eQTL_covariates/{}.v8.covariates.txt'.format(tissue),
    index_col=0, sep='\t') for tissue in tissues
}

###################
# load expression #
###################
print('loading expression...')
expression = pickle.load(open('/work-zfs/abattle4/karl/gp_fine_mapping/notebooks/chr22_expression', 'rb'))
expression = pd.concat(expression, sort=True)
expression.reset_index(inplace=True)
expression.rename(columns={'level_0': 'tissue'}, inplace=True)
expression.set_index(['gene_id', 'tissue'], inplace=True)

##########
# filter #
##########
print('filtering...')
variant_ids = np.intersect1d(
    cis_variants.index[~np.any(np.isnan(cis_variants.loc[:, samples].values), 1)],
    np.unique(associations.variant_id)
)
associations = associations[associations.variant_id.isin(variant_ids)]

# get zscores and LD matrix
z_scores = associations.pivot('tissue', 'variant_id', 'z_score')
nan_filter = np.all(~np.isnan(z_scores), 0)
p_filter = (associations.pivot('tissue', 'variant_id', 'pval_nominal').min(0) < 0.5)

z_scores = z_scores.iloc[:, nan_filter.values & p_filter.values]
cis_variants = cis_variants.loc[variant_ids].iloc[nan_filter.values & p_filter.values]
cis_variants = cis_variants.loc[:, samples]
##############
# compute LD #
##############
LD = np.corrcoef(cis_variants.values)
X = (cis_variants - cis_variants.mean(1)[:, None]).values / \
    np.sqrt(np.var(cis_variants, 1))[:, None]
data = {
    'X': X,
    'LD': LD,
    'Y': expression.loc[gene, samples].values,
    'zscores': z_scores.values,
    'covariates': covariates,
    'variant_ids': cis_variants.index.values,
    'sample_ids': samples,
    'tissue_ids': tissues
}
pickle.dump(data, open(snakemake.output[0], 'wb'))
"""