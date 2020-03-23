import pandas as pd
import numpy as np
import pickle

gene_expression = pd.read_csv(snakemake.input.expression, sep='\t', index_col=0)

genotype = pd.read_csv(snakemake.input.genotype, sep=' ')
genotype = genotype.set_index('IID').iloc[:, 5:]
genotype = (genotype - genotype.mean(0)) / genotype.std(0)

# drop individuals that do not have recorded expression
gene_expression = gene_expression.loc[:, ~np.all(np.isnan(gene_expression), 0)]

# filter down to relevant individuals
genotype = genotype.loc[gene_expression.columns]

# filter out snps with nans
genotype = genotype.loc[:, ~np.any(np.isnan(genotype), 0)]

covariates = {}
for tissue in gene_expression.index.values:
    covariates[tissue] = pd.read_csv(
        '/work-zfs/abattle4/lab_data/GTEx_v8/'
        'ciseQTL/GTEx_Analysis_v8_eQTL_covariates/{}.v8.covariates.txt'.format(tissue),
        sep='\t', index_col=0
    )

X = genotype.values.T
X = (X - X.mean(1)[:, None]) # / np.clip(X.std(1)[:, None], 1e-10, 1e10)

data = {
    'X': X,
    'Y': gene_expression.values,
    'covariates': covariates,
    'snp_ids': genotype.columns.values,
    'sample_ids': genotype.index.values,
    'tissue_ids': gene_expression.index.values
}

pickle.dump(data, open(snakemake.output[0], 'wb'))