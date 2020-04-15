import pandas as pd
import numpy as np
import pickle

def make_gtex_genotype_data_dict(expression_path, genotype_path, params):
    # load expression
    gene_expression = pd.read_csv(expression_path, sep='\t', index_col=0)

    #load genotype
    genotype = pd.read_csv(genotype_path, sep=' ')
    genotype = genotype.set_index('IID').iloc[:, 5:]

    # center, mean immpute
    genotype = (genotype - genotype.mean(0))
    genotype = genotype.fillna(0)

    # standardize
    if params.standardize == True:
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
    covariates = pd.read_csv('/work-zfs/abattle4/karl/cosie_analysis/output/GTEx/covariates.csv', sep='\t', index_col=[0, 1])
    covariates = covariates.loc[gene_expression.index]
    covariates = covariates.loc[:, genotype.index.values]

    X = genotype.values.T

    data = {
        'X': X,
        'Y': gene_expression.values,
        'covariates': covariates,
        'snp_ids': genotype.columns.values,
        'sample_ids': genotype.index.values,
        'tissue_ids': gene_expression.index.values
    }

    return data


data = make_gtex_genotype_data_dict(
    snakemake.input.expression,
    snakemake.input.genotype,
    snakemake.params
)

pickle.dump(data, open(snakemake.output[0], 'wb'))
