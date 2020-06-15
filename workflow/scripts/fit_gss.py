import pandas as pd
import pickle
import numpy as np
from coloc.independent_model_ss import IndependentFactorSER as GSS
from coloc.misc import *

def init_gss(data, K=10, p=1.0):
    print('initializing genotype model')
    gss = GSS(
        X=center_mean_impute(
            data.genotype_gtex.loc[data.expression.columns]).values.T,
        Y=data.expression.values,
        K=K,
        covariates=data.covariates.loc[:, data.expression.columns],
        snp_ids=data.common_snps,
        tissue_ids=data.expression.index.values,
        sample_ids=data.expression.columns.values
    )
    
    gss.prior_activity = np.ones(K) * p
    return gss

#data = pickle.load(open(snakemake.input[0], 'rb'))
data = load_gene_data(snakemake.wildcards.gene)

K = snakemake.params.k
p = snakemake.params.pi0

gss = init_gss(data, K, p)

fit_args = {
    'max_iter': 300,
    'update_covariate_weights': True,
    'update_weights': True,
    'update_pi': True,
    'ARD_weights': True,
    'update_variance': True,
    'verbose': True
}
print('training model')
for arg in fit_args:
    print('\t{}: {}'.format(arg, fit_args[arg]))
gss.fit(**fit_args, update_active=False)
gss.fit(**fit_args, update_active=True)


print('model fit:\n\titers:{}\n\tELBO:{}\n\trun-time:{}'.format(
    len(gss.elbos), gss.elbos[-1], gss.run_time))

print('saving model')
compute_records(gss)
strip_and_dump(gss, snakemake.output[0])
