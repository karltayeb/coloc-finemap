import numpy as np
import pandas as pd
from coloc.independent_model2 import IndependentFactorSER as M
from coloc.independent_model_summary import IndependentFactorSER as M2
from coloc.misc import *

genotype_path = snakemake.input.genotype
gene = snakemake.wildcards.gene

gencode = pd.read_csv(
    'output/GTEx/protein_coding_autosomal_egenes.txt', sep='\t')
gene = gencode.loc[gencode.gene == gene]
chr_num = gene.iloc[0, 0]
tss = gene.iloc[0, 1]
gene_name = gene.iloc[0, 3]

#load genotype
genotype = pd.read_csv(genotype_path, sep=' ')
genotype = genotype.set_index('IID').iloc[:, 5:]
# center, mean immpute
genotype = (genotype - genotype.mean(0))
genotype = genotype.fillna(0)
# standardize
genotype = genotype / genotype.std(0, ddof=0)
X = genotype.values.T

def make_simulation(X, T, pve, sparsity):
    # select snps
    N, M = X.shape
    causal_snps = np.random.choice(N, 5)
    true_effects = np.zeros((T, N))
    true_effects[:, causal_snps] = \
        np.random.binomial(1, sparsity, ((T, 5))) \

    tissue_variance = np.array([
        compute_sigma2(X, te, pve) for te in true_effects
    ])

    #simulate expression
    expression = (true_effects @ X) + \
        np.random.normal(size=(T, M)) * np.sqrt(tissue_variance)[:, None]
    expression = expression - expression.mean(1)[:, None]
    pos = np.array([int(x.split('_')[1]) for x in genotype.columns.values])

    # simulate expression 
    data = {
        'X': X,
        'Y': expression,
        'snp_ids': genotype.columns.values,
    }

    sim_info = {
        'causal_snps': causal_snps[true_effects[:, causal_snps].sum(0) > 0],
        'true_effects': true_effects,
        'tissue_variance': tissue_variance,
        'expression': expression
    }
    return data, sim_info

def compute_sigma2(X, true_effect, pve):
    var = np.var(true_effect @ X)
    sigma2_t = var/pve - var
    if sigma2_t == 0:
        # if variance is 0, there were no causal variants-- dont care what the variance is
        sigma2_t = 1.0
    return sigma2_t

T = int(snakemake.wildcards.t)
pve = float(snakemake.wildcards.pve) / 100
print('generating data')
print('\tT={}, pve={}'.format(T, pve))
data, info = make_simulation(X, T=T, pve=pve, sparsity=0.1)
model = M(**data, K=10)

print('fitting full model')
fit_args = {
    'max_iter': 300,
    'update_covariate_weights': False,
    'update_weights': True,
    'update_pi': True,
    'ARD_weights': True,
    'update_variance': True,
    'verbose': True
}
model.fit(**fit_args)

#save_model
print('saving model')
compute_records(model)
strip_and_dump(model, snakemake.output.model)

pickle.dump(info, open(snakemake.output.info, 'wb'))

base_path = snakemake.output[0][:-len('.model')]
print('generating scores and variant file')
try:
    component_scores(model).to_json('{}.scores'.format(base_path))
    gene = base_path.split('/')[-2]
    make_variant_report(model, gene).to_csv('{}.variants.bed'.format(base_path), sep='\t')
except Exception:
    print('There was an error generating secondary files')
