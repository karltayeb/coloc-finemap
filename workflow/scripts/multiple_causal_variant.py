import pickle
import numpy as np
import pandas as pd
from utils import compute_sigma2, get_cis_variants

def generate_data(X, sparsity=0.1, pve=0.1):
    """
    X: genotype data
    pve: percent variance explained by genotype (if there are causal snps)

    simulate 7 tissues with 2 causal variants
    2 variants are in varying degree of LD with one another
    first tissue has no causal variants
    for multiple tissues we generate independent, normally distirbuted
    effect sizes for each tissue

    there are two causal variants approximately (ld)
    """
    N, M = X.shape
    T = 10
    # normalize and compute LD
    X = (X - X.mean(1)[:, None]).values / np.sqrt(np.var(X, 1))[:, None]
    LD = X @ X.T / M

    # select snps
    causal_snps = np.random.choice(N, 5)
    true_effects = np.zeros((10, X.shape[0]))
    true_effects[:, causal_snps] = \
        np.random.binomial(1, sparsity, ((10, 5))) * np.random.normal((10, 5))

    tissue_variance = np.array([
        compute_sigma2(X, te, pve) for te in true_effects])

    #simulate expression
    expression = (true_effects @ X) + \
        np.random.normal(size=(T, M)) * np.sqrt(tissue_variance)[:, None]
    expression = expression - expression.mean(1)[:, None]

    # generate summary statistics
    beta = expression @ X.T / M
    se = np.sqrt(tissue_variance[:, None] / M)
    data = {
        'LD': LD,
        'X': X,
        'Y': expression,
        'betas': beta,
        'standard_errors': se,
        'zscores': beta / se,
        'causal_snps': causal_snps,
        'true_effects': true_effects,
        'tissue_variance': tissue_variance
    }
    return data

gencode = pd.read_csv(
    '/work-zfs/abattle4/lab_data/genomic_annotation_data/'
    'gencode.v19.genes.v6p.patched_contigs_TSS.bed', sep='\t')
gene = gencode.loc[gencode.iloc[:, 3] == snakemake.wildcards.gene]
chr_num = gene.iloc[0, 0]
tss = gene.iloc[0, 1]
gene_name = gene.iloc[0, 3]

cis_variants = pd.read_csv(snakemake.input[0], index_col=0)
data = generate_data(cis_variants,
    float(snakemake.wildcards.sparsity),
    float(snakemake.wildcards.pve))
pickle.dump(data, open(snakemake.output[0], 'wb'))

