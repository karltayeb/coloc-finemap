import pickle
import numpy as np
import pandas as pd

def compute_sigma2(X, true_effect, pve):
    var = np.var(true_effect @ X)
    sigma2_t = var/pve - var
    if sigma2_t == 0:
        # if variance is 0, there were no causal variants-- dont care what the variance is
        sigma2_t = 1.0
    return sigma2_t


def get_cis_variants(genotype, tss, window=500000, size=1000):
    pos = np.array([int(x.split('_')[1]) for x in genotype.index.values])
    #cis_variants = genotype.iloc[np.abs(pos - tss) < window]
    cis_variants = genotype.iloc[
        np.arange(genotype.shape[0])[np.argsort(np.abs(pos - tss))[:size]]
    ]
    return cis_variants

def generate_data(X, linkage=0.0, pve=0.1):
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
    T = 7
    # normalize and compute LD
    X = (X - X.mean(1)[:, None]).values / np.sqrt(np.var(X, 1))[:, None]
    LD = X @ X.T / M

    # select snps
    snp1 = np.random.choice(N)
    snp2 = np.argmin(np.abs((np.abs((LD + np.eye(LD.shape[0])*10)[snp1]) - linkage)))
    causal_snps = np.array([snp1, snp2])
    true_effects = np.zeros((7, X.shape[0]))
    true_effects[[1, 2, 3], snp1] = 1.0
    true_effects[[4, 5, 6], snp2] = 1.0

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
        'tissue_variance': tissue_variance,
        'ld': LD[snp1, snp2]
    }
    return data


gene = 'ENSG00000167280.12'
linkage = 0.0
pve = 0.1
output = "output/simulation/single_causal_variant/{}/ld_{:.2f}_pve_{:.2f}_data".format(gene, linkage, pve)


gencode = pd.read_csv('/work-zfs/abattle4/lab_data/genomic_annotation_data/gencode.v19.genes.v6p.patched_contigs_TSS.bed', sep='\t')
gene = gencode.loc[gencode.iloc[:, 3] == gene]
chr_num = gene.iloc[0, 0]
tss = gene.iloc[0, 1]
gene_name = gene.iloc[0, 3]

cis_variants = pd.read_csv('output/genotypes/{}_cis_variants'.format(gene_name), index_col=0)
data = generate_data(cis_variants, float(linkage), float(pve))
pickle.dump(data, open(output, 'wb'))

