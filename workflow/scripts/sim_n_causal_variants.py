import numpy as np
import pandas as pd
from coloc.independent_model2 import IndependentFactorSER as M
from coloc.independent_model_summary import IndependentFactorSER as M2
import scipy as sp
import pickle

def make_simulation(genotype, U, V, T, pve):
    """
    genotype is genotype
    U is the number of causal variants per tissue
    V is the total number of causal variants
    T is the number of tissues
    pve is the percent variance explained by genotype
    """
    print('generating data')
    print('\tT={}, \tU={}, \tV={},\tpve={}'.format(T, U, V, pve))
    
    # select snps
    NN = genotype.shape[1]
    start = np.random.choice(NN-1000)
    X = genotype.values.T[start:start+1000]
    snp_ids = genotype.columns.values[start:start+1000]

    N, M = X.shape
    # sample causal snps
    causal_snps = np.random.choice(N, V)
    true_effects = np.zeros((T, N))

    for t in range(T):
        # select U causal snps to use in tissue t
        causal_in_t = np.random.choice(causal_snps, U)
        true_effects[t, causal_in_t] = np.random.normal(size=causal_in_t.size)

    tissue_variance = np.array([
        compute_sigma2(X, te, pve) for te in true_effects
    ])

    #simulate expression
    expression = (true_effects @ X) + \
        np.random.normal(size=(T, M)) * np.sqrt(tissue_variance)[:, None]
    expression = expression - expression.mean(1)[:, None]

    # simulate expression
    data = {
        'X': X,
        'Y': expression,
        'snp_ids': snp_ids,
    }

    sim_info = {
        'causal_snps': causal_snps[np.abs(true_effects[:, causal_snps]).sum(0) > 0],
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


def compute_records(model):
    """
    save the model with data a weight parameters removed
    add 'mini_weight_measn' and 'mini_weight_vars' to model
    the model can be reconstituted in a small number of iterations
    """
    credible_sets, purity = model.get_credible_sets(0.999)
    active = np.array([purity[k] > 0.1 for k in range(model.dims['K'])])

    try:
        snps = np.unique(np.concatenate([
            credible_sets[k] for k in range(model.dims['K']) if active[k]]))
    except Exception:
        snps = np.unique(np.concatenate([
            credible_sets[k][:5] for k in range(model.dims['K'])]))
    mask = np.isin(model.snp_ids, snps)

    wv = model.weight_vars[:, :, mask]
    wm = model.weight_means[:, :, mask]

    records = {
        'active': active,
        'purity': purity,
        'credible_sets': credible_sets,
        'EXz': model.pi @ model.X,
        'mini_wm': wm,
        'mini_wv': wv,
        'snp_subset': mask
    }
    model.records = records


def strip_and_dump(model, path, save_data=False):
    """
    save the model with data a weight parameters removed
    add 'mini_weight_measn' and 'mini_weight_vars' to model
    the model can be reconstituted in a small number of iterations
    """
    compute_records(model)

    # purge precompute
    for key in model.precompute:
        model.precompute[key] = {}
    model.__dict__.pop('weight_means', None)
    model.__dict__.pop('weight_vars', None)
    if not save_data:
        model.__dict__.pop('X', None)
        model.__dict__.pop('Y', None)
        model.__dict__.pop('covariates', None)
    pickle.dump(model, open(path, 'wb'))

######### GET GENOTYPE
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

T = int(snakemake.wildcards.t)
U = int(snakemake.wildcards.snps_per_tissue)
V = U*3

pve = float(snakemake.wildcards.pve) / 100
data, info = make_simulation(genotype, U, V, T, pve)

pickle.dump(info, open(snakemake.output.info, 'wb'))
pickle.dump(data, open(snakemake.output.data, 'wb'))


##### TRAIN CAFEH GENOTYPE
model = M(**data, K=10)
print('fitting full model')
fit_args = {
    'max_iter': 300,
    'update_covariate_weights': False,
    'update_weights': True,
    'update_pi': True,
    'ARD_weights': True,
    'update_variance': True,
    'verbose': False
}
model.fit(**fit_args)
print('model fit:\n\titers:{}\n\tELBO:{}\n\trun-time:{}'.format(len(model.elbos), model.elbos[-1], model.run_time))

strip_and_dump(model, snakemake.output.model, False)

