import numpy as np
import pandas as pd
from coloc.independent_model2 import IndependentFactorSER as M
from coloc.independent_model_summary import IndependentFactorSER as M2
from coloc.misc import *
import scipy as sp

def make_simulation(genotype, T, pve, sparsity):
    # select snps
    NN = genotype.shape[1]
    start = np.random.choice(NN-1000)
    X = genotype.values.T[start:start+1000]
    snp_ids = genotype.columns.values[start:start+1000]

    N, M = X.shape
    causal_snps = np.random.choice(N, 5)
    true_effects = np.zeros((T, N))
    true_effects[:, causal_snps] = \
        np.random.binomial(1, sparsity, ((T, 5))) \

    if snakemake.params.sample_effects:
        true_effects = true_effects \
            * np.random.normal(size=true_effects.shape)
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
        'snp_ids': snp_ids,
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
pve = float(snakemake.wildcards.pve) / 100
print('generating data')
print('\tT={}, pve={}'.format(T, pve))
data, info = make_simulation(genotype, T=T, pve=pve, sparsity=0.1)

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
    'verbose': True
}
model.fit(**fit_args)

#save_model
base_path = snakemake.output.model
print('generating scores and variant file')
try:
    component_scores(model).to_json('{}.scores'.format(base_path))
    gene = base_path.split('/')[-2]
    make_variant_report(model, gene).to_csv('{}.variants.bed'.format(base_path), sep='\t')
except Exception:
    print('There was an error generating secondary files')

print('saving model')
compute_records(model)
strip_and_dump(model, snakemake.output.model, save_data=False)

### RUN SUSIE
print('training susie')
susie = M(**data, K=100)
a = sp.linalg.block_diag(*[np.ones((5)) * 1e10 for t in range(20)])
a = 1 / (a + 1e-10)
susie.a = a

susie.b = np.ones((susie.dims['T'], susie.dims['K']))
susie.weight_precision_a = susie.a

fit_args = {
    'max_iter': 2,
    'update_covariate_weights': True,
    'update_weights': True,
    'update_pi': True,
    'ARD_weights': True,
    'update_variance': True,
    'verbose': True,
}
susie.fit(**fit_args)

fit_args = {
    'max_iter': 20,
    'update_covariate_weights': True,
    'update_weights': True,
    'update_pi': True,
    'ARD_weights': True,
    'update_variance': True,
    'verbose': True
}

for components in np.arange(5*T).reshape(-1, 5):
    susie.fit(**fit_args, components=components)

base_path = snakemake.output.susie
print('generating scores and variant file')
try:
    component_scores(susie).to_json('{}.scores'.format(base_path))
    gene = base_path.split('/')[-2]
    make_variant_report(susie, gene).to_csv('{}.variants.bed'.format(base_path), sep='\t')
except Exception:
    print('There was an error generating secondary files')

print('saving susie')
compute_records(susie)
strip_and_dump(susie, snakemake.output.susie, save_data=False)
