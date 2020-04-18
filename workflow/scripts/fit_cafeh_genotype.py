import pandas as pd
import pickle
import numpy as np
from coloc.independent_model2 import IndependentFactorSER

def strip_and_dump(model, path):
    """
    save the model with data a weight parameters removed
    add 'mini_weight_measn' and 'mini_weight_vars' to model
    the model can be reconstituted in a small number of iterations
    """
    PIP = 1 - np.exp(np.log(1 - model.pi + 1e-10).sum(0))
    mask = (PIP > 0.01)
    wv = model.weight_vars[:, :, mask]
    wm = model.weight_means[:, :, mask]

    credible_sets, purity = model.get_credible_sets(0.99)
    active = np.array([purity[k] > 0.5 for k in range(model.dims['K'])])
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

    model.__dict__.pop('precompute', None)
    model.__dict__.pop('weight_vars', None)
    model.__dict__.pop('weight_means', None)
    model.__dict__.pop('X', None)
    model.__dict__.pop('Y', None)
    model.__dict__.pop('covariates', None)
    pickle.dump(model, open(path, 'wb'))

def component_scores(model):
    purity = model.get_credible_sets(0.99)[1]
    active = np.array([purity[k] > 0.5 for k in range(model.dims['K'])])
    if active.sum() > 0:
        mw = model.weight_means
        mv = model.weight_vars
        pi = model.pi
        scores = np.einsum('ijk,jk->ij', mw / np.sqrt(mv), model.pi)
        weights = pd.DataFrame(
            scores[:, active],
            index = model.tissue_ids,
            columns = np.arange(model.dims['K'])[active]
        )
    else:
        weights = pd.DataFrame(
            np.zeros((model.dims['T'], 1)),
            index = model.tissue_ids
        )
    return weights

def make_variant_report(model, gene):
    PIP = 1 - np.exp(np.log(1 - model.pi + 1e-10).sum(0))
    purity = model.get_credible_sets(0.99)[1]
    active = np.array([purity[k] > 0.5 for k in range(model.dims['K'])])

    if active.sum() == 0:
        active[0] = True

    pi = pd.DataFrame(model.pi.T, index=model.snp_ids)
    cset_alpha = pd.concat(
        [pi.iloc[:, k].sort_values(ascending=False).cumsum() - pi.iloc[:, k]
         for k in np.arange(model.dims['K']) if purity[k] > 0.5],
        sort=False, axis=1
    )

    most_likely_k = np.argmax(pi.values[:, active], axis=1)
    most_likely_p = np.max(pi.values[:, active], axis=1)
    most_likely_cset_alpha = cset_alpha.values[np.arange(pi.shape[0]), most_likely_k]

    A = pd.DataFrame(
        [PIP, most_likely_k, most_likely_p, most_likely_cset_alpha],
        index=['PIP','k', 'p', 'min_alpha'], columns=model.snp_ids).T

    A.loc[:, 'chr'] = [x.split('_')[0] for x in A.index]
    A.loc[:, 'start'] = [int(x.split('_')[1]) for x in A.index]
    A.loc[:, 'end'] = A.loc[:, 'start'] + 1
    A.reset_index(inplace=True)
    A.loc[:, 'variant_id'] = A.loc[:, 'index'].apply(lambda x: '_'.join(x.split('_')[:-1]))
    A.loc[:, 'ref'] = A.loc[:, 'index'].apply(lambda x: x.split('_')[-1])


    A.loc[:, 'gene_id'] = gene
    A = A.set_index(['chr', 'start', 'end'])
    return A.loc[:, ['gene_id', 'variant_id', 'ref', 'PIP', 'k', 'p', 'min_alpha']]

data = pickle.load(open(snakemake.input[0], 'rb'))
model = IndependentFactorSER(**data, K=snakemake.params.k)

fit_args = {
    'max_iter': 3,
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
model.fit(**fit_args)

base_path = snakemake.output[0][:-len('.model')]

print('generating scores and variant file')
try:
    component_scores(model).to_json('{}.scores'.format(base_path))
    gene = base_path.split('/')[-2]
    make_variant_report(model, gene).to_csv('{}.variants.bed'.format(base_path), sep='\t')
except Exception:
    print('There was an error generating secondary files')

print('saving model')
strip_and_dump(model, snakemake.output[0])
