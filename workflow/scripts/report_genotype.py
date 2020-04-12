import numpy as np
import pandas as pd
import pickle
from coloc.independent_model import IndependentFactorSER

import os
import glob

import matplotlib.pyplot as plt

def load_compact_model():
    #get data
    model_dict = pickle.load(open(
        snakemake.input.model, 'rb'))

    #load genotype
    genotype = pd.read_csv(snakemake.input.genotype, sep=' ')
    genotype = genotype.set_index('IID').iloc[:, 5:]

    # mean imputation
    genotype = (genotype - genotype.mean(0))
    genotype = genotype.fillna(0)
    genotype = genotype.iloc[:, model_dict['snps_in_cs']]

    # load expression
    expression = pd.read_csv(snakemake.input.expression, sep='\t', index_col=0)

    # filter down to relevant individuals
    genotype = genotype.loc[expression.columns]

    # filter down to common individuals
    individuals = np.intersect1d(genotype.index.values, expression.columns.values)
    genotype = genotype.loc[individuals]
    expression = expression.loc[:, individuals]

    data = {
        'X': genotype.values.T,
        'Y': expression.values
    }

    # load model
    model = IndependentFactorSER(**data, K=1)
    model_dict['dims']['N'] = model.dims['N']
    model_dict['pi'] = model_dict['pi'][:, model_dict['snps_in_cs']]
    model_dict['snp_ids'] = model_dict['snp_ids'][model_dict['snps_in_cs']]
    model_dict['prior_pi'] = model_dict['prior_pi'][model_dict['snps_in_cs']]
    model_dict.pop('precompute', None)
    model.__dict__.update(model_dict)
    return model

def report_component_scores(model):
    active = np.array([model.purity[k] > 0.1 for k in range(model.dims['K'])])
    if active.sum() > 0:
        mw = model.weight_means
        mv = model.weight_vars
        pi = model.pi
        scores = np.einsum('ijk,jk->ij', np.abs(mw) / np.sqrt(mv), model.pi)
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
    weight_json = weights.to_json()
    with open(snakemake.output.scores, 'w') as f:
        f.write(weight_json)

def report_credible_set(model):
    active = np.array([model.purity[k] > 0.1 for k in range(model.dims['K'])])

    if active.sum() > 0:
        pi = pd.DataFrame(model.pi.T, index=model.snp_ids)
        min_cset_alpha = pd.concat(
            [pi.iloc[:, k].sort_values(ascending=False).cumsum() - pi.iloc[:, k]
             for k in np.arange(model.dims['K']) if model.purity[k] > 0.5],
            sort=False, axis=1
        ).min(1)
    else:
        min_cset_alpha = []

    gene = snakemake.output.csets.split('/')[-2]
    with open(snakemake.output.csets, 'w') as f:
        for row in min_cset_alpha.reset_index().values:
            variant, val = row
            chrom, pos = variant.split('_')[:2]
            pos = int(pos)
            line = '{}\t{}\t{}\t{}\t{}'.format(chrom, pos, pos+1, gene, val)
            print(line, file=f)

def report_expected_weights(model, path, gene):
    active = np.array([model.purity[k] > 0.1 for k in range(model.dims['K'])])
    weights = pd.DataFrame(
        model.get_expected_weights()[:, active],
        index = model.tissue_ids,
        columns = np.arange(model.dims['K'])[active]
    )
    weight_json = weights.to_json()
    with open('{}/genotype.expected_weights'.format(path), 'w') as f:
        f.write(weight_json)

def report_ard_precision(model, path, gene):
    active = np.array([model.purity[k] > 0.1 for k in range(model.dims['K'])])
    weights = pd.DataFrame(
        model.prior_precision[:, active],
        index = model.tissue_ids,
        columns = np.arange(model.dims['K'])[active]
    )

    weight_json = weights.to_json()
    with open('{}/genotype.ard_precision'.format(path), 'w') as f:
        f.write(weight_json)
        
def report_ard_precision(model, path, gene):
    active = np.array([model.purity[k] > 0.1 for k in range(model.dims['K'])])
    weights = pd.DataFrame(
        model.prior_precision[:, active],
        index = model.tissue_ids,
        columns = np.arange(model.dims['K'])[active]
    )

    weight_json = weights.to_json()
    with open('{}/genotype.ard_precision'.format(path), 'w') as f:
        f.write(weight_json)


model = load_compact_model()
report_credible_set(model)
report_component_scores(model)
