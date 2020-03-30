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

    cols = ['IID']
    cols.extend(list(model_dict['snp_ids']))
    genotype = pd.read_csv(snakemake.input.genotype, sep=' ', index_col=0, usecols=cols)
    genotype = genotype.loc[model_dict['sample_ids']]
    genotype = genotype.iloc[:, model_dict['snps_in_cs']]
    expression = pd.read_csv(snakemake.input.expression, sep='\t', index_col=0)
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
    model.__dict__.update(model_dict)
    return model

def report_component_scores(model):
    active = np.array([model.purity[k] > 0.1 for k in range(model.dims['K'])])
    mw = model.weight_means
    mv = model.weight_vars
    pi = model.pi
    scores = ((np.abs(mw) / np.sqrt(mv)) * pi[None]).sum(-1)

    if active.sum() > 0:
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

def report_component_scores(model):
    active = np.array([model.purity[k] > 0.1 for k in range(model.dims['K'])])
    if active.sum() > 0:
        mw = model.weight_means
        mv = model.weight_vars
        pi = model.pi

        diags = np.stack([model._get_diag(t) for t in range(model.dims['T'])])
        a = np.clip(model.prior_precision[:, :, None], 1e-10, 1e10)
        b = (diags / model.tissue_variance[:, None])[:, None]
        s2 = 1 / (a + b)
        scores = ((np.abs(mw) / np.sqrt(s2)) * pi[None]).sum(-1)
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
            [pi.iloc[:, k].sort_values(ascending=False).cumsum() for k in np.arange(model.dims['K'])[active]],
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
