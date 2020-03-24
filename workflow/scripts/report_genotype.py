import numpy as np
import pandas as pd
import pickle
from coloc.independent_model import IndependentFactorSER

import os
import glob

import matplotlib.pyplot as plt

def load_compact_model(path, gene):
    #get data
    model_dict = pickle.load(open(
        '{}/genotype.model'.format(path), 'rb'))

    cols = ['IID']
    cols.extend(list(model_dict['snp_ids']))
    genotype = pd.read_csv('{}/{}.raw'.format(path, gene), sep=' ', index_col=0, usecols=cols)
    genotype = genotype.loc[model_dict['sample_ids']]
    genotype = genotype.iloc[:, model_dict['snps_in_cs']]
    
    expression = pd.read_csv('{}/{}.expression'.format(path, gene), sep='\t', index_col=0)
    data = {
        'X': genotype.values.T,
        'Y': expression.values
    }

    # load model
    model = IndependentFactorSER(**data, K=20)
    model_dict.pop('dims')
    model_dict['pi'] = model_dict['pi'][:, model_dict['snps_in_cs']]
    model_dict['snp_ids'] = model_dict['snp_ids'][model_dict['snps_in_cs']]
    model_dict['prior_pi'] = model_dict['prior_pi'][model_dict['snps_in_cs']]
    model.__dict__.update(model_dict)
    return model


def report_expected_weights(model, path, gene):
    active = np.array([model.purity[k] > 0.1 for k in range(20)])
    weights = pd.DataFrame(
        model.get_expected_weights()[:, active],
        index = model.tissue_ids,
        columns = np.arange(20)[active]
    )
    weight_json = weights.to_json()
    with open('{}/genotype.expected_weights'.format(path), 'w') as f:
        f.write(weight_json)
        
        
def report_ard_precision(model, path, gene):
    active = np.array([model.purity[k] > 0.1 for k in range(20)])
    weights = pd.DataFrame(
        model.prior_precision[:, active],
        index = model.tissue_ids,
        columns = np.arange(20)[active]
    )

    weight_json = weights.to_json()
    with open('{}/genotype.ard_precision'.format(path), 'w') as f:
        f.write(weight_json)
        
def report_ard_precision(model, path, gene):
    active = np.array([model.purity[k] > 0.1 for k in range(20)])
    weights = pd.DataFrame(
        model.prior_precision[:, active],
        index = model.tissue_ids,
        columns = np.arange(20)[active]
    )

    weight_json = weights.to_json()
    with open('{}/genotype.ard_precision'.format(path), 'w') as f:
        f.write(weight_json)
        
def report_credible_set(model, path, gene):
    active = np.array([model.purity[k] > 0.1 for k in range(20)])
    pi = pd.DataFrame(model.pi.T, index=model.snp_ids)
    min_cset_alpha = pd.concat(
        [pi.iloc[:, k].sort_values(ascending=False).cumsum() for k in np.arange(20)[active]],
        sort=False, axis=1
    ).min(1)

    with open('{}/genotype.cset.bed'.format(path), 'w') as f:
        for row in min_cset_alpha.reset_index().values:
            variant, val = row
            chrom, pos = variant.split('_')[:2]
            pos = int(pos)
            line = '{}\t{}\t{}\t{}\t{}'.format(chrom, pos, pos+1, gene, val)
            print(line, file=f)

df = pd.read_csv('../../output/GTEx/protein_coding_autosomal_egenes.txt', sep='\t')
paths = df.apply(lambda x: '../../output/GTEx/{}/{}'.format(x.chromosome, x.gene), axis=1)

for i, path in enumerate(paths):
    gene = path.split('/')[-1]
    print('{}: {} ({:.2f}%)'.format(i, gene, 100 * float(i / paths.shape[0]))) 
    try:
        model = load_compact_model(path, gene)
        report_credible_set(model, path, gene)
        report_expected_weights(model, path, gene)
        report_ard_precision(model, path, gene)
    except Exception:
        print('\tsomething went wrong')
