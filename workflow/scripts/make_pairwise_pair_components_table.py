import pickle
import numpy as np
import pandas as pd
from coloc.ard_ser import MVNFactorSER
from coloc.independent_model import IndependentFactorSER
from utils import *

def load_data(data_path):
    """
    load models and data
    """
    return pickle.load(open(data_path, 'rb'))


def load_model(data, genotype_model_path=None, summary_model_path=None):
    path = genotype_model_path or summary_model_path
    t1 = int(path.split('/')[-1].split('_')[1])
    t2 = int(path.split('/')[-1].split('_')[3])

    sub_data = {
        'X': data['X'],
        'LD': data['LD'],
        'Y': data['Y'][[t1, t2]],
        'betas': data['betas'][[t1, t2]],
        'standard_errors': data['standard_errors'][[t1, t2]],
        'zscores': data['zscores'][[t1, t2]],
        'causal_snps': data['causal_snps'],
        'true_effects': data['true_effects'][[t1, t2]],
        'tissue_variance': data['tissue_variance'][[t1, t2]],
        'ld' : data['ld']
    }

    if genotype_model_path is not None:
        model = IndependentFactorSER(np.zeros((1, 1)), np.zeros((1, 1)), 1)
        assign(model, pickle.load(open(path, 'rb')))
        model.X = sub_data['X']
        model.Y = sub_data['Y']
    else:
        model = MVNFactorSER(np.zeros((1, 1)), np.zeros((1, 1)), 1)
        assign(model, pickle.load(open(path, 'rb')))
        model.X = sub_data['X']
        model.Y = sub_data['zscores']
    return model, sub_data

data_path = 'output/simulation/single_causal_variant/pve_0.1/ld_0.8/gene_ENSG00000262000.1/data'
summary_model_paths = [
    'output/simulation/single_causal_variant/pve_0.1/ld_0.8/gene_ENSG00000262000.1/pairwise_summary/t1_0_t2_2_model_summary',
    'output/simulation/single_causal_variant/pve_0.1/ld_0.8/gene_ENSG00000262000.1/pairwise_summary/t1_0_t2_5_model_summary',
    'output/simulation/single_causal_variant/pve_0.1/ld_0.8/gene_ENSG00000262000.1/pairwise_summary/t1_1_t2_2_model_summary',
    'output/simulation/single_causal_variant/pve_0.1/ld_0.8/gene_ENSG00000262000.1/pairwise_summary/t1_1_t2_5_model_summary',
    'output/simulation/single_causal_variant/pve_0.1/ld_0.8/gene_ENSG00000262000.1/pairwise_summary/t1_4_t2_2_model_summary',
    'output/simulation/single_causal_variant/pve_0.1/ld_0.8/gene_ENSG00000262000.1/pairwise_summary/t1_4_t2_5_model_summary']
output = 'output/simulation/single_causal_variant/pve_0.1/ld_0.8/gene_ENSG00000262000.1/pairwise_summary/pairs_summary'

summary_pairs = []
data = load_data(data_path)
key = '/'.join(data_path.split('/')[5:-1])

sub_summary_paths = [x for x in snakemake.input.summary_model_paths if key in x]
for summary_model_path in sub_summary_paths:
    model, sub_data = load_model(data, summary_model_path=summary_model_path)
    df = make_table(model, sub_data)
    pairs = pair_coloc(df.loc[df.active == 1])
    if pairs.size > 0:
        summary_pairs.append(pairs)
summary_pairs = pd.concat(summary_pairs)
summary_pairs.to_csv(snakemake.output.summary_output, index=False, sep='\t')
