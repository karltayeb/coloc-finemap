import pickle
import numpy as np
import pandas as pd
from scipy.stats import norm
from coloc.ard_ser import MVNFactorSER
from coloc.independent_model import IndependentFactorSER
from make_tissue_pair_component import assign, load_models_and_data, matched_labels, make_table, pair_coloc

data_path = snakemake.input.data_path
genotype_model_path = snakemake.input.genotype_model_path
summary_model_path = snakemake.input.summary_model_path


def load_data(data_path):
    """
    load models and data
    """
    return pickle.load(open(data_path, 'rb'))


def load_model(data, genotype_model_path=None, summary_model_path=None):
    path = genotype_model_path or summary_model_path
    t1 = int(path.split('/')[-1].split('_')[1])
    t2 = int(path.split('/')[-1].split('_')[3])

    if genotype_model_path is not None:
        model = IndependentFactorSER(np.zeros((1, 1)), np.zeros((1, 1)), 1)
        assign(model, pickle.load(open(genotype_model_path, 'rb')))
        model.X = data['X']
        model.Y = data['Y'][[t1, t2]]
    else:
        model = MVNFactorSER(np.zeros((1, 1)), np.zeros((1, 1)), 1)
        assign(model, pickle.load(open(summary_model_path, 'rb')))
        model.X = data['X']
        model.Y = data['zscores'][[t1, t2]]
    return model

summary_keys = []
for data_path in snakemake.input.data_paths:
    data = load_data(data_path)
    key = '/'.join(data_path.split('/')[5:-1])

    sub_summary_paths = [x for x in snakemake.input.summary_model_paths if key in x]
    for summary_model_path in sub_summary_paths:
        model = load_model(data, summary_model_path=summary_model_path)
        df = make_table(model, data)
        pairs = pair_coloc(df.loc[df.active == 1])
        if pairs.size > 0:
            summary_pairs.append(pairs)
    summary_pairs = pd.concat(summary_pairs)
    summary_pairs.to_csv(snakemake.output.summary_output, index=False, sep='/t')
