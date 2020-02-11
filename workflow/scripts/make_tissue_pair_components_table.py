import pickle
import numpy as np
import pandas as pd
from scipy.stats import norm
from coloc.ard_ser import MVNFactorSER
from coloc.independent_model import IndependentFactorSER
from utils import *

def load_models_and_data(data_path, genotype_model_path, summary_model_path):
    """
    load models and data
    """
    genotype_model = IndependentFactorSER(np.zeros((1, 1)), np.zeros((1, 1)), 1)
    summary_model = MVNFactorSER(np.zeros((1, 1)), np.zeros((1, 1)), 1)

    data = pickle.load(open(data_path, 'rb'))
    assign(genotype_model, pickle.load(open(genotype_model_path, 'rb')))
    assign(summary_model, pickle.load(open(summary_model_path, 'rb')))

    genotype_model.X = data['X']
    genotype_model.Y = data['Y']

    summary_model.X = data['LD']
    summary_model.Y = data['zscores']
    return genotype_model, summary_model, data


if __name__ == '__main__':
    g, s, data = load_models_and_data(
        snakemake.input.data_path, snakemake.input.genotype_model_path, snakemake.input.summary_model_path)
    df = make_table(g, data)
    pairs = pair_coloc(df.loc[df.active==1])
    pairs.to_csv(snakemake.output.genotype_output, index=False, sep='\t')

    df = make_table(s, data)
    pairs = pair_coloc(df.loc[df.active==1])
    pairs.to_csv(snakemake.output.summary_output, index=False, sep='\t')