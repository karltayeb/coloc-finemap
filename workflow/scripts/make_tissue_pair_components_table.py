import pickle
import numpy as np
import pandas as pd
from scipy.stats import norm
from coloc.cafeh import CAFEH
from utils import *

def load_models_and_data(data_path, summary_model_path):
    """
    load models and data
    """
    summary_model = CAFEH(np.zeros((1, 1)), np.zeros((1, 1)), 1)
    data = pickle.load(open(data_path, 'rb'))
    assign(summary_model, pickle.load(open(summary_model_path, 'rb')))

    summary_model.X = data['LD']
    summary_model.Y = data['zscores']
    return summary_model, data

s, data = load_models_and_data(
    snakemake.input.data_path, snakemake.input.summary_model_path)

df = make_table(s, data)
pairs = pair_coloc(df, data)
pairs.to_csv(snakemake.output.summary_output, index=False, sep='\t')