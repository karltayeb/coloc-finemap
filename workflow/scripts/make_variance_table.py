import pickle
import itertools
import pandas as pd
import numpy as np
import glob
from utils import load_model

data = pickle.load(
    open(snakemake.input.data, 'rb')
)

model_dict = pickle.load(
    open(snakemake.input.model_path, 'rb')
)

ard_var = model_dict['prior_varance']

table = []
labels = (data['true_effects'] @ data['true_effects'].T !=0)
shared = ((data['true_effects']!=0).astype(int) @ (data['true_effects'].T !=0).astype(int))

for t1, t2 in itertools.combinations(np.arange(ard_var.shape[0]), 2):
    max_min_se = np.max(np.minimum(
        ard_var[t1], ard_var[t2]
    ))
    table.append({
        't1': t1,
        't2': t2,
        'max_min_se': max_min_se,
        'shared_causal_snps': shared[t1, t2],
        'label': labels[t1, t2]
    })
table = pd.DataFrame(table)
table.to_csv(snakemake.output[0])