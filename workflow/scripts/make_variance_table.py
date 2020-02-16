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
    open(snakemake.input.model, 'rb')
)

import pdb; pdb.set_trace()
model = load_model(data, summary_model_path=snakemake.input.model)
_, purity = model.get_credible_sets()
active = np.array([k for k in range(model.dims['K'])if purity[k] > 0.7])

ard_var = 1 / model_dict['prior_precision']

data['pi']

#todo: actually load the model, filter components
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
table.to_csv(snakemake.output[0], sep='\t')
