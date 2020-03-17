from coloc.cafeh import CAFEH
import pandas as pd
import numpy as np
from scipy.linalg import block_diag
import pickle

data = pickle.load(open(snakemake.input[0], 'rb'))


X = data['X']
Y = data['Y']
tissue_ids = data['tissue_ids']

n = CAFEH(**data, K=5*tissue_ids.size)
prior_precision = 1 / (block_diag(*[np.ones(5) for _ in range(tissue_ids.size)]) + 1e-10)
n.prior_precision = prior_precision

n.fit(max_iter=200, update_active=False, update_weights=True, update_pi=True, ARD_weights=True, verbose=True)
path = '/'.join(snakemake.output[0].split('/')[:-1])
name = snakemake.output[0].split('/')[-1]
n.save(path, name)


for i, tissue_id in tissue_ids:
    n = CAFEH(X=X, Y=Y[i][None], tissue_ids=[tissue_id], snp_ids=data['snp_ids'], K=5)
    n.fit(max_iter=200, update_active=False, update_weights=True,
          update_pi=True, ARD_weights=True, verbose=True)
    path = '/'.join(snakemake.output[0].split('/')[:-1])
    name = snakemake.output[0].split('/')[-1]
    n.save(path, name)
