import pickle
import pandas as pd
import numpy as np
from coloc.independent_model import IndependentFactorSER
import snakemake

data = pickle.load(open(snakemake.input[0], 'rb'))
g = IndependentFactorSER(X=data['X'], Y=data['Y'], K=10)
g.fit(max_iter=100, update_active=False, update_weights=True, update_pi=True, 
      ARD_weights=True, update_variance=True, verbose=True)
path = '/'.join(snakemake.output[0].split('/')[:-1])
name = snakemake.output[0].split('/')[-1]
g.save(path, name)
