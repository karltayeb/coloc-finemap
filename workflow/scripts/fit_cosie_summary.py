from coloc.ard_ser import MVNFactorSER
import pandas as pd
import numpy as np
import pickle
import pdb; pdb.set_trace()

data = pickle.load(open(snakemake.input[0], 'rb'))
n = MVNFactorSER(X=data['LD'], Y=data['zscores'], K=10)
n.fit(max_iter=200, update_active=False, update_weights=True, update_pi=True, ARD_weights=True, verbose=True)
path = '/'.join(snakemake.output[0].split('/')[:-1])
name = snakemake.output[0].split('/')[-1]
n.save(path, name)
