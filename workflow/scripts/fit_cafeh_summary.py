from coloc.ard_ser import MVNFactorSER
import pandas as pd
import numpy as np
import pickle

data = pickle.load(open(snakemake.input[0], 'rb'))

Y = data['zscores']
if 't1' in snakemake.output[0]:
    t1 = int(snakemake.wildcards.tissue1)
    t2 = int(snakemake.wildcards.tissue2)
    Y = Y[[t1, t2]]

kwargs = {
    'snp_ids': data['variant_ids'],
    'tissue_ids': data['tissue_ids']
}

LD = data['LD']
n = MVNFactorSER(X=LD, Y=Y, K=20, **kwargs)
n.fit(max_iter=200, update_active=False, update_weights=True, update_pi=True, ARD_weights=True, verbose=True)
path = '/'.join(snakemake.output[0].split('/')[:-1])
name = snakemake.output[0].split('/')[-1]
n.save(path, name)