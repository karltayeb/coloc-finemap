import pickle
import pandas as pd
import numpy as np
from sklearn import linear_model

data = pickle.load(open(snakemake.input[0], 'rb'))

genotype = pd.DataFrame(data['X'], columns=data['sample_ids'], index=data['variant_ids'])
covariates = data['covariates']
tissue_ids = data['tissue_ids']

LD  = []
for t_id in tissue_ids:
    lm = linear_model.LinearRegression()
    lm.fit(covariates[t_id].T,
           genotype.loc[:, genotype.columns.isin(covariates[t_id].columns)].T)
    lm.predict(covariates[t_id].T).shape
    LD_tissue = (genotype.loc[:, genotype.columns.isin(
        covariates[t_id].columns)].T - lm.predict(covariates[t_id].T)).corr()
    LD.append(LD_tissue.values)
LD = np.stack(LD)

data['LD'] = LD
pickle.dump(data, open(snakemake.output[0], 'wb'))