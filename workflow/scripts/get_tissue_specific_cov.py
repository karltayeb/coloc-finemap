import pickle
import pandas as pd
import numpy as np
import os

# load data
ld = pd.read_csv(snakemake.input.ld, sep='\t', header=None).values
associations = pd.read_csv(snakemake.input.associations, sep='\t', index_col=0)
snplist = np.squeeze(pd.read_csv(snakemake.input.snps, header=None).values)

# filter out nans
associations = associations.loc[:, ~np.any(np.isnan(associations), axis=0)]
mask = ~np.any(np.isnan(ld), axis=1)
ld = ld[mask][:, mask]
snplist = snplist[mask]

# filter to intersection
intersect = np.intersect1d(snplist, associations.columns.values)
associations = associations.loc[:, intersect]
mask = np.isin(snplist, intersect)
ld = ld[mask][:, mask]
snplist = snplist[mask]

alpha = 0.05
LD = np.array([
    (1 - alpha) * data['X'] + alpha * np.outer(y, y) for y in data['Y']
])

data = {
    'X': LD,
    'Y': associations.values,
    'tissue_ids': associations.index.values,
    'snp_ids': associations.columns.values
}

pickle.dump(data, open(snakemake.output[0], 'wb'))
