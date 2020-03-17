import pickle
import numpy as np

data = pickle.load(open(snakemake.input[0], 'rb'))
#alpha = snakemake.params.alpha
alpha = 0.1
"""
LD = np.array([
    (1 - alpha) * data['LD'] + alpha * np.outer(y, y) for y in data['zscores']
])
"""
zscore_cov = np.cov(data['zscores'].T)
data['LD'] = (1 - alpha) * data['LD'] + alpha * zscore_cov
pickle.dump(data, open(snakemake.output[0], 'wb'))
