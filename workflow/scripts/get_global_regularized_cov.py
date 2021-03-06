import pickle
import numpy as np

data = pickle.load(open(snakemake.input[0], 'rb'))
#alpha = snakemake.params.alpha
alpha = 0.1
zscore_cov = np.cov(data['Y'].T)
data['X'] = (1 - alpha) * data['X'] + alpha * zscore_cov
pickle.dump(data, open(snakemake.output[0], 'wb'))
