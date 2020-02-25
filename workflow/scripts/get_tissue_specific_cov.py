import pickle
import numpy as np

data = pickle.load(open(snakemake.input[0], 'rb'))
alpha = snakemake.params.alpha
if alpha is None:
    alpha = 1 / (data['X'].shape[1] + 1)

print('Making tissue specific covariance with alpha={}'.format(alpha))
LD = np.array([
    (1 - alpha) * data['LD'] + alpha * np.outer(y, y) for y in data['zscores']
])

data['LD'] = LD
pickle.dump(data, open(snakemake.output[0], 'wb'))
