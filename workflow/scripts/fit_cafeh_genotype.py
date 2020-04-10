import pickle
import numpy as np
from coloc.independent_model2 import IndependentFactorSER

data = pickle.load(open(snakemake.input[0], 'rb'))
model = IndependentFactorSER(**data, K=snakemake.params.k)
model.fit(
    max_iter=500,
    update_covariate_weights=True,
    update_weights=True,
    update_pi=True,
    ARD_weights=True,
    update_variance=True,
    verbose=True
)

PIP = 1 - np.exp(np.log(1 - model.pi + 1e-10).sum(0))
snps_in_cs = (PIP > 1e-2)

# only retain relevant snps
wm = model.weight_means[:, :, snps_in_cs]
wv = model.weight_vars[:, :, snps_in_cs]
pi = model.pi[:, snps_in_cs]

# make save dict
save_dict = model.__dict__
save_dict.pop('X')
save_dict.pop('Y')
save_dict.pop('covariates')
save_dict.pop('sample_covariate_map')
save_dict.pop('precompute')

save_dict['weight_means'] = wm
save_dict['weight_vars'] = wv
save_dict['snps_in_cs'] = snps_in_cs

# save pickle
pickle.dump(save_dict, open(snakemake.output[0], 'wb'))
