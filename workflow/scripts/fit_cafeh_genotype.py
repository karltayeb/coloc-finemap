import pickle
import numpy as np
from coloc.independent_model import IndependentFactorSER

"""
data = pickle.load(open(snakemake.input[0], 'rb'))

Y = data['Y']
if 't1' in snakemake.output[0]:
    t1 = int(snakemake.wildcards.tissue1)
    t2 = int(snakemake.wildcards.tissue2)
    Y = Y[[t1, t2]]

kwargs = {
    'covariates': data['covariates'] if 'covariates' in data else None,
    'snp_ids': data['variant_ids'] if 'variant_ids' in data else None,
    'tissue_ids': data['tissue_ids'] if 'tissue_ids' in data else None,
    'sample_ids': data['sample_ids'] if 'sample_ids' in data else None
}

g = IndependentFactorSER(X=data['X'], Y=Y, K=20, **kwargs)
g.fit(max_iter=100, update_active=False, update_weights=True, update_pi=True, 
      ARD_weights=True, update_variance=True, verbose=True)
path = '/'.join(snakemake.output[0].split('/')[:-1])
name = snakemake.output[0].split('/')[-1]
g.save(path, name)
"""
data = pickle.load(open(snakemake.input[0], 'rb'))
model = IndependentFactorSER(**data, K=40)
model.fit(max_iter=500, verbose=True, ARD_weights=True)

# get broad cs-- contains most information for model
cs, p = model.get_credible_sets(0.999)
active = np.array([p[k] > 0.01 for k in range(20)])

# boolean mask for relevant snps
snps_in_cs = np.unique(np.concatenate([cs[k] for k in range(20) if active[k]]))
snps_in_cs = np.isin(model.snp_ids, snps_in_cs)

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

save_dict['weight_means'] = wm
save_dict['weight_vars'] = wv
save_dict['snps_in_cs'] = snps_in_cs

# save pickle
pickle.dump(save_dict, open(snakemake.output[0], 'wb'))
