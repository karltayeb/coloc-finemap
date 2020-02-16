import pickle
from coloc.independent_model import IndependentFactorSER

data = pickle.load(open(snakemake.input[0], 'rb'))

Y = data['Y']
if 't1' in snakemake.output[0]:
    t1 = int(snakemake.wildcards.tissue1)
    t2 = int(snakemake.wildcards.tissue2)
    Y = Y[[t1, t2]]

kwargs = {
    'covariates': data['covariates'],
    'snp_ids': data['variant_ids'],
    'tissue_ids': data['tissue_ids'],
    'sample_ids': data['sample_ids']
}

g = IndependentFactorSER(X=data['X'], Y=Y, K=10, **kwargs)
g.fit(max_iter=100, update_active=False, update_weights=True, update_pi=True, 
      ARD_weights=True, update_variance=True, verbose=True)
path = '/'.join(snakemake.output[0].split('/')[:-1])
name = snakemake.output[0].split('/')[-1]
g.save(path, name)
