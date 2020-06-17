import pandas as pd
import pickle
import numpy as np
from coloc.independent_model_ss import IndependentFactorSER as GSS
from coloc.cafeh_ss import CAFEH as CSS
from coloc.misc import *
from coloc.simulation import *
from coloc.covariance import *
from itertools import product

superpop2samples = pickle.load(open('output/superpop2samples_1kG', 'rb'))

ld_functions = {
    'sample': sample_ld,
    'eur_ld': eur_ld,
    'asn_ld': asn_ld,
    'afr_ld': afr_ld,
    'reference': refernce_ld,
    'z': z_ld,
    'lw_sample': ledoit_wolf_sample_ld,
    'lw_refence': ledoit_wolf_reference_ld,
    'lw_z': ledoit_wolf_z_ld,
    'ref_z': ref_z_ld,
    'z3': z3_ld
}

def smooth_betas(data, ld, epsilon=1.0):
    """
    return a copy of data with smoothed effect sizes
    beta_sooth = SRS(SRS + epsilonS^2)^{-1} beta
    """
    Bs = []
    R = ld_functions['ld'](data)
    for i in range(data.S.shape[0]):
        S = np.diag(data.S.iloc[i].values)
        B = data.B.iloc[i].values
        SRS = S @ R @ S
        Bs.append(SRS @ np.linalg.solve(SRS + epsilon * S**2, B))
    data_smooth = deepcopy(data)
    data_smooth.B = Bs
    return data_smooth

def init_css(data, K=10, ld='sample', pi0=1.0, dispersion=1.0, epsilon=0.0):
    # TODO epsilon--- smooth expression?
    init_args = {
        'LD': ld_functions[ld](data),
        'B': data.B.values,
        'S': data.S.values,
        'K': K,
        'snp_ids': data.B.columns.values,
        'tissue_ids': data.B.index.values
    }
    name = 'simid-{}_gene-{}_k-{}_pi0-{}_d-{}_e-{}_ld-{}.css'.format(
        data.id, data.gene, K, pi0, dispersion, epsilon, ld)
    print('initializing summary stat model')
    css = CSS(**init_args)
    css.prior_activity = np.ones(K) * pi0
    css.tissue_precision_b = np.ones(css.dims['T']) * dispersion
    css.name = name
    return css

def init_gss(data, K=10, p=1.0):
    print('initializing genotype model')
    name = 'sim-{}_gene-{}_k-{}_pi-{}_ld-{}.gss'.format(
        data.id, data.gene, K, p, ld_type)
    gss = GSS(
        X=center_mean_impute(data.genotype_gtex).values.T,
        Y=data.expression.values,
        K=K,
        covariates=None,
        snp_ids=data.common_snps,
        tissue_ids=data.expression.index.values,
        sample_ids=data.expression.columns.values
    )
    gss.prior_activity = np.ones(K) * p
    gss.name = name
    return gss

ld_functions = {
    'sample': sample_ld,
    'eur1kG': eur_ld,
    'asn1kG': asn_ld,
    'afr1kG': afr_ld,
    'reference1kG': refernce_ld,
    'z': z_ld,
    'lwrefence': ledoit_wolf_reference_ld,
    'lwz': ledoit_wolf_z_ld,
    'refz': ref_z_ld,
    'z3': z3_ld
}

#gss = load(snakemake.input[0])
gene = snakemake.wildcards.gene
sim_spec = pd.read_csv('output/sim/ld/sim_spec.txt', sep='\t')
sim_data, sim_expression = load_sim_from_model_data(gene, sim_spec)

# make model_spec
pi0s = [0.1]
Ks = [10]
ld_types = list(ld_functions.keys())
dispersions = [0.5, 1.0, 5.0]
epsilons = [0.0]
model_spec = pd.DataFrame(
    list(product(ld_types, Ks, pi0s, dispersions, epsilons)),
    columns=['ld', 'K', 'pi0', 'dispersion', 'epsilon'])

# fit CSS to simulation data
fit_args = {
    'update_weights': True,
    'update_pi': True,
    'ARD_weights': True,
    'update_variance': False,
    'verbose': False,
    'max_iter': 50
}

pickle.dump(sim_expression, open(snakemake.output[0], 'wb'))

bp = '/'.join(snakemake.output[0].split('/')[:-1])
for i, row in model_spec.iterrows():
    print('model spec:')
    for arg in row.to_dict():
        print('\t{}: {}'.format(arg, row.to_dict()[arg]))

    #TODO dont fit if we already have it
    name = 'simid-{}_gene-{}_k-{}_pi0-{}_d-{}_e-{}_ld-{}.css'.format(
        sim_data.id, sim_data.gene,
        row.to_dict()['K'], row.to_dict()['pi0'],
        row.to_dict()['dispersion'], row.to_dict()['epsilon'],
        row.to_dict()['ld'])
    save_path = bp + '/' + name
    if not os.path.isfile(save_path):
        css = init_css(sim_data, **row.to_dict())
        print('fitting model')
        css.fit(**fit_args, update_active=False)
        css.fit(**fit_args, update_active=True)
        compute_records_css(css)
        print('saving model to {}'.format(save_path))
        strip_and_dump(css, save_path)
        rehydrate_model(css)
    else:
        print('{} alread fit'.format(name))

pi0s = [0.1]
Ks = [10]
ld_types = list(ld_functions.keys())
dispersions = [1.0]
epsilons = [0.1, 1.0, 10.0]
model_spec = pd.DataFrame(
    list(product(ld_types, Ks, pi0s, dispersions, epsilons)),
    columns=['ld', 'K', 'pi0', 'dispersion', 'epsilon'])


smoothed_data = {}
for epsilon in model_spec.epsilon.unique():
    smoothed_data[epsilon] = smooth_betas(sim_data, epsilon)

for i, row in model_spec.iterrows():
    print('model spec:')
    for arg in row.to_dict():
        print('\t{}: {}'.format(arg, row.to_dict()[arg]))

    name = 'simid-{}_gene-{}_k-{}_pi0-{}_d-{}_e-{}_ld-{}.css'.format(
        sim_data.id, sim_data.gene,
        row.to_dict()['K'], row.to_dict()['pi0'],
        row.to_dict()['dispersion'], row.to_dict()['epsilon'],
        row.to_dict()['ld'])
    save_path = bp + '/' + name
    if not os.path.isfile(save_path):
        css = init_css(smoothed_data[row.epsilon], **row.to_dict())
        print('fitting model')
        css.fit(**fit_args, update_active=False)
        css.fit(**fit_args, update_active=True)
        compute_records_css(css)
        print('saving model to {}'.format(save_path))
        strip_and_dump(css, save_path)
        rehydrate_model(css)
    else:
        print('{} alread fit'.format(name))

# TODO fit GSS

"""
# fit GSS to simulation data
gss_sim = init_gss(sim_data, 10, 0.1)
fit_args = {
    'max_iter': 300,
    'update_covariate_weights': True,
    'update_weights': True,
    'update_pi': True,
    'ARD_weights': True,
    'update_variance': True,
    'verbose': False
}
print('training model')
for arg in fit_args:
    print('\t{}: {}'.format(arg, fit_args[arg]))
gss_sim.fit(**fit_args, update_active=False)
gss_sim.fit(**fit_args, update_active=True)
compute_records_gss(gss_sim)
print('saving model to {}'.format(snakemake.output[0]))
strip_and_dump(gss_sim, snakemake.output[0], False)
"""
