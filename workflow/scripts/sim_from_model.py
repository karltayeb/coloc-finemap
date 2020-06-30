import pandas as pd
import pickle
import numpy as np
from coloc.independent_model_ss import IndependentFactorSER as GSS
from coloc.cafeh_ss import CAFEH as CSS
from coloc.misc import *
from coloc.simulation import *
from coloc.covariance import *
from itertools import product
import sys

with open(snakemake.output[0], 'w') as f:
    print('Begin logging...', file=f)

#log = open(snakemake.output[0], "a")
#sys.stdout = log

superpop2samples = pickle.load(open(
    '/work-zfs/abattle4/karl/cosie_analysis/output/superpop2samples_1kG', 'rb'))

def cov2corr(X):
    """
    scale covariance matrix to correlaton matrix
    """
    diag = np.sqrt(np.diag(X))
    return (1/diag[:, None]) * X * (1/diag[None])

sample_ld = lambda data: np.corrcoef(data.data.X.T)
refernce_ld = lambda data: np.corrcoef(data.data.X1kG.T)
z_ld = lambda data: np.corrcoef(data.summary_stats.B.T.values / np.sqrt(data.summary_stats.V.T.values))
ledoit_wolf_sample_ld = lambda data: cov2corr(covariance.ledoit_wolf(data.data.X)[0])
ledoit_wolf_reference_ld = lambda data: cov2corr(covariance.ledoit_wolf(data.data.X1kG)[0])
ledoit_wolf_z_ld = lambda data: cov2corr(covariance.ledoit_wolf(data.summary_stats.B.values / np.sqrt(data.summary_stats.V.values))[0])
z3_ld = lambda data: z_ld(data)**3

eur_ld = lambda data: np.corrcoef(
    center_mean_impute(
        data.data.genotype_1kG.loc[superpop2samples['EUR']]).values.T
    + np.random.normal(scale=1e-10, size=(1000, superpop2samples['EUR'].size))
)
asn_ld = lambda data: np.corrcoef(
    center_mean_impute(
        data.data.genotype_1kG.loc[superpop2samples['ASN']]).values.T
    + np.random.normal(scale=1e-10, size=(1000, superpop2samples['ASN'].size))

)
afr_ld = lambda data: np.corrcoef(
    center_mean_impute(
        data.data.genotype_1kG.loc[superpop2samples['AFR']]).values.T
    + np.random.normal(scale=1e-10, size=(1000, superpop2samples['AFR'].size))

)
sea_ld = lambda data: np.corrcoef(
    center_mean_impute(
        data.data.genotype_1kG.loc[superpop2samples['SEA']]).values.T
    + np.random.normal(scale=1e-10, size=(1000, superpop2samples['SEA'].size))

)
amr_ld = lambda data: np.corrcoef(
    center_mean_impute(
        data.data.genotype_1kG.loc[superpop2samples['AMR']]).values.T
    + np.random.normal(scale=1e-10, size=(1000, superpop2samples['AMR'].size))
)

def ref_z_ld(data, alpha=None):
    """
    mix reference ld and zscore ld
    """
    if alpha is None:
        alpha = data.data.X1kG.shape[0] / (data.data.X1kG.shape[0] + data.summary_stats.B.shape[0])
    return alpha * refernce_ld(data) + (1 - alpha) * z_ld(data)

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

def average_ld(data, ld1, ld2, alpha=None):
    if alpha is None:
        alpha = data.data.X1kG.shape[0] / (data.data.X1kG.shape[0] + data.summary_stats.B.shape[0])
    return alpha * ld1(data) + (1 - alpha) * ld2(data)

def tril2symm(values):
    A = np.eye(1000) * 0.5
    A[np.tril_indices(1000, k=-1)] = values
    A = A + A.T
    return A

weighted_mean = lambda x1, x2, w1, w2: 1 / (w1 + w2) * (w1 * x1 + w2 * x2)

def average_ld(R1, R2, w1, w2):
    """
    weighted entrywise average of 2 ld matrices
    """
    return weighted_mean(R1, R2, w1, w2)

def fisher_average_ld(R1, R2, w1, w2):
    """
    weighted average of fisher-z-transformed correlations
    """
    tril = np.tril_indices(R1.shape[0], k=-1)
    z1 = np.arctanh(np.clip(R1[tril], -1+1e-5, 1-1e-5))
    z2 = np.arctanh(np.clip(R2[tril], -1+1e-5, 1-1e-5))
    A = np.tanh(weighted_mean(z1, z2, w1, w2))
    return tril2symm(A)

def data2ld(data, population='GTEx', ld='sample', **kwargs):
    superpop2samples = pickle.load(open(
        '/work-zfs/abattle4/karl/cosie_analysis/output/superpop2samples_1kG', 'rb'))

    if population == 'GTEx':
        X = data.data.X
    elif population != '1kG':
        X = data.data.genotype_1kG.loc[superpop2samples[population]].values
    else:
        X = data.data.genotype_1kG.values
    R_sample = np.corrcoef(X.T)
    w_sample = X.shape[0]
    
    R_z = z_ld(data)
    w_z = data.summary_stats.B.shape[0]
    
    if ld == 'sample':
        return R_sample
    elif ld == 'meanz':
        return average_ld(R_sample, R_z, w_sample, w_z)
    elif ld == 'fisherz':
        return fisher_average_ld(R_sample, R_z, w_sample, w_z)
    else:
        raise NotImplementedError

def smooth_betas(data, ld='sample', epsilon=0.0, **kwargs):
    """
    return a copy of data with smoothed effect sizes
    beta_sooth = SRS(SRS + epsilonS^2)^{-1} beta
    """
    if np.isclose(epsilon, 0):
        return data
    Bs = []
    R = ld_functions[ld](data)
    for i in range(data.summary_stats.S.shape[0]):
        S = np.diag(data.summary_stats.S.iloc[i].values)
        B = data.summary_stats.B.iloc[i].values
        SRS = S @ R @ S
        Bs.append(SRS @ np.linalg.solve(SRS + epsilon * S**2, B))
    data_smooth = deepcopy(data)
    data_smooth.summary_stats.B = pd.DataFrame(np.stack(Bs), columns=data.summary_stats.B.columns)
    return data_smooth

def init_gss(sim, update_variance=False, K=10, pi0=1.0, **kwargs):
    # TODO epsilon--- smooth expression?
    init_args = {
        'X': sim.data.X.T,
        'Y': sim.simulation.expression.values,
        'K': K,
        'snp_ids': sim.summary_stats.B.columns.values,
        'tissue_ids': sim.summary_stats.B.index.values
    }
    fit_args = {
        'update_weights': True,
        'update_pi': True,
        'ARD_weights': True,
        'update_variance': update_variance,
        'verbose': False,
        'max_iter': 50
    }
    name = kwargs['model_key']
    print('initializing summary stat model')
    gss = GSS(**init_args)
    gss.prior_activity = np.ones(K) * pi0
    gss.name = name
    return gss, fit_args

def fit_gss(sim, model_spec):
    print('fitting model')
    gss, fit_args = init_gss(sim, **model_spec.dropna().to_dict())
    gss.fit(**fit_args, update_active=False)
    gss.fit(**fit_args, update_active=True)
    compute_records_gss(gss)
    print('saving model to {}'.format(save_path))
    strip_and_dump(gss, save_path)

def init_css(sim, K=10, ld='sample', pi0=1.0, dispersion=1.0, **kwargs):
    # TODO epsilon--- smooth expression?
    init_args = {
        'LD': ld_functions[ld](sim),
        'B': sim.summary_stats.B.values,
        'S': sim.summary_stats.S.values,
        'K': K,
        'snp_ids': sim.summary_stats.B.columns.values,
        'tissue_ids': sim.summary_stats.B.index.values
    }
    fit_args = {
        'update_weights': True,
        'update_pi': True,
        'ARD_weights': True,
        'update_variance': False,
        'verbose': True,
        'max_iter': 50
    }
    name = kwargs['model_key']
    print('initializing summary stat model')
    css = CSS(**init_args)
    css.prior_activity = np.ones(K) * pi0
    css.tissue_precision_b = np.ones(css.dims['T']) * dispersion
    css.name = name
    return css, fit_args

def fit_css(data, model_spec):
    css, fit_args = init_css(data, **model_spec.to_dict())
    print('fitting model')
    css.fit(**fit_args, update_active=False)
    css.fit(**fit_args, update_active=True)
    compute_records_css(css)
    print('saving model to {}'.format(save_path))
    strip_and_dump(css, save_path)


def init_css2(sim, K=10, ld='sample', pi0=1.0, dispersion=1.0, **kwargs):
    # TODO epsilon--- smooth expression?
    init_args = {
        'LD': data2ld(sim, **kwargs),
        'B': sim.summary_stats.B.values,
        'S': sim.summary_stats.S.values,
        'K': K,
        'snp_ids': sim.summary_stats.B.columns.values,
        'tissue_ids': sim.summary_stats.B.index.values
    }
    fit_args = {
        'update_weights': True,
        'update_pi': True,
        'ARD_weights': True,
        'update_variance': False,
        'verbose': False,
        'max_iter': 50
    }
    name = kwargs['model_key']
    print('initializing summary stat model')
    css = CSS(**init_args)
    css.prior_activity = np.ones(K) * pi0
    css.tissue_precision_b = np.ones(css.dims['T']) * dispersion
    css.name = name
    return css, fit_args

def fit_css2(data, save_path, model_spec):
    css, fit_args = init_css2(data, **model_spec.to_dict())
    print('fitting model')
    css.fit(**fit_args, update_active=False)
    css.fit(**fit_args, update_active=True)
    compute_records_css(css)
    print('saving model to {}'.format(save_path))
    strip_and_dump(css, save_path)


sim_spec = pd.read_csv(snakemake.input[0], sep='\t')
model_spec = pd.read_csv(snakemake.input[1], sep='\t')
spec = sim_spec[sim_spec.sim_id == snakemake.wildcards.sim_id].iloc[0]
sim_data = load_sim_data(spec)

bp = '/'.join(snakemake.output[0].split('/')[:-1])
for _, ms in model_spec.iterrows():
    name = ms.model_key
    save_path = bp + '/' + name
    if not os.path.isfile(save_path):
        smoothed_data = smooth_betas(sim_data, **ms.dropna().to_dict())
        if ms.model == 'gss':
            fit_gss(smoothed_data, ms.dropna())
        if ms.model == 'css':
            fit_css2(smoothed_data, save_path, ms.dropna())
    else:
        print('{} alread fit'.format(name))

"""
bp = '/'.join(snakemake.output[0].split('/')[:-1])
for _, ms in model_spec.iterrows():
    name = ms.model_key
    save_path = bp + '/' + name
    if not os.path.isfile(save_path):
        smoothed_data = smooth_betas(sim_data, **ms.dropna().to_dict())
        if ms.model == 'gss':
            fit_gss(smoothed_data, ms.dropna())
        if ms.model == 'css':
            fit_css(smoothed_data, ms.dropna())
    else:
        print('{} alread fit'.format(name))
"""

