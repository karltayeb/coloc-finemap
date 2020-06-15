import pandas as pd
import pickle
import numpy as np
from coloc.independent_model_ss import IndependentFactorSER as GSS
from coloc.cafeh_ss import CAFEH as CSS

from coloc.misc import *
from coloc.simulation import *
from coloc.covariance import *

def init_css(data, ld_type, K=10, p=1.0, d=1.0):
    init_args = {
        'LD': ld_functions[ld_type](data),
        'B': data.B.values,
        'S': data.S.values,
        'K': K,
        'snp_ids': data.B.columns.values,
        'tissue_ids': data.B.index.values
    }
    name = 'sim-{}_gene-{}_k-{}_pi-{}_d-{}.ld-{}.css'.format(
        data.id, data.gene, K, p, d, ld_type)
    print('initializing summary stat model')
    css = CSS(**init_args)
    css.prior_activity = np.ones(K) * p
    css.tissue_precision_b = np.ones(css.dims['T']) * d
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
    'lw_sample': ledoit_wolf_sample_ld,
    'lw_refence': ledoit_wolf_reference_ld,
    'lw_z': ledoit_wolf_z_ld,
    'ref_z': ref_z_ld,
    'z3': z3_ld
}

#gss = load(snakemake.input[0])
gene = snakemake.wildcards.gene
gss = load('output/GTEx/{}/{}/{}.k20.pi01.gss'.format(
    get_chromosome(gene), gene, gene
))
sim_data = sim_from_model(gene, gss, thin=True, sim_id=snakemake.wildcards.simid)

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


# fit CSS to simulation data
fit_args = {
'update_weights': True,
'update_pi': True,
'ARD_weights': True,
'update_variance': False,
'verbose': False,
'max_iter': 50
}

bp = '/'.join(snakemake.output[0].split('/')[:-1])
for ld_type in ld_functions:
    print(ld_type)
    css = init_css(sim_data, ld_type, 20, 0.1, 1.0)
    css.fit(**fit_args, update_active=False)
    css.fit(**fit_args, update_active=True)

    compute_records_css(css)
    save_path = bp + '/' + css.name
    print('saving model to {}'.format(save_path))
    strip_and_dump(css, save_path)
    rehydrate_model(css)
