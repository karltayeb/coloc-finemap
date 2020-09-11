import numpy as np
import json
import glob
import pandas as pd
from collections import defaultdict
from os.path import isfile

GTEx_genotype_model = pd.read_csv('output/requests/GTEx_cafeh_genotype_ss.txt', header=None).iloc[:, 0].values
rule fit_gtex_cafeh_genotype_ss:
    input:
        expand('{path}', path=GTEx_genotype_model)

GTEx_genotype = pd.read_csv('output/requests/GTEx_get_gtex_genotype.txt', header=None).iloc[:, 0].values
rule get_all_gtex_genotype:
    input:
        expand('{path}', path=GTEx_genotype)

#GTEx_genotype_variant_report = pd.read_csv('output/requests/GTEx_cafeh_genotype_ss_variant_reports.txt', header=None).iloc[:, 0].values
paths = open('output/requests/GTEx_cafeh_genotype_ss_variant_reports.txt', 'r').read().split('\n')
GTEx_genotype_variant_report = [x for x in paths if not isfile(x)]
rule gtex_cafeh_genotype_variant_reports:
    input:
        expand('{path}', path=GTEx_genotype_variant_report)

cad_requests = pd.read_csv('output/CAD/requests.txt', header=None).iloc[:, 0].values
rule cad_gtex:
    input:
        expand('{path}', path=cad_requests)

UKBB_request = pd.read_csv('output/UKBB/individual_phenotype_requests.txt', header=None).iloc[:, 0].values
rule ukbb_gtex_individual:
    input:
        expand('{path}', path=UKBB_request)

def get_paths(request):
    """
    get list of tile from request file
    """
    paths = open(request, 'r').read().strip().split('\n')
    paths = [x for x in paths if not isfile(x)]
    return paths

rule terminal_rule:
    input:
        expand('{path}', path=get_paths(config['request']))

"""
BOGgenes = pd.read_csv(
    'output/GTEx/BOGgenes.txt', sep='\t', index_col=None)
get_path = lambda row: 'output/GTEx/{}/{}/{}.associations'.format(row.chromosome, row.gene, row.gene)
association_paths = [get_path(row) for _, row in BOGgenes.iterrows()]
rule get_BOG_associations:
    input:
        association_paths

css_1kG_paths = np.loadtxt('output/GTEx/gss_css_1kG.txt', dtype=str)
rule gss_css_1kG_gtex:
    input:
        expand('{path}', path=css_1kG_paths[:200])

fixedvar_paths = [x[:-7] + 'fixedvar.gss' for x in css_1kG_paths]
rule fit_gss_fixedvar:
    input:
        expand('{path}', path=fixedvar_paths[:200])

genes = pd.read_csv('output/sim/ld/1kgenes.txt', sep='\t', header=None).values[:, 0]
rule get_genotypes_for_sim:
	input:
		expand('{path}.raw', path=genes),
		expand('{path}.1kG.raw', path=genes),
		expand('{path}.snp2rsid.json', path=genes),
		expand('{path}.associations', path=genes)


gss_1kG_paths = np.loadtxt('output/requests/1kgenes.k20.pi01.gss.txt', dtype=str)
rule gss_1kG_random_k20_pi01:
	input:
		expand('{path}', path=gss_1kG_paths)

cad_known_coloc = np.loadtxt('output/requests/CAD_known_coloc.txt', dtype=str)
rule cad_known_coloc:
	input:
		expand('{path}', path=cad_known_coloc)

cad_n_egenes = np.loadtxt('output/requests/CAD_intersecting_egene.txt', dtype=str)
rule cad_n_egenes:
	input:
		expand('{path}', path=cad_n_egenes)
"""