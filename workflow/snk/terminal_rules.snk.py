import numpy as np
import json
import glob
import pandas as pd
from collections import defaultdict

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
