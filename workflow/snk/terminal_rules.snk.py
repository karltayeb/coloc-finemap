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
rule gss_css_1kG_gtex:
    input:
        expand('{path}', path=fixedvar_paths[:200])