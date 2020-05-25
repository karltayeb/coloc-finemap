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