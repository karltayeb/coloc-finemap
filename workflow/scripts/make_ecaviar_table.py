import pickle
import itertools
import pandas as pd
import numpy as np
import glob

dfs = []
key = []
for p in snakemake.input.caviar_posteriors:
    dfs.append(pd.read_csv(p, sep='\t', usecols=[2]))
    key.append(p.split('_')[-2])

df = pd.concat(dfs, keys=key, axis=1)

table = []
for t1, t2 in itertools.combinations(np.arange(len(key)), 2):
    clpp = (df.loc[:, 't{}'.format(t1)] * df.loc[:, 't{}'.format(t2)]).max()[0]
    table.append({
        't1': t1,
        't2': t2,
        'CLPP': clpp
    })
table = pd.DataFrame(table)
table.to_csv(snakemake.output[0])