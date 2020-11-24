import pickle
import numpy as np
import pandas as pd
from cafeh.cafeh_ss import CAFEH
from cafeh.misc import plot_components

import sys
sys.path.append('/work-zfs/abattle4/karl/cosie_analysis/utils/')
from misc import *

def _get_minalpha(pi):
    """
    report the minimum alpha value to include this snp in cs
    """
    argsort = np.flip(np.argsort(pi))
    resort = np.argsort(argsort)
    cumsum = np.cumsum(pi[argsort])
    minalpha = np.roll(cumsum, 1)
    minalpha[0] = 0
    return minalpha[resort]

def get_minalpha(model):
    return  pd.DataFrame(
        np.array([_get_minalpha(model.pi[k]) for k in range(model.dims['K'])]).T,
        index=model.snp_ids
    )

# load a model
gene = snakemake.wildcards.gene
model = pickle.load(open(snakemake.input.model, 'rb'))
model._decompress_model()

study_pip = model.get_study_pip().T
table = study_pip.reset_index().melt(id_vars='index').rename(columns={
    'index': 'variant_id',
    'variable': 'study',
    'value': 'pip' 
})

v2r = load_var2rsid(gene)
table.loc[:, 'rsid'] = table.variant_id.apply(lambda x: v2r.get(x, '-'))

top_component = pd.Series(model.pi.argmax(0), index=model.snp_ids).to_dict()
table.loc[:, 'top_component'] = table.variant_id.apply(lambda x: top_component.get(x))

minalpha = get_minalpha(model).to_dict()
table.loc[:, 'alpha'] = [minalpha.get(k).get(v) for k, v in zip(table.top_component.values, table.variant_id.values)]

rank = pd.DataFrame({k: np.argsort(np.flip(np.argsort(model.pi[k]))) for k in range(model.dims['K'])}, index=model.snp_ids).to_dict()
table.loc[:, 'rank'] = [rank.get(k).get(v) for k, v in zip(table.top_component.values, table.variant_id.values)]

active = pd.DataFrame(model.active, index=model.study_ids)
active.loc['all'] = (model.active.max(0) > 0.5).astype(int)
active = active.to_dict()
table.loc[:, 'p_active'] = [active.get(k).get(s) for k, s in zip(table.top_component.values, table.study.values)]

pi = pd.Series(model.pi.max(0), index=model.snp_ids).to_dict()
table.loc[:, 'pi'] = table.variant_id.apply(lambda x: pi.get(x))

table.loc[:, 'chr'] = get_chr(gene)
table.loc[:, 'start'] = table.variant_id.apply(lambda x: int(x.split('_')[1]))
table.loc[:, 'end'] = table.start + 1
table.loc[:, 'gene'] = gene

table.loc[:, ['chr', 'start', 'end', 'variant_id', 'rsid', 'study', 'gene', 'pip', 'pi', 'top_component', 'p_active', 'alpha']]

table = table.loc[:, ['chr', 'start', 'end', 'variant_id', 'rsid', 'study', 'pip', 'top_component', 'p_active', 'pi', 'alpha', 'rank']]
small_table = table[table.p_active > 0.5].sort_values(by=['chr', 'start'])

# add effect size and variance
study2idx = {s: i for i, s in enumerate(model.study_ids)}
var2idx = {s: i for i, s in enumerate(model.snp_ids)}
small_table.loc[:, 'effect'] = [model.weight_means[study2idx.get(s), c, var2idx.get(v)] for s, v, c in zip(
    small_table.study, small_table.variant_id, small_table.top_component)]
small_table.loc[:, 'effect_var'] = [model.weight_vars[study2idx.get(s), c, var2idx.get(v)] for s, v, c in zip(
    small_table.study, small_table.variant_id, small_table.top_component)]

small_table.to_csv(snakemake.output.report, sep='\t', index=None)

