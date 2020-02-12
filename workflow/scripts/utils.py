import pickle
import itertools
import numpy as np
import pandas as pd
from scipy.stats import norm
from coloc.ard_ser import MVNFactorSER
from coloc.independent_model import IndependentFactorSER

def compute_sigma2(X, true_effect, pve):
    var = np.var(true_effect @ X)
    sigma2_t = var/pve - var
    if sigma2_t == 0:
        # if variance is 0, there were no causal variants-- dont care what the variance is
        sigma2_t = 1.0
    return sigma2_t


def get_cis_variants(genotype, tss, window=500000, size=1000):
    pos = np.array([int(x.split('_')[1]) for x in genotype.index.values])
    #cis_variants = genotype.iloc[np.abs(pos - tss) < window]
    cis_variants = genotype.iloc[
        np.arange(genotype.shape[0])[np.argsort(np.abs(pos - tss))[:size]]
    ]
    return cis_variants

    
def assign(obj, dictionary):
    for key in dictionary.keys():
        obj.__dict__[key] = dictionary[key]


def matched_labels(model, data):
    """
    match components to true causal effects
    true_labels: [T, K] boolean indicatin if tissue t should use matched component k
    pse [T, K]: probabaility of a sign error from posterior distribtuion
    active [K]: does each component pass purity threshold
    matched [K]: is each component matched to a true causal variant (in 95% cs)
    num_tissues: the number of tissues that should be in a matched component
    """
    credible_sets, purity = model.get_credible_sets()
    active_sets = np.array([k for k in range(model.dims['K']) if purity[k] > 0.5])

    true_labels = []
    matched = []
    num_tissues = []
    for k in range(model.dims['K']):
        if k in active_sets:
            causal_snp = data['causal_snps'][np.isin(data['causal_snps'], credible_sets[k])]
            true_labels.append(np.any(data['true_effects'][:, causal_snp] != 0, 1))
            if causal_snp.size > 0:
                matched.append(k)
        else:
            true_labels.append(np.zeros(model.dims['T']).astype(bool))
        num_tissues.append(true_labels[-1].sum())

    true_labels = np.array(true_labels)
    pse = (norm.cdf(-np.abs(model.weight_means / np.sqrt(model.weight_vars))) * model.pi[None]).sum(-1)

    active = np.isin(np.arange(model.dims['K']), active_sets)
    matched = np.isin(np.arange(model.dims['K']), np.array(matched))
    return true_labels, pse, active, matched, np.array(num_tissues)


def make_table(model, data):
    true_labels, pse, a, m, n = matched_labels(model, data)
    df = pd.concat([
        pd.DataFrame(pse).T,
        pd.DataFrame(np.stack([a, m, n]).T, columns=['active', 'matched', 'num_tissues'])], axis=1)

    df.index.name = 'component'

    df = df.reset_index()
    df = df.melt(id_vars=['active', 'matched', 'num_tissues', 'component'], value_name='p_sign_error')
    df = df.rename({'variable':'tissue'}, axis=1)
    df.loc[:, 'label'] = true_labels.T.flatten()
    df.loc[:, 'ard_variance'] = 1 / model.prior_precision.flatten()
    return df

def pair_coloc(df, data):
    tissue_colocalize = data['true_effects'] @ data['true_effects'].T
    pair_results = []
    num_tissues = np.unique(df.tissue).size
    for t1, t2 in itertools.combinations(np.arange(num_tissues), 2):
        for k in np.unique(df.loc[df.active==1].component):
            d = {
                't1': t1,
                't2': t2,
                't1_p_sign_error': df.loc[(df.tissue == t1) & (df.component == k)].p_sign_error.iloc[0],
                't2_p_sign_error': df.loc[(df.tissue == t2) & (df.component == k)].p_sign_error.iloc[0],
                't1_ard_variance': df.loc[(df.tissue == t1) & (df.component == k)].ard_variance.iloc[0],
                't2_ard_variance': df.loc[(df.tissue == t2) & (df.component == k)].ard_variance.iloc[0],
                'k': k,
                'matched': np.any(df.loc[df.component == k].matched == 1),
                'label': df.loc[(df.tissue == t1) & (df.component == k)].label.iloc[0] * \
                    df.loc[(df.tissue == t2) & (df.component == k)].label.iloc[0]
            }
            pair_results.append(d)
        if np.unique(df.loc[df.active==1].component).size == 0:
            d = {
                't1': t1,
                't2': t2,
                't1_p_sign_error': 0.5,
                't2_p_sign_error': 0.5,
                't1_ard_variance': 0.0,
                't2_ard_variance': 0.0,
                'k': -1,
                'matched': False,
                'label': bool(tissue_colocalize[t1, t2])
            }
            pair_results.append(d)
    return pd.DataFrame(pair_results)