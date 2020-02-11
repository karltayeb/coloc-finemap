import pickle
import numpy as np
import pandas as pd
from scipy.stats import norm
from coloc.ard_ser import MVNFactorSER
from coloc.independent_model import IndependentFactorSER


def assign(obj, dictionary):
    for key in dictionary.keys():
        obj.__dict__[key] = dictionary[key]


def load_models_and_data(data_path, genotype_model_path, summary_model_path):
    """
    load models and data
    """
    genotype_model = IndependentFactorSER(np.zeros((1, 1)), np.zeros((1, 1)), 1)
    summary_model = MVNFactorSER(np.zeros((1, 1)), np.zeros((1, 1)), 1)

    data = pickle.load(open(data_path, 'rb'))
    assign(genotype_model, pickle.load(open(genotype_model_path, 'rb')))
    assign(summary_model, pickle.load(open(summary_model_path, 'rb')))

    genotype_model.X = data['X']
    genotype_model.Y = data['Y']

    summary_model.X = data['LD']
    summary_model.Y = data['zscores']
    return genotype_model, summary_model, data


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

def pair_coloc(df):
    pair_results = []
    num_tissues = np.unique(df.tissue).size
    for t1 in range(0, num_tissues):
        for t2 in range(t1, num_tissues):
            for k in np.unique(df.component):
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
    return pd.DataFrame(pair_results)

g, s, data = load_models_and_data(
    snakemake.input.data_path, snakemake.input.genotype_model_path, snakemake.input.summary_model_path)
df = make_table(g, data)
pairs = pair_coloc(df.loc[df.active==1])
pairs.to_csv(snakemake.output.genotype_output, index=False, sep='\t')

df = make_table(s, data)
pairs = pair_coloc(df.loc[df.active==1])
pairs.to_csv(snakemake.output.summary_output, index=False, sep='\t')