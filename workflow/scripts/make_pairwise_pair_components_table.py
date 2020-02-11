import pickle
import numpy as np
import pandas as pd
from scipy.stats import norm
from coloc.ard_ser import MVNFactorSER
from coloc.independent_model import IndependentFactorSER
from make_tissue_pair_components_table import assign, load_models_and_data, matched_labels, make_table, pair_coloc

def load_data(data_path):
    """
    load models and data
    """
    return pickle.load(open(data_path, 'rb'))


def load_model(data, genotype_model_path=None, summary_model_path=None):
    path = genotype_model_path or summary_model_path
    t1 = int(path.split('/')[-1].split('_')[1])
    t2 = int(path.split('/')[-1].split('_')[3])

    if genotype_model_path is not None:
        model = IndependentFactorSER(np.zeros((1, 1)), np.zeros((1, 1)), 1)
        assign(model, pickle.load(open(genotype_model_path, 'rb')))
        model.X = data['X']
        #model.Y = data['Y']# [[t1, t2]]
    else:
        model = MVNFactorSER(np.zeros((1, 1)), np.zeros((1, 1)), 1)
        assign(model, pickle.load(open(summary_model_path, 'rb')))
        model.X = data['X']
        #model.Y = data['zscores']# [[t1, t2]]
    return model, t1, t2

if __name__ == '__main__':
    summary_output = 'output/simulation/single_causal_variant/pve_0.2/ld_0.8/gene_ENSG00000262000.1/pairwise_summary/pairs_summary'
    data_paths = ['output/simulation/single_causal_variant/pve_0.2/ld_0.8/gene_ENSG00000262000.1/data']
    summary_model_paths = ['output/simulation/single_causal_variant/pve_0.2/ld_0.8/gene_ENSG00000262000.1/pairwise_summary/t1_0_t2_2_model_summary', 
        'output/simulation/single_causal_variant/pve_0.2/ld_0.8/gene_ENSG00000262000.1/pairwise_summary/t1_0_t2_5_model_summary', 
        'output/simulation/single_causal_variant/pve_0.2/ld_0.8/gene_ENSG00000262000.1/pairwise_summary/t1_1_t2_2_model_summary',
        'output/simulation/single_causal_variant/pve_0.2/ld_0.8/gene_ENSG00000262000.1/pairwise_summary/t1_1_t2_5_model_summary',
        'output/simulation/single_causal_variant/pve_0.2/ld_0.8/gene_ENSG00000262000.1/pairwise_summary/t1_4_t2_2_model_summary',
        'output/simulation/single_causal_variant/pve_0.2/ld_0.8/gene_ENSG00000262000.1/pairwise_summary/t1_4_t2_5_model_summary']

    summary_keys = []
    for data_path in data_paths:
        data = load_data(data_path)
        key = '/'.join(data_path.split('/')[5:-1])

        sub_summary_paths = [x for x in summary_model_paths if key in x]
        for summary_model_path in sub_summary_paths:
            model, t1, t2 = load_model(data, summary_model_path=summary_model_path)
            df = make_table(model, data)

            import pdb; pdb.set_trace()
            #df = df.loc[df.tissue in [t1, t2]]
            pairs = pair_coloc(df.loc[df.active == 1])
            if pairs.size > 0:
                summary_pairs.append(pairs)
        summary_pairs = pd.concat(summary_pairs)
        summary_pairs.to_csv(summary_output, index=False, sep='/t')
