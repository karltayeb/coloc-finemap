import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
import json
from cafeh.kls import categorical_kl
from types import SimpleNamespace
import os
import random
import string
from copy import deepcopy

from types import SimpleNamespace
import subprocess
import copy
import json
import os
import numpy as np
import pandas as pd

gc26 = pd.read_csv(
    '/work-zfs/abattle4/lab_data/annotation/gencode.v26/gencode.v26.annotation.gene.txt', sep='\t')

def get_chr(gene):
    return gc26.set_index('gene_id').chr.to_dict().get(gene, None)


def get_tss(gene):
    return gc26.set_index('gene_id').start_pos.to_dict().get(gene, None)


def plink_get_genotype(gene, bfile, save_path):
    tss = get_tss(gene)
    cmd = ' '.join(
        ['plink',
         '--bfile', bfile,
         '--chr', get_chr(gene)[3:],
         '--from-bp', str(int(np.maximum(tss-1e6, 0))),
         '--to-bp', str(int(tss+1e6)),
         '--maf', '0.01',
         '--geno', '0.1',
         '--recode', 'A',
         '--keep-allele-order',
         '--snps-only', '--write-snplist', '--allow-no-sex',
         '--out', save_path])
    return cmd


def load_var2rsid(gene):
    v2rp = '/work-zfs/abattle4/karl/cosie_analysis/output/GTEx/{}/{}/{}.snp2rsid'.format(
        get_chr(gene), gene, gene)
    v2r = json.load(open(v2rp, 'r'))
    return v2r


def load_gtex_genotype(gene, use_rsid=False):
    """
    load gtex genotype for variants in 1Mb window of gene tss
    @param gene: ensemble gene id
    @param use_rsid: boolean to convert variant id to rsid
    """
    gp = '/work-zfs/abattle4/karl/cosie_analysis/output/GTEx/{}/{}/{}.raw'.format(
        get_chr(gene), gene, gene)
    v2r = load_var2rsid(gene)

    print('loading gtex genotypes...')
    genotype = pd.read_csv(gp, sep=' ')
    genotype = genotype.set_index('IID').iloc[:, 5:]

    # recode genotypes
    coded_snp_ids = np.array([x.strip() for x in genotype.columns])
    snp_ids = {x: '_'.join(x.strip().split('_')[:-1]) for x in coded_snp_ids}
    genotype.rename(columns=snp_ids, inplace=True)

    if use_rsid:
        genotype.rename(columns=v2r, inplace=True)
    return genotype


def load_gtex_associations(gene):
    """
    ap = '/work-zfs/abattle4/karl/cosie_analysis/output/GTEx/{}/{}/{}.associations'.format(
        get_chr(gene), gene, gene)
    v2rp = '/work-zfs/abattle4/karl/cosie_analysis/output/GTEx/{}/{}/{}.snp2rsid.json'.format(
        get_chr(gene), gene, gene)
    """
    v2r = load_var2rsid(gene)
    ap = '/work-zfs/abattle4/karl/cosie_analysis/output/GTEx/{}/{}/{}.associations'.format(
        get_chr(gene), gene, gene)

    associations = pd.read_csv(ap, index_col=0)
    associations.loc[:, 'rsid'] = associations.variant_id.apply(lambda x: v2r.get(x, '-'))
    associations.loc[:, 'pos'] = associations.variant_id.apply(lambda x: int(x.split('_')[1]))
    associations.loc[:, 'ref'] = associations.variant_id.apply(lambda x: x.split('_')[2])
    associations.loc[:, 'alt'] = associations.variant_id.apply(lambda x: x.split('_')[3])
    associations.loc[:, 'sample_size'] = (associations.ma_count / associations.maf / 2)
    associations.loc[:, 'S'] = np.sqrt(
        associations.slope**2/associations.sample_size + associations.slope_se**2)
    return associations


def load_gtex_expression(gene):
    """
    load expression, drop unexpressed individuals
    """
    # load expression
    ep = '/work-zfs/abattle4/karl/cosie_analysis/output/GTEx/{}/{}/{}.expression'.format(get_chr(gene), gene, gene)
    gene_expression = pd.read_csv(ep, sep='\t', index_col=0)
    # drop individuals that do not have recorded expression
    gene_expression = gene_expression.loc[:, ~np.all(np.isnan(gene_expression), 0)]
    return gene_expression


def make_snp_format_table(gene):
    ep = '/work-zfs/abattle4/karl/cosie_analysis/output/GTEx/{}/{}/{}.expression'.format(get_chr(gene), gene, gene)
    gp = '/work-zfs/abattle4/karl/cosie_analysis/output/GTEx/{}/{}/{}.raw'.format(get_chr(gene), gene, gene)
    gp1kG = '/work-zfs/abattle4/karl/cosie_analysis/output/GTEx/{}/{}/{}.1kG.raw'.format(get_chr(gene), gene, gene)
    ap = '/work-zfs/abattle4/karl/cosie_analysis/output/GTEx/{}/{}/{}.associations'.format(get_chr(gene), gene, gene)
    v2r = load_var2rsid(gene)


    with open(gp, 'r') as f:
        snps = f.readline().strip().split()[6:]

    with open(gp1kG, 'r') as f:
        rsids = f.readline().strip().split()[6:]

    summary_stat_snps = np.unique(pd.read_csv(ap, sep=',', usecols=[3]).variant_id)
    summary_stat_snps = {snp: True for snp in summary_stat_snps}
    vid_codes = {'_'.join(x.split('_')[:-1]): x.split('_')[-1] for x in snps}
    rsid_codes = {x.split('_')[0]: x.split('_')[1] for x in rsids}
    table = []
    for vid in vid_codes:
        ref = vid.split('_')[-2]
        rsid = v2r.get(vid, '-')
        table.append({
            'variant_id': vid,
            'rsid': v2r.get(vid, '-'),
            'ref': ref,
            'flip_gtex': ref != vid_codes.get(vid, '-'),
            'flip_1kG': ref != rsid_codes.get(rsid, '-'),
            'in_1kG': rsid_codes.get(rsid, False) != False,
            'has_test': summary_stat_snps.get(vid, False)
        })
    return pd.DataFrame(table)




# DEPRECATE
def get_common_snps(gene):
    table = make_snp_format_table(gene)
    return table[table.has_test & table.in_1kG].rsid.values


def load_old(model_path):
    model = pickle.load(open(model_path, 'rb'))
    rehydrate_model(model)
    model.name = model_path.split('/')[-1]
    return model

def load(path):
    model = pickle.load(open(path, 'rb'))
    model._decompress_model()
    return model


def need_to_flip(variant_id):
    _, _, major, minor, _, ref = variant_id.strip().split('_')
    if minor != ref:
        return True
    else:
        return False


flip = lambda x: (x-1)*-1 + 1

def load_1kG_genotype(gene):
    gp1kG = '/work-zfs/abattle4/karl/cosie_analysis/output/GTEx/{}/{}/{}.1kG.raw'.format(get_chr(gene), gene, gene)
    table = make_snp_format_table(gene)

    genotype = pd.read_csv(gp1kG, sep=' ')
    genotype = genotype.set_index('IID').iloc[:, 5:]

    # recode genotypes
    coded_snp_ids = np.array([x.strip() for x in genotype.columns])
    snp_ids = {x: '_'.join(x.strip().split('_')[:-1]) for x in coded_snp_ids}
    genotype.rename(columns=snp_ids, inplace=True)

    flip_1kG = table[table.flip_1kG & (table.rsid != '-')].rsid.values
    flip_1kG = np.intersect1d(flip_1kG, genotype.columns)
    genotype.loc[:, flip_1kG] = genotype.loc[:, flip_1kG].applymap(flip)
    return genotype


def load_gtex_summary_stats(gene):
    ap = '/work-zfs/abattle4/karl/cosie_analysis/output/GTEx/{}/{}/{}.associations'.format(get_chr(gene), gene, gene)
    v2r = load_var2rsid(gene)

    associations = pd.read_csv(ap)
    associations.loc[:, 'sample_size'] = (associations.ma_count / associations.maf / 2)
    Ba = associations.pivot('study', 'variant_id', 'slope')
    Va = associations.pivot('study', 'variant_id', 'slope_se')**2
    n = associations.pivot('study', 'variant_id', 'sample_size')
    Sa = np.sqrt(Ba**2/n + Va)

    [x.rename(columns=v2r, inplace=True) for x in [Ba, Sa, Va, n]];
    return Ba, Sa, Va, n


def load_genotype(genotype_path, flip=False):
    """
    fetch genotype
    flip codings that are need to be flipped
    set snp_ids to be consistent with gtex
    """
    #load genotype
    genotype = pd.read_csv(genotype_path, sep=' ')
    genotype = genotype.set_index('IID').iloc[:, 5:]

    # recode genotypes
    coded_snp_ids = np.array([x.strip() for x in genotype.columns])
    snp_ids = {x: '_'.join(x.strip().split('_')[:-1]) for x in coded_snp_ids}
    ref = {'_'.join(x.strip().split('_')[:-1]): x.strip().split('_')[-1] for x in coded_snp_ids}
    genotype.rename(columns=snp_ids, inplace=True)

    if flip:
        flips = np.array([need_to_flip(vid) for vid in coded_snp_ids])
        genotype.iloc[:, flips] = genotype.iloc[:, flips].applymap(flip) 
    return genotype, ref


def load_gtex_residual_expression(gene):
    """
    load residual expression as dataframe
    """
    expression = load_gtex_expression(gene)
    covariates = pd.read_csv(
        '/work-zfs/abattle4/karl/cosie_analysis/output/GTEx/covariates.csv',
        sep='\t', index_col=[0, 1])

    residual_expression = {}
    for study in expression.index.values:
        samples = expression.columns.values[~np.isnan(expression.loc[study])]
        beta = np.linalg.pinv(covariates.loc[study, samples].T) @ expression.loc[study, samples]
        residual_expression[study] = \
            expression.loc[study, samples] - covariates.loc[study, samples].T @ beta

    residual_expression = pd.DataFrame(residual_expression).T
    return residual_expression.loc[expression.index].loc[:, expression.columns]


def center_mean_impute(genotype):
    """
    center columns of dataframe
    fill na with 0 (mean imputation)
    """
    X = genotype - genotype.mean()
    X = X.fillna(0)
    return X

def load_gene_data(gene, thin=False):
    # Load GTEx and 1kG genotype
    # flip genotype encoding to be consistent with GTEx associations
    print('loading genotypes...')
    genotype = load_gtex_genotype(gene)
    genotype1kG = load_1kG_genotype(gene)

    print('loading expression...')
    expression = load_gtex_expression(gene)

    # load GTEx summary stats
    print('loading associations...')
    B, S, V, n = load_gtex_summary_stats(gene)

    # filter down to list of snps present in GTEx and 1kG
    print('filtering down to common snps')
    common_snps = get_common_snps(gene)

    if thin:
        d = int(common_snps.size / 1000)
        common_snps = common_snps[::d]
        common_snps = common_snps[:1000]

    genotype = genotype.loc[:, common_snps]
    genotype1kG = genotype1kG.loc[:, common_snps]
    B = B.loc[:, common_snps]
    S = S.loc[:, common_snps]
    V = V.loc[:, common_snps]
    n = n.loc[:, common_snps]

    X = center_mean_impute(genotype).values
    X1kG = center_mean_impute(genotype1kG).values

    covariates = pd.read_csv(
        '/work-zfs/abattle4/karl/cosie_analysis/output/GTEx/covariates.csv', sep='\t', index_col=[0, 1])
    covariates = covariates.loc[expression.index].loc[:, genotype.index.values]

    return SimpleNamespace(**{
        'genotype_1kG': genotype1kG.loc[:, common_snps],
        'genotype_gtex': genotype.loc[:, common_snps],
        'X': X, 'X1kG': X1kG, 'expression': expression,
        'B': B, 'S': S, 'V':V, 'n':n, 'common_snps': common_snps,
        'gene': gene, 'id': gene, 'covariates': covariates
    })


def compute_sigma2(X, true_effect, pve):
    var = np.var(true_effect @ X.T)
    sigma2_t = var/pve - var
    if sigma2_t == 0:
        # if variance is 0, there were no causal variants-- dont care what the variance is
        sigma2_t = 1.0
    return sigma2_t


def linregress(y, X):
    diag = np.einsum('ij,ij->i', X.T, X.T)
    betas = y @ X / diag
    var = np.var(y[:, None] - betas * X, 0) / diag
    s2 = betas**2 / y.size + var
    return betas, var, s2


def make_gtex_genotype_data_dict(expression_path, genotype_path, standardize=False, flip=False):
    # load expression
    gene_expression = load_gtex_expression(expression_path)
    genotype, ref = load_genotype(genotype_path, flip)

    # center, mean immpute
    genotype = (genotype - genotype.mean(0))
    genotype = genotype.fillna(0)

    # standardize
    if standardize:
        genotype = genotype / genotype.std(0)

    # filter down to common individuals
    individuals = np.intersect1d(genotype.index.values, gene_expression.columns.values)
    genotype = genotype.loc[individuals]
    gene_expression = gene_expression.loc[:, individuals]


    # load covariates
    covariates = pd.read_csv('/work-zfs/abattle4/karl/cosie_analysis/output/GTEx/covariates.csv', sep='\t', index_col=[0, 1])
    covariates = covariates.loc[gene_expression.index]
    covariates = covariates.loc[:, genotype.index.values]
    X = genotype.values.T
    data = {
        'X': X,
        'Y': gene_expression.values,
        'covariates': covariates,
        'snp_ids': genotype.columns.values,
        'sample_ids': genotype.index.values,
        'study_ids': gene_expression.index.values
    }
    return data

def compute_summary_stats(data):
    B = {}
    V = {}
    S = {}
    for i, study in enumerate(data['study_ids']):
        cov = data['covariates'].loc[study]
        mask = ~np.isnan(cov.iloc[0])
        cov = cov.values[:, mask]
        y = data['Y'][i, mask]
        X = data['X'][:, mask]

        #H = cov.T @ np.linalg.solve(cov @ cov.T, cov)
        H = (np.linalg.pinv(cov) @ cov)
        yc = y - y @ H
        Xc = X - X @ H
        # prep css data
        B[study] = (Xc @ yc) / np.einsum('ij,ij->i', Xc, Xc)
        r = yc - B[study][:, None]*Xc
        V[study] = np.einsum('ij,ij->i', r, r) / np.einsum('ij,ij->i', Xc, Xc) / (yc.size)
        S[study] = np.sqrt(B[study]**2/yc.size + V[study])

    B = pd.DataFrame(B, index=data['snp_ids']).T
    V = pd.DataFrame(V, index=data['snp_ids']).T
    S = pd.DataFrame(S, index=data['snp_ids']).T
    return B, S, V


def rehydrate_model(model):
    model.weight_means = np.zeros((model.dims['T'],model.dims['K'],model.dims['N']))
    model.weight_vars = np.ones((model.dims['T'],model.dims['K'],model.dims['N']))

    model.weight_means[:, :, model.records['snp_subset']] = model.records['mini_wm']
    model.weight_vars[:, :, model.records['snp_subset']] = model.records['mini_wv']


def load_model(model_path, expression_path=None, genotype_path=None, load_data=False):
    if expression_path is None:
        gene = model_path.split('/')[-2]
        base_path = '/'.join(model_path.split('/')[:-1])
        expression_path = '{}/{}.expression'.format(base_path, gene)
    if genotype_path is None:
        gene = model_path.split('/')[-2]
        base_path = '/'.join(model_path.split('/')[:-1])
        genotype_path = '{}/{}.raw'.format(base_path, gene)

    model = pickle.load(open(model_path, 'rb'))
    rehydrate_model(model)

    if load_data:
        data = make_gtex_genotype_data_dict(expression_path, genotype_path)
        model.__dict__.update(data)
    return model


def compute_records(model):
    """
    save the model with data a weight parameters removed
    add 'mini_weight_measn' and 'mini_weight_vars' to model
    the model can be reconstituted in a small number of iterations
    """
    PIP = 1 - np.exp(np.log(1 - model.pi + 1e-10).sum(0))
    mask = (PIP > 0.01)
    wv = model.weight_vars[:, :, mask]
    wm = model.weight_means[:, :, mask]

    credible_sets, purity = model.get_credible_sets(0.99)
    active = np.array([purity[k] > 0.5 for k in range(model.dims['K'])])
    records = {
        'active': active,
        'purity': purity,
        'credible_sets': credible_sets,
        'EXz': model.pi @ model.X,
        'mini_wm': wm,
        'mini_wv': wv,
        'snp_subset': mask
    }
    model.records = records


def compute_records_gss(model):
    """
    save the model with data a weight parameters removed
    add 'mini_weight_measn' and 'mini_weight_vars' to model
    the model can be reconstituted in a small number of iterations
    """
    credible_sets, purity = model.get_credible_sets(0.999)
    active = model.active.max(0) > 0.5
    try:
        snps = np.unique(np.concatenate([
            credible_sets[k] for k in range(model.dims['K']) if active[k]]))
    except Exception:
        snps = np.unique(np.concatenate([
            credible_sets[k][:5] for k in range(model.dims['K'])]))
    mask = np.isin(model.snp_ids, snps)

    wv = model.weight_vars[:, :, mask]
    wm = model.weight_means[:, :, mask]

    records = {
        'active': active,
        'purity': purity,
        'credible_sets': credible_sets,
        'mini_wm': wm,
        'mini_wv': wv,
        'snp_subset': mask
    }
    model.records = records


def compute_records_css(model):
    """
    save the model with data a weight parameters removed
    add 'mini_weight_measn' and 'mini_weight_vars' to model
    the model can be reconstituted in a small number of iterations
    """
    credible_sets, purity = model.get_credible_sets(0.999)
    active = model.active.max(0) > 0.5
    try:
        snps = np.unique(np.concatenate([
            credible_sets[k] for k in range(model.dims['K']) if active[k]]))
    except Exception:
        snps = np.unique(np.concatenate([
            credible_sets[k][:5] for k in range(model.dims['K'])]))
    mask = np.isin(model.snp_ids, snps)

    wv = model.weight_vars[:, :, mask]
    wm = model.weight_means[:, :, mask]

    records = {
        'active': active,
        'purity': purity,
        'credible_sets': credible_sets,
        'mini_wm': wm,
        'mini_wv': wv,
        'snp_subset': mask
    }
    model.records = records


def strip_and_dump(model, path, save_data=False):
    """
    save the model with data a weight parameters removed
    add 'mini_weight_measn' and 'mini_weight_vars' to model
    the model can be reconstituted in a small number of iterations
    """
    # purge precompute
    for key in model.precompute:
        model.precompute[key] = {}
    model.__dict__.pop('weight_means', None)
    model.__dict__.pop('weight_vars', None)
    if not save_data:
        model.__dict__.pop('X', None)
        model.__dict__.pop('Y', None)
        model.__dict__.pop('covariates', None)
        model.__dict__.pop('LD', None)
    pickle.dump(model, open(path, 'wb'))


def compute_pip(model):
    active = model.records['active']
    return 1 - np.exp(np.log(1 - model.pi + 1e-10).sum(0))


def component_scores(model):
    purity = model.get_credible_sets(0.99)[1]
    active = np.array([purity[k] > 0.5 for k in range(model.dims['K'])])
    if active.sum() > 0:
        mw = model.weight_means
        mv = model.weight_vars
        pi = model.pi
        scores = np.einsum('ijk,jk->ij', mw / np.sqrt(mv), model.pi)
        weights = pd.DataFrame(
            scores[:, active],
            index = model.study_ids,
            columns = np.arange(model.dims['K'])[active]
        )
    else:
        weights = pd.DataFrame(
            np.zeros((model.dims['T'], 1)),
            index = model.study_ids
        )
    return weights


def make_variant_report(model, gene):
    """
    report variant level infromation from model
    """
    PIP = 1 - np.exp(np.log(1 - model.pi + 1e-10).sum(0))
    purity = model.get_credible_sets(0.99)[1]
    active = np.array([purity[k] > 0.5 for k in range(model.dims['K'])])
    if active.sum() == 0:
        active[0] = True

    pi = pd.DataFrame(model.pi.T, index=model.snp_ids)
    cset_alpha = pd.concat(
        [pi.iloc[:, k].sort_values(ascending=False).cumsum() - pi.iloc[:, k]
         for k in np.arange(model.dims['K']) if active[k]],
        sort=False, axis=1
    )

    most_likely_k = np.argmax(pi.values[:, active], axis=1)
    most_likely_p = np.max(pi.values[:, active], axis=1)
    most_likely_cset_alpha = cset_alpha.values[np.arange(pi.shape[0]), most_likely_k]

    A = pd.DataFrame(
        [PIP, most_likely_k, most_likely_p, most_likely_cset_alpha],
        index=['PIP','k', 'p', 'min_alpha'], columns=model.snp_ids).T

    A.loc[:, 'chr'] = [x.split('_')[0] for x in A.index]
    A.loc[:, 'start'] = [int(x.split('_')[1]) for x in A.index]
    A.loc[:, 'end'] = A.loc[:, 'start'] + 1
    A.reset_index(inplace=True)
    A.loc[:, 'variant_id'] = A.loc[:, 'index'].apply(lambda x: '_'.join(x.split('_')[:-1]))
    A.loc[:, 'ref'] = A.loc[:, 'index'].apply(lambda x: x.split('_')[-1])


    A.loc[:, 'gene_id'] = gene
    A = A.set_index(['chr', 'start', 'end'])
    return A.loc[:, ['gene_id', 'variant_id', 'ref', 'PIP', 'k', 'p', 'min_alpha']]


def plot_components(self, thresh=0.5, save_path=None, show=True):
    """
    plot inducing point posteriors, weight means, and probabilities
    """
    weights = self.get_expected_weights()
    active_components = self.active.max(0) > 0.5
    if not np.any(active_components):
        active_components[0] = True

    fig, ax = plt.subplots(1, 3, figsize=(18, 4))
    sns.heatmap(self.active[:, active_components], ax=ax[0],
                cmap='Blues', xticklabels=np.arange(active_components.size)[active_components])
    sns.heatmap(self.get_expected_weights()[:, active_components], ax=ax[1],
                cmap='RdBu', center=0, xticklabels=np.arange(active_components.size)[active_components])

    for k in np.arange(self.dims['K'])[active_components]:
        ax[2].scatter(
            np.arange(self.dims['N'])[self.pi.T[:, k] > 2/self.dims['N']],
            self.pi.T[:, k][self.pi.T[:, k] > 2/self.dims['N']],
            alpha=0.5, label='k{}'.format(k))
    ax[2].scatter(np.arange(self.dims['N']), np.zeros(self.dims['N']), alpha=0.0)
    ax[2].set_title('pi')
    ax[2].set_xlabel('SNP')
    ax[2].set_ylabel('probability')
    ax[2].legend(bbox_to_anchor=(1.04,1), loc="upper left")


def kl_components(m1, m2):
    """
    pairwise kl of components for 2 models
    """
    kls = np.array([[
        categorical_kl(m1.pi[k1], m2.pi[k2])
        + categorical_kl(m2.pi[k2], m1.pi[k1])
        for k1 in range(m1.dims['K'])] for k2 in range(m2.dims['K'])])
    return kls


def kl_heatmap(m1, m2, thresh=0.5):
    a1 = m1.active.max(0) > thresh
    if not np.any(a1):
        a1[0] = True
    a2 = m2.active.max(0) > thresh
    if not np.any(a2):
        a2[0] = True
    Q = kl_components(m1, m2).T
    sns.heatmap(Q[a1][:, a2],
                yticklabels=np.arange(a1.size)[a1],
                xticklabels=np.arange(a2.size)[a2],
                vmin=0, vmax=20, cmap='Greys_r',
                linewidths=0.1, linecolor='k'
               )
    plt.title('Component KL')
    plt.xlabel(m2.name)
    plt.ylabel(m1.name)


def average_ld(m1, m2, L):
    Q = np.zeros((m1.dims['K'], m2.dims['K']))
    if L is None:
        return Q
    for k1 in range(m1.dims['K']):
        for k2 in range(m2.dims['K']):
            s1 = np.random.choice(m1.dims['N'], 10, replace=True, p=m1.pi[k1])
            s2 = np.random.choice(m2.dims['N'], 10, replace=True, p=m2.pi[k2])
            Q[k1, k2] = np.einsum('ms,ms->s', L[:, s1],  L[:, s2]).mean()
    return Q


def average_ld_heatmap(m1, m2, L):
    a1 = m1.active.max(0) > 0.5
    if not np.any(a1):
        a1[0] = True
    a2 = m2.active.max(0) > 0.5
    if not np.any(a2):
        a2[0] = True

    Q = average_ld(m1, m2, L)
    sns.heatmap(Q[a1][:, a2],
                yticklabels=np.arange(a1.size)[a1],
                xticklabels=np.arange(a2.size)[a2],
                center=0, cmap='RdBu_r'
               )
    plt.title('Average LD')
    plt.xlabel(m2.name)
    plt.ylabel(m1.name)


def average_r2(m1, m2, L):
    Q = np.zeros((m1.dims['K'], m2.dims['K']))
    for k1 in range(m1.dims['K']):
        for k2 in range(m2.dims['K']):
            s1 = np.random.choice(m1.dims['N'], 10, replace=True, p=m1.pi[k1])
            s2 = np.random.choice(m2.dims['N'], 10, replace=True, p=m2.pi[k2])
            Q[k1, k2] = (np.einsum('ms,ms->s', L[:, s1],  L[:, s2])**2).mean()
    return Q


def average_r2_heatmap(m1, m2, L):
    a1 = m1.active.max(0) > 0.5
    if not np.any(a1):
        a1[0] = True
    a2 = m2.active.max(0) > 0.5
    if not np.any(a2):
        a2[0] = True

    Q = average_ld(m1, m2, L)
    sns.heatmap(Q[a1][:, a2],
                yticklabels=np.arange(m1.dims['K'])[a1],
                xticklabels=np.arange(m2.dims['K'])[a2],
                vmin=0, vmax=1, cmap='Greys',
                linewidths=0.1, linecolor='k',
               )
    plt.title('Average R2')
    plt.xlabel(m2.name)
    plt.ylabel(m1.name)

kl_sum = lambda A1, A2, k1, k2: np.sum(
    [categorical_kl(A1[:, t, k1], A2[:, t, k2]) for t in range(A1.shape[1] )])

def active_kl(m1, m2):
    A1 = np.stack([m1.active, 1 - m1.active])
    A2 = np.stack([m2.active, 1 - m2.active])
    kls = np.array([[
        kl_sum(A1, A2, k1, k2) + kl_sum(A2, A1, k2, k1)
        for k1 in range(m1.dims['K'])] for k2 in range(m2.dims['K'])])
    return kls


def _active_overlap(a, b):
    """
    active_in_both / active_in_either
    return 0 if none active
    """ 
    a = a > 0.5
    b = b > 0.5
    active_in_both = (a & b).sum()
    active_in_either = (a | b).sum()
    if active_in_either > 0:
        return active_in_both / active_in_either
    return 0


def active_overlap(m1, m2):
    """
    active_in_both / active_in_either
    return 0 if none active
    """ 
    overlap = np.array([[
        _active_overlap(a, b)
        for a in m1.active.T] for b in m2.active.T])
    return overlap

    a = m1.active[:, 0] > 0.5
    b = m1.active[:, 0] > 0.5
    active_in_both = (a & b).sum()
    active_in_either = (a | b).sum()
    if active_in_either > 0:
        return active_in_both / active_in_either
    else:
        return 0


def active_overlap_heatmap(m1, m2):
    a1 = m1.active.max(0) > 0.5
    if not np.any(a1):
        a1[0] = True
    a2 = m2.active.max(0) > 0.5
    if not np.any(a2):
        a2[0] = True
    Q = active_overlap(m1, m2).T
    sns.heatmap(Q[a1][:, a2],
                yticklabels=np.arange(m1.dims['K'])[a1],
                xticklabels=np.arange(m2.dims['K'])[a2],
                vmin=0, vmax=1, cmap='Greys',
                linewidths=0.1, linecolor='k'
               )
    plt.title('Component Activity Overlap')
    plt.xlabel(m2.name)
    plt.ylabel(m1.name)


def active_kl_heatmap(m1, m2, thresh=0.5):
    a1 = m1.active.max(0) > thresh
    if not np.any(a1):
        a1[0] = True
    a2 = m2.active.max(0) > thresh
    if not np.any(a2):
        a2[0] = True
    Q = active_kl(m1, m2).T
    sns.heatmap(Q[a1][:, a2],
                yticklabels=np.arange(m1.dims['K'])[a1],
                xticklabels=np.arange(m2.dims['K'])[a2],
                vmin=0, vmax=20, cmap='Greys_r',
                linewidths=0.1, linecolor='k'
               )
    plt.title('Component Activity KL')
    plt.xlabel(m2.name)
    plt.ylabel(m1.name)


def comparison_heatmaps(m1, m2, L):
    fig, ax = plt.subplots(1, 3, figsize=(16, 4))
    plt.sca(ax[0])
    active_overlap_heatmap(m1, m2)
    plt.sca(ax[1])
    kl_heatmap(m1, m2)
    plt.sca(ax[2])
    average_r2_heatmap(m1, m2, L)
