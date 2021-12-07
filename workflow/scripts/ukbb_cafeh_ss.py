import pysam
import json
import pandas as pd
import numpy as np
import sys
sys.path.append('/work-zfs/abattle4/karl/cosie_analysis/utils/')
from misc import load_gtex_genotype, load_var2rsid, center_mean_impute
from cafeh.cafeh_ss import CAFEH as CSS
from cafeh.fitting import weight_ard_active_fit_procedure, fit_all
from cafeh.model_queries import summary_table, coloc_table

COLUMNS = [
    'tissue', 'chr', 'pos', 'ref', 'alt',
    'rsid', 'variant_id',
    'sample_size', 'slope', 'slope_se',
    'S', 'z', 'zS', 'pval_nominal'
]

def cast(s):
    try:
        return ast.literal_eval(s)
    except Exception:
        return s

def load_lookup(phenotype=None):
    df = pd.read_csv(
        '/work-zfs/abattle4/marios/GTEx_v8/coloc/aggregate_phenotype_loci', sep='\t', header=None)
    df.columns = ['chrom', 'start', 'end', 'rsid', 'phenotype']
    df = df.set_index(['rsid'])
    
    continuous_phenotypes = [
        'BMI', 'Diastolic_blood_pressure',
        'Systolic_blood_pressure', 'LDL',
        'Lipoprotein_A', 'QRS_duration', 'Ventricular_rate'
    ]
    df['source'] = df.phenotype.apply(lambda x: 'UKBB' if x not in continuous_phenotypes else 'UKBB_continuous')
    if phenotype is not None:
        df = df[df.phenotype == phenotype]
    return df


def load_bim():
    return pd.read_csv(snakemake.input.bim, sep = '\t', header=None)\
        .rename(columns={
            0: 'chrom',
            1: 'variant_id',
            2: 'cm',
            3: 'pos',
            4: 'ref',
            5: 'alt'})

def load_gtex_genotype2(phenotype, locus, use_rsid=False):
    """
    load gtex genotype for variants in 1Mb window of gene tss
    @param gene: ensemble gene id
    @param use_rsid: boolean to convert variant id to rsid
    """
    lookup = load_lookup(phenotype)
    
    d = lookup.loc[locus].to_dict()
    d['locus'] = locus
    d['geno_source'] = 'gtex'
    gp = 'output/GWAS_only/{source}/{phenotype}/{chrom}/{locus}/{phenotype}.{locus}.{geno_source}.raw'.format(**d)
    v2rp = 'output/GWAS_only/{source}/{phenotype}/{chrom}/{locus}/{phenotype}.{locus}.{geno_source}.snp2rsid'.format(**d)
    v2r = json.load(open(v2rp, 'r'))

    print('loading gtex genotypes...')
    genotype = pd.read_csv(gp, sep=' ')
    genotype = genotype.set_index('IID').iloc[:, 5:]

    # recode genotypes
    coded_snp_ids = np.array([x.strip() for x in genotype.columns])
    snp_ids = {x: '_'.join(x.strip().split('_')[:-1]) for x in coded_snp_ids}
    genotype.rename(columns=snp_ids, inplace=True)

    if use_rsid:
        genotype.rename(columns=v2r, inplace=True)
    return genotype, v2r


def load_phecode_gwas(phenotype, locus, rel=''):
    ukbb = pysam.TabixFile(rel + 'output/UKBB/{}/{}.tsv.bgz'.format(phenotype, phenotype))
    lookup = load_lookup(phenotype)
    chrom = int(lookup.loc[locus].chrom[3:])
    left = lookup.loc[locus].start
    right = lookup.loc[locus].end
    if phenotype in ['Uterine_polyp', 'Uterine_leiomyoma']:
        chrom = "{0:0=2d}".format(chrom)

    print(chrom, left, right)
    #lines = ukbb.fetch(chrom, left, right)
    lines = ukbb.fetch(chrom, 0, 1e20)

    header = [
        'CHROM', 'POS', 'ID', 'REF', 'ALT',
        'ac', 'af', 'num_cases', 'num_controls',
        'beta', 'sebeta', 'Tstat', 'pval', 'pval_SAIGE_NoSPA',
        'Is_Converged', 'varT', 'varTstar']

    df = pd.DataFrame(
        list(map(cast, x.strip().split('\t')) for x in lines),
        columns=header
    )


    df = df.apply(pd.to_numeric, errors='ignore')
    df = df.loc[:,~df.columns.duplicated()]

    df.beta = pd.to_numeric(df.beta, errors='coerce')
    df.sebeta = pd.to_numeric(df.sebeta, errors='coerce')
    df.Tstat = pd.to_numeric(df.Tstat, errors='coerce')
    df.pval = pd.to_numeric(df.pval, errors='coerce')

    df.rename(columns={
        'CHROM': 'chr',
        'POS': 'pos',
        'REF': 'ref',
        'ALT': 'alt',
        'ID': 'rsid',
        'beta': 'slope',
        'sebeta': 'slope_se'
    }, inplace=True)

    # TODO: Effective Sample size?
    df.loc[:, 'sample_size'] = df.num_cases + df.num_controls

    df.loc[:, 'tissue'] = phenotype
    df.loc[:, 'z'] = df.slope / df.slope_se
    df.loc[:, 'zS'] = np.sqrt((df.z**2 / df.sample_size) + 1)
    df.loc[:, 'S'] = np.sqrt((df.slope**2 / df.sample_size) + df.slope_se**2)

    df = df.loc[:, COLUMNS]
    return df


def load_ukbb_gwas(phenotype, locus, rel = ''):
    header = [
        'variant', 'chr', 'pos', 'ref', 'alt', 'rsid',
        'varid', 'consequence', 'consequence_category',
        'info', 'call_rate', 'AC', 'AF', 'minor_allele', 
        'minor_AF', 'p_hwe', 'n_called', 'n_not_called',
        'n_hom_ref', 'n_het', 'n_hom_var', 'n_non_ref',
        
        'r_heterozygosity', 'r_het_hom_var', 'r_expected_het_frequency',
        'variant', 'minor_allele', 'minor_AF',
        'low_confidence_variant', 'n_complete_samples',
        'AC', 'ytx', 'beta', 'se', 'tstat', 'pval'
    ]

    header_cc = [
        'variant', 'chr', 'pos', 'ref', 'alt', 'rsid',
        'varid', 'consequence', 'consequence_category',
        'info', 'call_rate', 'AC', 'AF', 'minor_allele', 
        'minor_AF', 'p_hwe', 'n_called', 'n_not_called',
        'n_hom_ref', 'n_het', 'n_hom_var', 'n_non_ref',
        
        'r_heterozygosity', 'r_het_hom_var', 'r_expected_het_frequency',
        'variant', 'minor_allele', 'minor_AF', 'expected_case_minor_AC',
        'low_confidence_variant', 'n_complete_samples',
        'AC', 'ytx', 'beta', 'se', 'tstat', 'pval'
    ]

    ukbb = pysam.TabixFile(rel + 'output/UKBB_continuous/{}/{}.tsv.bgz'.format(phenotype, phenotype))
    
    lookup = load_lookup(phenotype)
    chrom = int(lookup.loc[locus].chrom[3:])
    left = lookup.loc[locus].start
    right = lookup.loc[locus].end

    print(chrom, left, right)
    #lines = ukbb.fetch(chrom, left, right)
    lines = ukbb.fetch(chrom, 0, 1e20)

    df = pd.DataFrame(
        list(map(cast, x.strip().split('\t')) for x in lines),
    )

    if df.shape[1] == len(header):
        df.columns = header
    else:
        df.columns = header_cc

    df = df.apply(pd.to_numeric, errors='ignore')
    df = df.loc[:,~df.columns.duplicated()]

    df.beta = pd.to_numeric(df.beta, errors='coerce')
    df.se = pd.to_numeric(df.se, errors='coerce')
    df.tstat = pd.to_numeric(df.tstat, errors='coerce')
    df.pval = pd.to_numeric(df.pval, errors='coerce')

    df.rename(columns={
        'varid': 'variant_id',
        'beta': 'slope',
        'se': 'slope_se',
        'pval': 'pval_nominal',
        'n_complete_samples': 'sample_size'
    }, inplace=True)
    df.loc[:, 'tissue'] = phenotype
    df.loc[:, 'z'] = df.slope / df.slope_se
    df.loc[:, 'zS'] = np.sqrt((df.z**2 / df.sample_size) + 1)
    df.loc[:, 'S'] = np.sqrt((df.slope**2 / df.sample_size) + df.slope_se**2)
    df = df[(df.low_confidence_variant == 'false')]

    df = df.loc[:, COLUMNS]
    return df

def make_table(model, locus, rsid2variant_id, bim):
    table = summary_table(model)
    bim = bim[~bim.rsid.duplicated()].set_index('rsid')

    # annotate table
    table['rsid'] = table.variant_id
    table['variant_id'] = table.rsid.apply(lambda x: rsid2variant_id.get(x, 'chr0_0_A_B_n'))
    table['chr'] = bim.loc[table.rsid].chrom # table.variant_id.apply(lambda x: (x.split('_')[0]))
    table['start'] = bim.loc[table.rsid].pos #table.variant_id.apply(lambda x: int(x.split('_')[1]))
    table['end'] = table.start + 1
    table['sentinal_snp'] = locus

    table = table.loc[:, ['chr', 'start', 'end', 'sentinal_snp',
                          'variant_id', 'rsid', 'study', 'pip',
                          'top_component', 'p_active', 'pi', 'alpha',
                          'rank', 'effect', 'effect_var']]
    return table


if __name__ == "__main__":
    locus = snakemake.wildcards.locus
    study = snakemake.wildcards.study
    phenotype = snakemake.wildcards.phenotype

    if study == 'UKBB':
        gwas = load_phecode_gwas(phenotype, locus)
    if study == 'UKBB_continuous':
        gwas = load_ukbb_gwas(phenotype, locus)

    # load genotype
    gtex_genotype, v2r = load_gtex_genotype2(phenotype, locus, use_rsid=True)
    gtex_genotype = gtex_genotype.loc[:,~gtex_genotype.columns.duplicated()]
    rsid2variant_id = {v: k for k, v in v2r.items()}
    
    bim = load_bim()
    bim['rsid'] = bim.variant_id.apply(lambda x: v2r.get(x, '-'))
    bim = bim[bim.rsid != '-']

    print('{} UKBB variants'.format(gwas.rsid.unique().size))
    print('{} reference variants'.format(bim.rsid.unique().size))

    # flip variants with swapped ref/alt alleles
    # remove variants with mismatched ref/alt
    print('harmonizing GWAS and GTEx')
    a = bim.set_index('rsid')[['ref', 'alt']]
    b = gwas[~gwas.rsid.duplicated()].set_index('rsid').loc[:, ['ref', 'alt']]
    anb = np.intersect1d(a.index, b.index)
    c = pd.concat([a.loc[anb], b.loc[anb]], axis=1, join='inner')
    correct = (c.iloc[:, 1] == c.iloc[:, 3]) & (c.iloc[:, 0] == c.iloc[:, 2])
    flipped = (c.iloc[:, 1] == c.iloc[:, 2]) & (c.iloc[:, 0] == c.iloc[:, 3])
    bad = ~(correct | flipped)
    print('Correct: {}, Flipped: {}, Bad {}'.format(correct.sum(), flipped.sum(), bad.sum()))

    gwas.loc[gwas.rsid.isin(flipped[flipped].index), 'slope'] \
        = gwas.loc[gwas.rsid.isin(flipped[flipped].index)].slope * -1

    shared_variants = c[~bad].index.values
    shared_variants = np.intersect1d(gtex_genotype.columns, shared_variants)
    shared_variants = np.intersect1d(b.index, shared_variants)

    # combine summary stat
    df = gwas
    df = df[df.rsid.isin(shared_variants)]

    # reshape summary stats for CAFEH
    B = df[~df.duplicated(['tissue', 'rsid'])].pivot('tissue', 'rsid', 'slope')
    S = df[~df.duplicated(['tissue', 'rsid'])].pivot('tissue', 'rsid', 'S')

    # remove variants with missing sumstats
    mask = ~(np.any(np.isnan(B), 0) | np.any(np.isnan(S), 0))
    B = B.loc[:, mask]
    S = S.loc[:, mask]

    variants = B.columns.values
    study_ids = B.index.values

    print('{} intersecting, fully observed variants'.format(variants.size))

    K = snakemake.params.K

    if snakemake.params.zscore:
        print('Using z scores...')
        B = B.values / S.values
        S = np.ones_like(B)
    else:
        print('Using effect sizes...')
        B = B.values
        S = S.values

    if snakemake.params.zld:
        LD = np.corrcoef(
            np.concatenate(
                [center_mean_impute(gtex_genotype.loc[:, variants]).values, B/S]),
            rowvar=False
        )
    else:
        LD = np.corrcoef(
            center_mean_impute(gtex_genotype.loc[:, variants]).values,
            rowvar=False
        )

    init_args = {
        'LD': LD,
        'B': B,
        'S': S,
        'K': K,
        'snp_ids': variants,
        'study_ids': study_ids,
        'tolerance': 1e-8
    }

    css = CSS(**init_args)
    css.prior_activity = np.ones(K) * 0.1
    css.weight_precision_b = np.ones_like(css.weight_precision_b) * 1

    print('fitting CAFEH')
    weight_ard_active_fit_procedure(css, max_iter=10, verbose=True)
    fit_all(css, max_iter=30, verbose=True)

    # save variant report
    table = make_table(css, locus, rsid2variant_id, bim)
    table.to_csv(snakemake.output.variant_report, sep='\t', index=False)

    css.save(snakemake.output.model)
