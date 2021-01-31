import pysam
import pandas as pd
import numpy as np
import sys
sys.path.append('/work-zfs/abattle4/karl/cosie_analysis/utils/')
from misc import load_gtex_genotype, load_var2rsid, center_mean_impute
from cafeh.cafeh_ss import CAFEH as CSS
from cafeh.fitting import weight_ard_active_fit_procedure, fit_all
from cafeh.model_queries import summary_table, coloc_table

sample_ld = lambda g: np.corrcoef(center_mean_impute(g), rowvar=False)

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

def load_cad_gwas(gene, rel = ''):
    gwas = pysam.TabixFile(rel + 'output/CAD/CAD/CAD.tsv.bgz')
    gc = pd.read_csv(rel + 'output/annotations/genes_hg19.bed', sep='\t')
    gc.loc[:, 'left'] = np.maximum(0, gc.start - 1e6)
    gc.loc[:, 'right'] = gc.end + 1e6

    gene2chr = gc.set_index('gene_id').chrom.to_dict()
    gene2left = gc.set_index('gene_id').left.to_dict()
    gene2right = gc.set_index('gene_id').right.to_dict()

    left = gene2left.get(gene)
    right = gene2right.get(gene)
    chrom = gene2chr.get(gene)
    print('range: ', chrom, left, right)
    lines = gwas.fetch(chrom, left, right)

    header =[
        '#CHR', 'POS', 'ID', 'REF', 'ALT',
        'ref_frequency', 'beta', 'beta_se',
        'p', 'sample_size', 'chromosome_number']
    df = pd.DataFrame(
        list(map(cast, x.strip().split('\t')) for x in lines),
        columns=header
    )
    
    df = df.apply(pd.to_numeric, errors='ignore')
    df = df.loc[:,~df.columns.duplicated()]

    df.beta = pd.to_numeric(df.beta, errors='coerce')
    df.beta_se = pd.to_numeric(df.beta_se, errors='coerce')
    df.p = pd.to_numeric(df.p, errors='coerce')
            
    df.rename(columns={
        '#CHR': 'chr',
        'POS': 'pos',
        'REF': 'ref',
        'ALT': 'alt',
        'ID': 'rsid',
        'beta': 'slope',
        'beta_se': 'slope_se',
        'p': 'pval_nominal'}, inplace=True)

    df.loc[:, 'ref'] = df.ref.str.upper()
    df.loc[:, 'alt'] = df.alt.str.upper()

    df.loc[:, 'tissue'] = 'CAD'
    df.loc[:, 'gene'] = gene
    df.loc[:, 'z'] = df.slope / df.slope_se
    df.loc[:, 'zS'] = np.sqrt((df.z**2 / df.sample_size) + 1)
    df.loc[:, 'S'] = df.slope_se # large sample approximation

    df = df.loc[:, COLUMNS]
    return df

def load_ukbb_gwas(phenotype, gene, rel = ''):

    gc = pd.read_csv(rel + 'output/annotations/genes_hg19.bed', sep='\t')
    gc.loc[:, 'left'] = np.maximum(0, gc.start - 1e6)
    gc.loc[:, 'right'] = gc.end + 1e6

    gene2chr = gc.set_index('gene_id').chrom.to_dict()
    gene2left = gc.set_index('gene_id').left.to_dict()
    gene2right = gc.set_index('gene_id').right.to_dict()


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

    #ukbb = pysam.TabixFile(snakemake.input.sumstats)
    ukbb = pysam.TabixFile(rel + 'output/UKBB_continuous/{}/{}.tsv.bgz'.format(phenotype, phenotype))
    chrom = int(gene2chr.get(gene)[3:])
    left = gene2left.get(gene)
    right = gene2right.get(gene)
    if (right - left) > 1e7:
        right = left + 1e7

    print(chrom, left, right)
    lines = ukbb.fetch(chrom, left, right)

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
    df.loc[:, 'gene'] = gene
    df.loc[:, 'z'] = df.slope / df.slope_se
    df.loc[:, 'zS'] = np.sqrt((df.z**2 / df.sample_size) + 1)
    df.loc[:, 'S'] = np.sqrt((df.slope**2 / df.sample_size) + df.slope_se**2)
    df = df[(df.low_confidence_variant == 'false')]

    df = df.loc[:, COLUMNS]
    return df

def load_phecode_gwas(phenotype, gene, rel=''):
    gc = pd.read_csv(rel + 'output/annotations/genes_hg19.bed', sep='\t')
    gc.loc[:, 'left'] = np.maximum(0, gc.start - 1e6)
    gc.loc[:, 'right'] = gc.end + 1e6

    gene2chr = gc.set_index('gene_id').chrom.to_dict()
    gene2left = gc.set_index('gene_id').left.to_dict()
    gene2right = gc.set_index('gene_id').right.to_dict()

    ukbb = pysam.TabixFile(rel + 'output/UKBB/{}/{}.tsv.bgz'.format(phenotype, phenotype))
    #ukbb = pysam.TabixFile(snakemake.input.sumstats)
    chrom = int(gene2chr.get(gene)[3:])
    left = gene2left.get(gene)
    right = gene2right.get(gene)
    if (right - left) > 1e7:
        right = left + 1e7
    if phenotype in ['Uterine_polyp', 'Uterine_leiomyoma']:
        chrom = "{0:0=2d}".format(chrom)

    print(chrom, left, right)
    lines = ukbb.fetch(chrom, left, right)


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
    df.loc[:, 'gene'] = gene
    df.loc[:, 'z'] = df.slope / df.slope_se
    df.loc[:, 'zS'] = np.sqrt((df.z**2 / df.sample_size) + 1)
    df.loc[:, 'S'] = np.sqrt((df.slope**2 / df.sample_size) + df.slope_se**2)

    df = df.loc[:, COLUMNS]
    return df

def load_grasp_gwas(phenotype, gene, rel=''):
    gc = pd.read_csv(rel + 'output/annotations/genes_hg19.bed', sep='\t')
    gc.loc[:, 'left'] = np.maximum(0, gc.start - 1e6)
    gc.loc[:, 'right'] = gc.end + 1e6

    gene2chr = gc.set_index('gene_id').chrom.to_dict()
    gene2left = gc.set_index('gene_id').left.to_dict()
    gene2right = gc.set_index('gene_id').right.to_dict()

    ukbb = pysam.TabixFile(rel + 'output/GRASP/{}/{}.tsv.bgz'.format(phenotype, phenotype))
    #ukbb = pysam.TabixFile(snakemake.input.sumstats)
    chrom = gene2chr.get(gene)
    left = gene2left.get(gene)
    right = gene2right.get(gene)
    if (right - left) > 1e7:
        right = left + 1e7
    if phenotype in ['Uterine_polyp', 'Uterine_leiomyoma']:
        chrom = "{0:0=2d}".format(chrom)

    print(chrom, left, right)
    lines = ukbb.fetch(chrom, left, right)

    header= [
        'CHR', 'POS', 'ID', 'REF', 'ALT',
        'ref_frequency', 'beta', 'beta_se',
        'p', 'sample_size', 'chromosome_number']

    df = pd.DataFrame(
        list(map(cast, x.strip().split('\t')) for x in lines),
        columns=header
    )

    df.rename(columns={
        'POS': 'pos',
        'REF': 'ref',
        'ALT': 'alt',
        'ID': 'rsid',
        'beta': 'slope',
        'beta_se': 'slope_se',
        'p': 'pval_nominal'}, inplace=True)

    df.loc[:, 'ref'] = df.ref.str.upper()
    df.loc[:, 'alt'] = df.alt.str.upper()

    df = df.apply(pd.to_numeric, errors='ignore')
    df = df.loc[:,~df.columns.duplicated()]
    df.slope = pd.to_numeric(df.slope, errors='coerce')
    df.slope_se = pd.to_numeric(df.slope_se, errors='coerce')
    df.pval_nominal = pd.to_numeric(df.pval_nominal, errors='coerce')


    df.loc[:, 'tissue'] = phenotype
    df.loc[:, 'S'] = np.sqrt((df.slope**2 / df.sample_size) + df.slope_se**2)
    df.loc[:, 'z'] = df.slope / df.slope_se
    df.loc[:, 'zS'] = np.sqrt((df.z**2 / df.sample_size) + 1)

    df = df.loc[:, COLUMNS]
    return df


def load_gtex_associations(gene, rel=''):
    """
    ap = '/work-zfs/abattle4/karl/cosie_analysis/output/GTEx/{}/{}/{}.associations'.format(
        gene2chr.get(gene), gene, gene)
    v2rp = '/work-zfs/abattle4/karl/cosie_analysis/output/GTEx/{}/{}/{}.snp2rsid.json'.format(
        gene2chr.get(gene), gene, gene)
    """
    gc = pd.read_csv(rel + 'output/annotations/genes_hg19.bed', sep='\t')
    gc.loc[:, 'left'] = np.maximum(0, gc.start - 1e6)
    gc.loc[:, 'right'] = gc.end + 1e6

    gene2chr = gc.set_index('gene_id').chrom.to_dict()

    v2r = load_var2rsid(gene)
    ap = '/work-zfs/abattle4/karl/cosie_analysis/output/GTEx/{}/{}/{}.associations'.format(
        gene2chr.get(gene), gene, gene)

    df = pd.read_csv(ap, index_col=0)
    df.loc[:, 'rsid'] = df.variant_id.apply(lambda x: v2r.get(x, '-'))
    df.loc[:, 'pos'] = df.variant_id.apply(lambda x: int(x.split('_')[1]))
    df.loc[:, 'ref'] = df.variant_id.apply(lambda x: x.split('_')[2])
    df.loc[:, 'alt'] = df.variant_id.apply(lambda x: x.split('_')[3])
    df.loc[:, 'sample_size'] = (df.ma_count / df.maf / 2)
    
    df.loc[:, 'S'] = np.sqrt(
        df.slope**2/df.sample_size + df.slope_se**2)
    df.loc[:, 'z'] = df.slope / df.slope_se
    df.loc[:, 'zS'] = np.sqrt((df.z**2 / df.sample_size) + 1)
    df.loc[:, 'S'] = np.sqrt((df.slope**2 / df.sample_size) + df.slope_se**2)
    df.loc[:, 'chr'] = gene2chr.get(gene)
    df = df.loc[:, COLUMNS]
    return df

def make_table(model, gene, rsid2variant_id, filter_variants=False):
    table = summary_table(model, filter_variants=filter_variants)

    # annotate table
    table.loc[:, 'rsid'] = table.variant_id
    table.loc[:, 'variant_id'] = table.rsid.apply(lambda x: rsid2variant_id.get(x, 'chr0_0_A_B_n'))
    table.loc[:, 'chr'] = table.variant_id.apply(lambda x: (x.split('_')[0]))
    table.loc[:, 'start'] = table.variant_id.apply(lambda x: int(x.split('_')[1]))
    table.loc[:, 'end'] = table.start + 1
    table.loc[:, 'gene'] = gene

    table = table.loc[:, ['chr', 'start', 'end', 'gene',
                          'variant_id', 'rsid', 'study', 'pip',
                          'top_component', 'p_active', 'pi', 'alpha',
                          'rank', 'effect', 'effect_var']]
    return table

if __name__ == "__main__":
    gene = snakemake.wildcards.gene
    study = snakemake.wildcards.study
    phenotype = snakemake.wildcards.phenotype

    # load summary stats
    gtex = load_gtex_associations(gene)
    rsid2variant_id = gtex.set_index('rsid').variant_id.to_dict()

    if study == 'UKBB':
        gwas = load_phecode_gwas(phenotype, gene)
    if study == 'UKBB_continuous':
        gwas = load_ukbb_gwas(phenotype, gene)
    if study == 'CAD':
        gwas = load_cad_gwas(gene)
    if study == 'GRASP':
        gwas = load_grasp_gwas(phenotype, gene)

    # load genotype
    gtex_genotype = load_gtex_genotype(gene, use_rsid=True)
    gtex_genotype = gtex_genotype.loc[:,~gtex_genotype.columns.duplicated()]


    print('{} GTEx variants'.format(gtex.rsid.unique().size))
    print('{} UKBB variants'.format(gwas.rsid.unique().size))

    # flip variants with swapped ref/alt alleles
    # remove variants with mismatched ref/alt
    print('harmonizing GWAS and GTEx')
    a = gtex[~gtex.rsid.duplicated()].set_index('rsid').loc[:, ['ref', 'alt']]
    b = gwas[~gwas.rsid.duplicated()].set_index('rsid').loc[:, ['ref', 'alt']]
    c = pd.concat([a, b], axis=1, join='inner')

    correct = (c.iloc[:, 1] == c.iloc[:, 3]) & (c.iloc[:, 0] == c.iloc[:, 2])
    flipped = (c.iloc[:, 1] == c.iloc[:, 2]) & (c.iloc[:, 0] == c.iloc[:, 3])
    bad = ~(correct | flipped)
    print('Correct: {}, Flipped: {}, Bad {}'.format(correct.sum(), flipped.sum(), bad.sum()))

    gwas.loc[gwas.rsid.isin(flipped[flipped].index), 'slope'] \
        = gwas.loc[gwas.rsid.isin(flipped[flipped].index)].slope * -1

    shared_variants = c[~bad].index.values

    # combine summary stat
    df = pd.concat([gtex, gwas])
    df = df[df.rsid.isin(shared_variants)]

    # reshape summary stats for CAFEH
    B = df[~df.duplicated(['tissue', 'rsid'])].pivot('tissue', 'rsid', 'slope')
    S = df[~df.duplicated(['tissue', 'rsid'])].pivot('tissue', 'rsid', 'S')


    model = pickle.load(open(snakemake.input.model, 'rb'))
    model._decompress_model()

    variants = model.snp_ids
    study_ids = model.study_ids

    #  remove variants with missing sumstats
    #mask = ~(np.any(np.isnan(B), 0) | np.any(np.isnan(S), 0))
    #B = B.loc[:, mask]
    #S = S.loc[:, mask]
    B = B.loc[:, variants]
    S = S.loc[:, variants]

    B = B.values / S.values
    S = np.ones_like(B)
    LD = sample_ld(gtex_genotype.loc[:, variants])


    model.B = B
    model.S = S
    model.LD = LD

    model.run_time
    model.fit(update_weights=True, update_pi=False, verbose=True, max_iter=2)

    table = make_table(model, gene, rsid2variant_id, filter_variants=False)
    table.to_csv(snakemake.output.variant_report, sep='\t', index=False)

