import pysam
import pandas as pd
import numpy as np
import sys
sys.path.append('/work-zfs/abattle4/karl/cosie_analysis/utils/')
from misc import load_gtex_genotype
from cafeh.cafeh_ss import CAFEH as CSS
from cafeh.fitting import weight_ard_active_fit_procedure, fit_all
from cafeh.model_queries import summary_table, coloc_table

sample_ld = lambda g: np.corrcoef(center_mean_impute(g), rowvar=False)

gc = pd.read_csv('output/annotations/gencode/gencode_v29_v19.tsv', sep='\t')
gene2tss = gc.set_index('gene_id').start_pos19.to_dict()

gc26 = pd.read_csv(
    '/work-zfs/abattle4/lab_data/annotation/gencode.v26/gencode.v26.annotation.gene.txt', sep='\t')
gene2chr = gc26.set_index('gene_id').chr.to_dict()


print('gene2chr', gene2chr.get(wildcards.gene))

def cast(s):
    try:
        return ast.literal_eval(s)
    except Exception:
        return s

def load_ukbb_gwas(phenotype, gene):
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

    ukbb = pysam.TabixFile(input.sumstats)
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

    df = df.loc[:, ['tissue', 'chr', 'pos', 'ref', 'alt', 'rsid', 'variant_id', 'slope', 'slope_se', 'S', 'z', 'zS']]
    return df

def load_phecode_gwas(phenotype, gene):
    ukbb = pysam.TabixFile(input.sumstats)
    chrom = int(gene2chr.get(gene)[3:])
    left = gene2left.get(gene)
    right = gene2right.get(gene)
    if (right - left) > 1e7:
        right = left + 1e7
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

    df = df.loc[:, ['tissue', 'chr', 'pos', 'ref', 'alt', 'rsid', 'variant_id', 'sample_size', 'slope', 'slope_se', 'S', 'z', 'zS']]
    return df

def load_gtex_associations(gene):
    """
    ap = '/work-zfs/abattle4/karl/cosie_analysis/output/GTEx/{}/{}/{}.associations'.format(
        gene2chr.get(gene), gene, gene)
    v2rp = '/work-zfs/abattle4/karl/cosie_analysis/output/GTEx/{}/{}/{}.snp2rsid.json'.format(
        gene2chr.get(gene), gene, gene)
    """
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
    df = df.loc[:, ['tissue', 'chr', 'pos', 'ref', 'alt', 'rsid', 'variant_id', 'slope', 'slope_se', 'S', 'z', 'zS']]
    return df

def make_table(model, gene, rsid2variant_id):
    table = summary_table(model)

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

gene = wildcards.gene
study = wildcards.study
phenotype = wildcards.phenotype

# load summary stats
gtex = load_gtex_associations(gene)
rsid2variant_id = gtex.set_index('rsid').variant_id.to_dict()

if study == 'UKBB':
    gwas = load_phecode_gwas(phenotype, gene)
if study == 'UKBB_continuous':
    gwas = load_ukbb_gwas(phenotype, gene)

# load genotype
gtex_genotype = load_gtex_genotype(gene, use_rsid=True)
gtex_genotype = gtex_genotype.loc[:,~gtex_genotype.columns.duplicated()]

# flip variants with swapped ref/alt alleles
# remove variants with mismatched ref/alt
print('harmonizing GWAS and GTEx')
a = gtex[~gtex.rsid.duplicated()].set_index('rsid').loc[:, ['ref', 'alt']]
b = gwas[~gwas.rsid.duplicated()].set_index('rsid').loc[:, ['ref', 'alt']]
c = pd.concat([a, b], axis=1, join='inner')

correct = (c.iloc[:, 1] == c.iloc[:, 3]) & (c.iloc[:, 0] == c.iloc[:, 2])
flipped = (c.iloc[:, 1] == c.iloc[:, 2]) & (c.iloc[:, 0] == c.iloc[:, 3])
bad = ~(correct | flipped)

gwas.loc[gwas.rsid.isin(flipped[flipped].index), 'z'] \
    = gwas.loc[gwas.rsid.isin(flipped[flipped].index)].z * -1

shared_variants = c[~bad].index.values

# combine summary stat
df = pd.concat([gtex, gwas])
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

print('{} GTEx variants'.format(gtex.rsid.unique().size))
print('{} UKBB variants'.format(gwas.rsid.unique().size))
print('{} intersecting, fully observed variants'.format(variants.size))

K = snakemake.params.K

if snakemake.params.zscore:
    B = B.values / S.values
    S = np.ones_like(B.values)
else:
    B = B.values
    S = S.values

init_args = {
    'LD': sample_ld(gtex_genotype.loc[:, variants]),
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
table = make_table(css, gene, rsid2variant_id)
table.to_csv(output.variant_report, sep='\t', index=False)

ct = coloc_table(css, phenotype, gene=gene)
ct.to_csv(output.coloc_report, sep='\t', index=False)

css.save(output.model)

