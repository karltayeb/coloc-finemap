from ukbb_cafeh_ss import *


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

    print('Using z scores...')
    B = B.values / S.values
    S = np.ones_like(B)

    z = (B/S)[0]
    X = center_mean_impute(gtex_genotype.loc[:, variants]).values
    X = (X / np.sqrt((X**2).sum(0))).T
    pd.Series(z, index=variants).to_csv(snakemake.output.z, sep='\t')
    pd.DataFrame(X, index=variants).to_csv(snakemake.output.X, sep='\t')

