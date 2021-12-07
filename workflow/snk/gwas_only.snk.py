
source2bfile = {
    'gtex': '/work-zfs/abattle4/marios/GTEx_v8/coloc/GTEx_all_genotypes',
    '1kg': '/work-zfs/abattle4/marios/annotations/1kG_plink/1000G_hg38_plink_merged'
}

wildcard_constraints:
    source="1kg|gtex"

rule get_gtex_genotype_for_gwas:
    output:
        snplist = \
        'output/GWAS_only/{study}/{phenotype}/{chr}/{locus}/{phenotype}.{locus}.{source}.snplist',
        genotype = 'output/GWAS_only/{study}/{phenotype}/{chr}/{locus}/{phenotype}.{locus}.{source}.raw',
        log = 'output/GWAS_only/{study}/{phenotype}/{chr}/{locus}/{phenotype}.{locus}.{source}.log'
    group: "g"
    params:
        bfile = lambda w: source2bfile.get(w)
    run:
        from utils.misc import plink_get_genotype
        import subprocess

        def load_lookup():
            df = pd.read_csv(
                '/work-zfs/abattle4/marios/GTEx_v8/coloc/aggregate_phenotype_loci', sep='\t', header=None)
            df.columns = ['chrom', 'start', 'end', 'rsid', 'phenotype']
            df = df.set_index(['chrom', 'rsid'])
            return df

        def plink_get_genotype_gwas_only(lookup, bfile, save_path):
            chrom = save_path.split('/')[-3]
            rsid = save_path.split('/')[-2]
            start, end = lookup.loc[(chrom, rsid)].values[0][:2]

            cmd = ' '.join(
                ['plink',
                 '--bfile', bfile,
                 '--chr', chrom[3:],
                 '--from-bp', str(start),
                 '--to-bp', str(end),
                 '--maf', '0.01',
                 '--geno', '0.1',
                 '--recode', 'A',
                 '--keep-allele-order',
                 '--snps-only', '--write-snplist', '--allow-no-sex',
                 '--out', save_path])
            return cmd

        lookup = load_lookup()
        q = len(params.source) + 5
        print(params.bfile)
        print(output.genotype[:-q])
        cmd = plink_get_genotype_gwas_only(lookup, params.bfile, output.genotype[:-q])
        print(cmd)
        shell(cmd)
        print('PLINK FINISHED RUNNING?')

rule snpid2rsid_for_gwas:
    input:
        'output/GWAS_only/{study}/{phenotype}/{chr}/{locus}/{phenotype}.{locus}.{source}.snplist'
    output:
       rsids = 'output/GWAS_only/{study}/{phenotype}/{chr}/{locus}/{phenotype}.{locus}.{source}.rsids',
       rsid_map = 'output/GWAS_only/{study}/{phenotype}/{chr}/{locus}/{phenotype}.{locus}.{source}.snp2rsid'
    group: "g"
    script:
        "../../workflow/scripts/variantid2rsid.py"


rule fit_gwas_z_cafeh:
    input:
        genotype_gtex = 'output/GWAS_only/{study}/{phenotype}/{chr}/{locus}/{phenotype}.{locus}.{source}.raw',
        sumstats='output/{study}/{phenotype}/{phenotype}.tsv.bgz',
        tabix_index='output/{study}/{phenotype}/{phenotype}.tsv.bgz.tbi',
        v2r = 'output/GWAS_only/{study}/{phenotype}/{chr}/{locus}/{phenotype}.{locus}.{source}.snp2rsid'
    output:
        variant_report=\
            'output/GWAS_only/{study}/{phenotype}/{chr}/{locus}/{phenotype}.{locus}.{source}.z.variant_report',
        model=\
            'output/GWAS_only/{study}/{phenotype}/{chr}/{locus}/{phenotype}.{locus}.{source}.z.css'
    params:
        K=20,
        zscore=True,
        zld=False
    group: 'report'
    script:
        '../../workflow/scripts/ukbb_cafeh_ss.py'

rule fit_gwas_z_cafeh_zld:
    input:
        genotype_gtex = 'output/GWAS_only/{study}/{phenotype}/{chr}/{locus}/{phenotype}.{locus}.{source}.raw',
        sumstats='output/{study}/{phenotype}/{phenotype}.tsv.bgz',
        tabix_index='output/{study}/{phenotype}/{phenotype}.tsv.bgz.tbi',
        v2r = 'output/GWAS_only/{study}/{phenotype}/{chr}/{locus}/{phenotype}.{locus}.{source}.snp2rsid'
    output:
        variant_report=\
            'output/GWAS_only/{study}/{phenotype}/{chr}/{locus}/{phenotype}.{locus}.{source}.zld.z.variant_report',
        model=\
            'output/GWAS_only/{study}/{phenotype}/{chr}/{locus}/{phenotype}.{locus}.{source}.zld.z.css'
    params:
        K=20,
        zscore=True,
        zld=True
    group: 'report'
    script:
        '../../workflow/scripts/ukbb_cafeh_ss.py'