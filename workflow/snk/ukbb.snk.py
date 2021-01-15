# UKBB_continuous

rule download_ukbb_pheno:
    input:
        manifest='output/UKBB_continuous/manifest.txt',
        pheno2manifest='output/UKBB_continuous/pheno2manifest'
    output:
        temp_sumstats='output/UKBB_continuous/{phenotype}/_{phenotype}.tsv.bgz'
    run:
        import pandas as pd
        import subprocess
        import json

        manifest = pd.read_csv(input.manifest, sep='\t').set_index(['Phenotype Code', 'Sex'])
        pheno2manifest = json.load(open(input.pheno2manifest, 'r'))

        print('downloading summary stats')
        phenotype = wildcards.phenotype
        source_file = manifest.loc[:, 'wget command']\
            .loc[pheno2manifest[phenotype]]\
            .values[0]\
            .split(' ')[1]

        cmd = 'wget {} -O {}'.format(source_file, output.temp_sumstats)
        print(cmd)
        subprocess.run(cmd, shell=True)

rule ukbb_build_index:
    input:
        var='output/UKBB_continuous/variants.tsv.bgz',
        sumstats_raw='output/UKBB_continuous/{phenotype}/_{phenotype}.tsv.bgz'
    output:
        sumstats='output/UKBB_continuous/{phenotype}/{phenotype}.tsv.bgz',
        tabix_index='output/UKBB_continuous/{phenotype}/{phenotype}.tsv.bgz.tbi'
    run:
        shell("paste <(zcat {input.var}) <(zcat {input.sumstats_raw}) | bgzip > {output.sumstats}")
        shell("tabix -s 2 -b 3 -e 3 -S 1 {output.sumstats}")

study2col_hits = {
    'UKBB_continuous': (36, 6),
    'UKBB': (13, 3),
    'CAD': (9, 3)
}
rule ukbb_get_hits:
    input:
        sumstats='output/{study}/{phenotype}/{phenotype}.tsv.bgz'
    output:
        hits='output/{study}/{phenotype}/{phenotype}.hits.txt'
    params:
        pcol = lambda wildcards: study2col_hits[wildcards.study][0],
        rsidcol = lambda wildcards: study2col_hits[wildcards.study][1]
    run:
        shell("grep -w -F -f <(zcat {input.sumstats} | awk '{{if(${params.pcol} < 1e-6){{print ${params.rsidcol}}}}}') /work-zfs/abattle4/marios/annotations/1kG_plink/1000G_hg38_plink_merged.bim > {output.hits}")

study2col_request = {
    'UKBB_continuous': 1,
    'UKBB': 0,
    'CAD': 0
}

rule ukbb_get_cis_genes:
    input:
        hits='output/{study}/{phenotype}/{phenotype}.hits.txt'
    output:
        request='output/{study}/{phenotype}/{phenotype}.genes.txt'
    run:
        import pandas as pd
        import numpy as np
        from tqdm import tqdm
        from glob import glob

        hits = pd.read_csv(input.hits, sep='\t', header=None)
        hits = np.array(['chr{}_{}'.format(c, p) for c, p in zip(hits.iloc[:, 0], hits.iloc[:, 3])])

        genes = []
        for file in tqdm(glob('/work-zfs/abattle4/lab_data/GTEx_v8/ciseQTL/GTEx_Analysis_v8_eQTL/*.v8.signif_variant_gene_pairs.txt')):
            sig = pd.read_csv(file, sep='\t')
            sig.loc[:, 'tissue'] = file.split('/')[-1].split('.')[0]
            mask = sig.variant_id.apply(lambda x: '_'.join(x.split('_')[:2])).isin(hits)
            sig = sig[mask]
            sig = sig[~sig.gene_id.duplicated(keep='first')]
            genes.append(sig)

        genes = pd.concat(genes)
        print(genes.gene_id.unique().size)

        results = pd.DataFrame([{
            'study': wildcards.study,
            'phenotype': wildcards.phenotype,
            'tissue': row.tissue,
            'chr': row.variant_id.split('_')[0],
            'pos': row.variant_id.split('_')[1],
            'gene_id': row.gene_id
        } for _, row in genes.iterrows()])
        results.to_csv(output[0], sep='\t')

rule ukbb_get_request:
    input:
        hits='output/{study}/{phenotype}/{phenotype}.hits.txt'
    output:
        request='output/{study}/{phenotype}/{phenotype}.request.txt'
    params:
        base_col = lambda wildcards: study2col_request[wildcards.study]
    run:
        import pandas as pd
        import numpy as np
        from tqdm import tqdm
        from glob import glob

        hits = pd.read_csv(input.hits, sep='\t', header=None)
        hits = np.array(['{}_{}'.format(c, p) for c, p in zip(hits.iloc[:, 0], hits.iloc[:, 3])])
    
        genes = []
        for file in tqdm(glob('/work-zfs/abattle4/lab_data/GTEx_v8/ciseQTL/GTEx_Analysis_v8_eQTL/*.v8.signif_variant_gene_pairs.txt')):
            sig = pd.read_csv(file, sep='\t')
            mask = sig.variant_id.apply(lambda x: '_'.join(x[3:].split('_')[:2])).isin(hits)
            sig = sig[mask]
            sig[~sig.gene_id.duplicated()]
            genes.append(sig)

        genes = pd.concat(genes)
        genes = genes[~genes.gene_id.duplicated()]

        template = 'output/{study}/{phe}/{chr}/{gene}/{gene}.{phe}.z.variant_report'
        requests = [template.format(study=wildcards.study,
                            phe=wildcards.phenotype,
                            chr=row.variant_id.split('_')[0], gene=row.gene_id)
                    for _, row in genes.iterrows()
                   ]

        with open(output.request, 'w') as f:
            for r in requests:
                print(r, file=f)


rule ukbb_gtex_cafeh:
    input:
        genotype_gtex = 'output/GTEx/{chr}/{gene}/{gene}.raw',
        associations = 'output/GTEx/{chr}/{gene}/{gene}.associations',
        sumstats='output/UKBB_continuous/{phenotype}/{phenotype}.tsv.bgz',
        tabix_index='output/UKBB_continuous/{phenotype}/{phenotype}.tsv.bgz.tbi',
        v2r = 'output/GTEx/{chr}/{gene}/{gene}.snp2rsid'
    output:
        variant_report='output/UKBB_continuous/{phenotype}/{chr}/{gene}/{gene}.{phenotype}.z.variant_report0',
        coloc_report='output/UKBB_continuous/{phenotype}/{chr}/{gene}/{gene}.{phenotype}.z.coloc_report0',
        model='output/UKBB_continuous/{phenotype}/{chr}/{gene}/{gene}.{phenotype}.z.css0'
    run:
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

        gc = pd.read_csv('output/annotations/genes_hg19.bed', sep='\t')
        gc.loc[:, 'left'] = np.maximum(0, gc.start - 1e6)
        gc.loc[:, 'right'] = gc.end + 1e6

        gene2chr = gc.set_index('gene_id').chrom.to_dict()
        gene2left = gc.set_index('gene_id').left.to_dict()
        gene2right = gc.set_index('gene_id').right.to_dict()


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
        study = 'UKBB_continuous'
        phenotype = wildcards.phenotype

        # load summary stats
        gtex = load_gtex_associations(gene)
        rsid2variant_id = gtex.set_index('rsid').variant_id.to_dict()

        gwas = load_ukbb_gwas(phenotype, gene)

        # combine summary stat
        shared_variants = np.intersect1d(gtex.rsid.unique(), gwas.rsid.unique())
        df = pd.concat([gtex, gwas])
        df = df[df.rsid.isin(shared_variants)]#load genotype
        gtex_genotype = load_gtex_genotype(gene, use_rsid=True)
        gtex_genotype = gtex_genotype.loc[:,~gtex_genotype.columns.duplicated()]

        # load summary stats
        gtex = load_gtex_associations(gene)
        gwas = load_ukbb_gwas(phenotype, gene)

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

        mask = ~(np.any(np.isnan(B), 0) | np.any(np.isnan(S), 0))
        B = B.loc[:, mask]
        S = S.loc[:, mask]

        variants = B.columns.values
        study_ids = B.index.values

        print('{} GTEx variants'.format(gtex.rsid.unique().size))
        print('{} UKBB variants'.format(gwas.rsid.unique().size))
        print('{} intersecting, fully observed variants'.format(variants.size))

        K = 20

        init_args = {
            'LD': sample_ld(gtex_genotype.loc[:, variants]),
            'B': B.values / S.values,
            'S': np.ones_like(B),
            'K': K,
            'snp_ids': variants,
            'study_ids': study_ids,
            'tolerance': 1e-8
        }
        css = CSS(**init_args)
        css.prior_activity = np.ones(K) * 0.1
        css.weight_precision_b = np.ones_like(css.weight_precision_b) * 1


        print('fit model with imputed z-score')
        weight_ard_active_fit_procedure(css, max_iter=10, verbose=True)
        fit_all(css, max_iter=30, verbose=True)

        # save variant report
        table = make_table(css, gene, rsid2variant_id)
        table.to_csv(output.variant_report, sep='\t', index=False)
        
        ct = coloc_table(css, phenotype, gene=gene)
        ct.to_csv(output.coloc_report, sep='\t', index=False)

        css.save(output.model)

# UKBB SAIGE
ukbb_phenotypes = pd.read_csv('output/UKBB/UKBB_phenotypes.txt', sep ='\t', header=None)

rule download_ukbb_saige_sumstats:
    output:
        sumstats='output/UKBB/{phenotype}/{phenotype}.tsv.bgz'
    params:
        phecode=lambda wildcards: ukbb_phenotypes.set_index(1).loc[wildcards.phenotype].iloc[0]
    shell:
        "wget ftp://share.sph.umich.edu/UKBB_SAIGE_HRC//PheCode_{params.phecode}_SAIGE_MACge20.txt.vcf.gz -O {output}"

rule ukbb_saige_build_index:
    input:
        sumstats='output/UKBB/{phenotype}/{phenotype}.tsv.bgz',
    output:
        tabix_index='output/UKBB/{phenotype}/{phenotype}.tsv.bgz.tbi'
    shell:
        "tabix -s 1 -b 2 -e 2 -S 1 {input}"

