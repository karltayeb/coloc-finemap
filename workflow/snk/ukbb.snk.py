# eQTL gene

rule download_eqtlgen:
    output:
        'output/eQTLGEN/eQTLGEN.txt.gz'
    run:
     shell('wget https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/cis-eQTLs_full_20180905.txt.gz -O {output}')
     # zcat eQTLGEN.txt.gz | awk '{print > "chr"$3"_eQTLGEN.tsv"}'

rule eqtlgen_build_index:
    input:
        'output/eQTLGEN/eQTLGEN.txt.gz'
    output:
        sumstats='output/eQTLGEN/eQTLGEN/eQTLGEN.tsv.bgz',
        tabix_index='output/eQTLGEN/eQTLGEN/eQTLGEN.tsv.bgz.tbi'
    run:
        shell("paste <(zcat {input.var}) <(zcat {input.sumstats_raw}) | bgzip > {output.sumstats}")
        shell("tabix -s 2 -b 3 -e 3 -S 1 {output.sumstats}")

rule run_eqtlgen_x_gtex:
    input:
        'output/eQTLGEN/{chr}_eQTLGEN.tsv'
    output:
        attempted_genes_path = 'output/eQTLGEN/eQTLGEN/{chr}/{chr}_attempted_genes.txt',
        variant_report_path = 'output/eQTLGEN/eQTLGEN/{chr}/{chr}_variant_report.txt',
        active_path = 'output/eQTLGEN/eQTLGEN/{chr}/{chr}_p_active.txt',
        error_path = 'output/eQTLGEN/eQTLGEN/{chr}/{chr}_error.txt',
    notebook:
        "notebooks/final_analysis/eQTLgen.ipynb"

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
        genes='output/{study}/{phenotype}/{phenotype}.genes.txt'
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

rule gwas_generate_requests:
    input:
        genes='output/{study}/{phenotype}/{phenotype}.genes.txt'
    output:
        request='output/{study}/{phenotype}/{phenotype}.{request}.requests.txt'
    run:
        from glob import glob
        import pandas as pd
        from tqdm import tqdm

        request = wildcards.request

        template = 'output/{study}/{phenotype}/{chr}/{gene}/{gene}.{phenotype}.{request}'

        df = pd.read_csv(input.genes, sep='\t', index_col=0)
        df = df[~df.gene_id.duplicated()]

        requests = []
        for _, row in df.iterrows():
            requests.append(template.format(
                study=row.study,
                phenotype=row.phenotype,
                gene=row.gene_id,
                chr=row.chr,
                request=request
            ))
        with open(output.request, 'w') as f:
            [print(r, file=f) for r in requests];

    
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


# UKBB SAIGE
ukbb_phenotypes = pd.read_csv(config['UKBB_phenotypes'], sep ='\t', header=None)

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


rule grasp_build_index:
    input:
        sumstats='output/GRASP/{phenotype}/{phenotype}.tsv.bgz',
    output:
        tabix_index='output/GRASP/{phenotype}/{phenotype}.tsv.bgz.tbi'
    shell:
        "tabix -s 1 -b 2 -e 2 -S 1 {input}"
