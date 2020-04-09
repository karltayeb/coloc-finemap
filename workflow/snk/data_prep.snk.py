import numpy as np
import json
import glob
import pandas as pd
from collections import defaultdict

gencode = pd.read_csv(
    'output/GTEx/protein_coding_autosomal_egenes.txt', sep='\t', index_col=0)

rule get_cis_variants:
    output:
        "output/genotypes/{gene}_cis_variants"
    script:
        "workflow/scripts/get_cis_variants.py"

tissues = [x.split('.')[0].split('/')[-1] for x in glob.glob(
    '/work-zfs/abattle4/lab_data/GTEx_v8/ciseQTL/GTEx_Analysis_v8_eQTL_all_associations/*allpairs.txt')]
rule build_indices:
    input:
        expand(
            'output/GTEx/index/{tissue}.association.index', tissue=tissues
        )


#####
# GATHER DATA
#####
rule get_gtex_associations:
    input:
        expand(
            'output/GTEx/index/{tissue}.association.index', tissue=tissues
        )
    output:
        'output/GTEx/gene_{gene}/{gene}.associations'
    script:
        '../../workflow/scripts/get_gtex_associations.py'

rule get_gtex_expression:
    input:
        expand(
            'output/GTEx/index/{tissue}.association.index', tissue=tissues
        )
    output:
        'output/GTEx/gene_{gene}/{gene}.expression'
    script:
        '../../workflow/scripts/get_gtex_expression.py'

rule get_gtex_ld:
    input:
        associations = 'output/GTEx/gene_{gene}/{gene}.associations'
    params:
        chrom = lambda wildcards: gencode.loc[wildcards.gene].chromosome,
        from_bp = lambda wildcards: gencode.loc[wildcards.gene].tss - 500000,
        to_bp = lambda wildcards: gencode.loc[wildcards.gene].tss + 500000
    output:
        'output/GTEx/gene_{gene}/{gene}.ld',
        'output/GTEx/gene_{gene}/{gene}.snplist',
        'output/GTEx/gene_{gene}/{gene}.raw'
    shell:
        'plink --bfile /work-zfs/abattle4/marios/GTEx_v8/coloc/GTEx_all_genotypes'
        ' --chr {params.chrom} --from-bp {params.from_bp} --to-bp {params.to_bp}  --maf 0.01 --r square'
        ' --out output/GTEx/gene_{wildcards.gene}/{wildcards.gene} --write-snplist --recodeA'

rule get_gtex_data:
    input:
        associations = 'output/GTEx/gene_{gene}/{gene}.associations',
        ld = 'output/GTEx/gene_{gene}/{gene}.ld',
        snps = 'output/GTEx/gene_{gene}/{gene}.snplist'
    output:
        "output/GTEx/gene_{gene}/data"
    wildcard_constraints:
        gene = "(?!\/)[^\/]+(?=\/)"
    script:
        "../../workflow/scripts/get_gtex_data.py"

rule get_gtex_genotype_data:
    input:
        expression = 'output/GTEx/gene_{gene}/{gene}.expression',
        genotype = 'output/GTEx/gene_{gene}/{gene}.raw',
    output:
        "output/GTEx/gene_{gene}/genotype_data"
    wildcard_constraints:
        gene = "(?!\/)[^\/]+(?=\/)"
    script:
        "../../workflow/scripts/get_gtex_genotype_model_data.py"


# GATHER DATA
#####
rule get_gtex_associations2:
    input:
        expand(
            'output/GTEx/index/{tissue}.association.index', tissue=tissues
        )
    output:
        'output/GTEx/{chr}/{gene}/{gene}.zscores',
        'output/GTEx/{chr}/{gene}/{gene}.betas',
        'output/GTEx/{chr}/{gene}/{gene}.standard_errors'
    script:
        '../../workflow/scripts/get_gtex_associations.py'

rule get_gtex_expression2:
    input:
        expand(
            'output/GTEx/index/{tissue}.association.index', tissue=tissues
        )
    output:
        'output/GTEx/{chr}/{gene}/{gene}.expression'
    script:
        '../../workflow/scripts/get_gtex_expression.py'

rule get_gtex_genotype2:
    params:
        chrom = lambda wildcards: gencode.loc[wildcards.gene].chromosome,
        from_bp = lambda wildcards: np.maximum(0, gencode.loc[wildcards.gene].tss - 1000000),
        to_bp = lambda wildcards: gencode.loc[wildcards.gene].tss + 1000000
    output:
        'output/GTEx/{chrom}/{gene}/{gene}.snplist',
        'output/GTEx/{chrom}/{gene}/{gene}.raw'
    shell:
        'plink --bfile /work-zfs/abattle4/marios/GTEx_v8/coloc/GTEx_all_genotypes'
        ' --chr {params.chrom} --from-bp {params.from_bp} --to-bp {params.to_bp}  --maf 0.01  --geno 0.1'
        ' --out output/GTEx/{wildcards.chrom}/{wildcards.gene}/{wildcards.gene} --write-snplist --recodeA'

rule get_gtex_ld2:
    params:
        chrom = lambda wildcards: gencode.loc[wildcards.gene].chromosome,
        from_bp = lambda wildcards: np.maximum(0, gencode.loc[wildcards.gene].tss - 1000000),
        to_bp = lambda wildcards: gencode.loc[wildcards.gene].tss + 1000000
    output:
        'output/GTEx/{params.chrom}/{gene}/{gene}.ld',
        'output/GTEx/{params.chrom}/{gene}/{gene}.snplist'
    shell:
        'plink --bfile /work-zfs/abattle4/marios/GTEx_v8/coloc/GTEx_all_genotypes'
        ' --chr {params.chrom} --from-bp {params.from_bp} --to-bp {params.to_bp}  --maf 0.01  --geno 0.1 --r square'
        ' --out output/GTEx/{wildcards.chrom}/{wildcards.gene}/{wildcards.gene} --write-snplist'

rule get_gtex_data2:
    input:
        associations = 'output/GTEx/{chr}/{gene}/{gene}.associations',
        ld = 'output/GTEx/{chr}/{gene}/{gene}.ld',
        snps = 'output/GTEx/{chr}/{gene}/{gene}.snplist'
    output:
        "output/GTEx/{chr}/{gene}/data"
    wildcard_constraints:
        gene = "(?!\/)[^\/]+(?=\/)"
    script:
        "../../workflow/scripts/get_gtex_data.py"

rule get_gtex_genotype_data2:
    input:
        expression = 'output/GTEx/{chr}/{gene}/{gene}.expression',
        genotype = 'output/GTEx/{chr}/{gene}/{gene}.raw',
    output:
        temp("output/GTEx/{chr}/{gene}/genotype.data")
    wildcard_constraints:
        gene = "(?!\/)[^\/]+(?=\/)"
    script:
        "../../workflow/scripts/get_gtex_genotype_model_data.py"


rule build_gene_seek_index:
    output:
        'output/GTEx/index/{tissue}.association.index'
    run:
        print('building index for {}'.format(wildcards.tissue))
        gene_start = {}
        with open('/work-zfs/abattle4/lab_data/GTEx_v8/ciseQTL/'
                  'GTEx_Analysis_v8_eQTL_all_associations/{}.allpairs.txt'.format(wildcards.tissue), 'r') as f:
            i = 0
            last_gene = ''
            for line in f:
                gene = line.split('\t')[0]

                if gene != last_gene:
                    gene_start[gene] = i
                i += len(line)
                last_gene = gene
        json.dump(gene_start, open(output[0], 'w'))

rule build_gene_tissue_map:
    input:
        expand(
            'output/GTEx/index/{tissue}.association.index', tissue=tissues
        )
    output:
        'output/GTEx/gene_tissue_map.json'
    run:
        association_indices = {
            x.split('/')[-1].split('.')[0]: json.load(open(x, 'r'))
            for x in input
        }

        gene_tissue_map = {}
        for tissue in association_indices:
            for gene in association_indices[tissue]:
                if gene in gene_tissue_map:
                    gene_tissue_map[gene].append(tissue)
                else:
                    gene_tissue_map[gene] = [tissue]
        json.dump(gene_tissue_map, open(output[0], 'w'))

rule make_bins:
    input:
        'output/enrichment/GTEx_maf_tss/GTEx_maf_tss.{suffix}'
    output:
        'output/enrichment/GTEx_maf_tss_binned/bins.{suffix}'
    run:
        maf_bins = np.linspace(0, 1, 51)
        tss_bins = np.linspace(-500000, 500000, 51)

        bins = {}
        for m in range(len(maf_bins)):
            bins[m] = {}
            for t in range(len(tss_bins)):
                bins[m][t] = defaultdict(list)

        with open(input[0], 'r') as f:
            for line in f:
                chromosome, variant, gene, maf, dtss = line.split('\t')
                if float(maf) > 0.01:
                    maf_bin = np.digitize(float(maf), maf_bins)
                    tss_bin = np.digitize(float(dtss), tss_bins)
                    bins[maf_bin][tss_bin][variant].append(gene)
        json.dump(bins, open(output[0], 'w'))

rule make_gene_variant_lookup:
    input:
        'output/enrichment/GTEx_maf_tss/GTEx_maf_tss.{suffix}'
    output:
        expand('output/enrichment/GTEx_maf_tss_lookup/chr{chr}/chr{chr}.lookup.{suffix}',
            chr=list(range(1, 23)), suffix='{suffix}')
    run:
        lookup = {'chr{}'.format(x): defaultdict(dict) for x in range(1, 23)}
        with open(input[0], 'r') as f:
            for line in f:
                chromosome, variant, gene, maf, dtss = line.strip().split('\t')
                if float(maf) > 0.01:
                    lookup[chromosome][variant][gene] = int(dtss)

        for out_path in output:
            chromosome = out_path.split('/')[-1].split('.')[0]
            json.dump(lookup[chromosome], open(out_path, 'w'))

rule get_tissue_specific_cov:
    input:
        associations = 'output/GTEx/gene_{gene}/{gene}.associations',
        ld = 'output/GTEx/gene_{gene}/{gene}.ld',
        snps = 'output/GTEx/gene_{gene}/{gene}.snplist'
    output:
        "output/GTEx/gene_{gene}/tissue_specific_cov/data"
    params:
        alpha = None
    script:
        "../../workflow/scripts/get_tissue_specific_cov.py"

rule get_regressed_genotype_cov:
    input:
        "output/{path}/data"
    output:
        "output/{path}/regressed_genotype_cov/data"
    script:
        "../../workflow/scripts/get_regressed_genotype_cov.py"

rule get_global_tissue_cov:
    input:
        "output/{path}/data"
    output:
        "output/{path}/global_regularized_cov/data"
    script:
        "../../workflow/scripts/get_global_regularized_cov.py"

rule make_maf_tss_table:
    input:
        'maf/{part}'
    output:
        'output/enrichment/GTEx_maf_tss/GTEx_maf_tss.{part}'
    script:
        '../../workflow/scripts/make_maf_tss_table.py'
