# intermediate rules
rule get_cis_variants:
    output:
        "output/genotypes/{gene}_cis_variants"
    script:
        "workflow/scripts/get_cis_variants.py"

rule get_gtex_data:
    input:
        "output/genotypes/{gene}_cis_variants"
    output:
        "output/GTEx/gene_{gene}/data"
    wildcard_constraints:
        gene = "(?!\/)[^\/]+(?=\/)"
    script:
        "workflow/scripts/get_gtex_data.py"


import glob
tissues = [x.split('.')[0].split('/')[-1] for x in glob.glob(
    '/work-zfs/abattle4/lab_data/GTEx_v8/ciseQTL/GTEx_Analysis_v8_eQTL_all_associations/*allpairs.txt')]
rule build_indices:
    input:
        expand(
            'output/GTEx/index/{tissue}.association.index', tissue=tissues
        )

rule get_gtex_associations:
    input:
        expand(
            'output/GTEx/index/{tissue}.association.index', tissue=tissues
        )
    output:
        'output/GTEx/gene_{gene}/{gene}.associations'
    script:
        'workflow/scripts/get_gtex_associations.py'

rule get_gtex_ld:
    input:
        'output/GTEx/gene_{gene}/{gene}.associations'
    output:
        'output/GTEx/gene_{gene}/{gene}.ld'
    shell:
        'plink --bfile /work-zfs/abattle4/marios/GTEx_v8/coloc/GTEx_all_genotypes'
        ' --chr chr1 --from-bp 13550 --to-bp 1631911  --maf 0.01 --r square'

rule build_gene_seek_index:
    output:
        'output/GTEx/index/{tissue}.association.index'
    run:
        from collections import defaultdict
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

rule get_tissue_specific_cov:
    input:
        "output/{path}/data"
    output:
        "output/{path}/tissue_specific_cov/data"
    params:
        alpha = None
    script:
        "workflow/scripts/get_tissue_specific_cov.py"

rule get_regressed_genotype_cov:
    input:
        "output/{path}/data"
    output:
        "output/{path}/regressed_genotype_cov/data"
    script:
        "workflow/scripts/get_regressed_genotype_cov.py"

rule get_global_tissue_cov:
    input:
        "output/{path}/data"
    output:
        "output/{path}/global_regularized_cov/data"
    script:
        "workflow/scripts/get_global_regularized_cov.py"

rule make_maf_tss_table:
    input:
        'maf/{part}'
    output:
        'output/enrichment/GTEx_maf_tss/GTEx_maf_tss.{part}'
    script:
        'workflow/scripts/make_maf_tss_table.py'
