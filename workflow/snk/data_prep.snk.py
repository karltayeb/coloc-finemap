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


rule grep_associations:
    output:
        'output/GTEx/gene_{gene}/associations/{tissue}.associations'
    run:
        'grep {gene} \
        /work-zfs/abattle4/lab_data/GTEx_v8/ciseQTL/GTEx_Analysis_v8_eQTL_all_associations{tissue}.allpairs.txt \
        > {output}'

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
