configfile: "config/config.yaml"

rule get_cis_variants:
    output:
        "output/genotypes/{gene}_cis_variants"
    script:
        "workflow/scripts/get_cis_variants.py"

rule simulate_single_causal_variant:
    input:
        "output/genotypes/{gene}_cis_variants"
    output:
        "output/simulation/single_causal_variant/pve_{pve}/ld_{linkage}/gene_{gene}/data"
    script:
        "workflow/scripts/single_causal_variant.py"

rule fit_cosie_summary:
    input:
        "output/simulation/single_causal_variant/pve_{pve}/ld_{linkage}/gene_{gene}/data"
    output:
        "output/simulation/single_causal_variant/pve_{pve}/ld_{linkage}/gene_{gene}/model_summary"
    script:
        "workflow/scripts/fit_cosie_summary.py"

rule fit_cosie_genotype:
    input:
        "output/simulation/single_causal_variant/pve_{pve}/ld_{linkage}/gene_{gene}/data"
    output:
        "output/simulation/single_causal_variant/pve_{pve}/ld_{linkage}/gene_{gene}/model_genotype"
    script:
        "workflow/scripts/fit_cosie_genotype.py"

rule run_coloc:
    input:
        "output/simulation/single_causal_variant/pve_{pve}/ld_{linkage}/gene_{gene}/data"
    output:
        "output/simulation/single_causal_variant/pve_{pve}/ld_{linkage}/gene_{gene}/coloc"
    script:
        "workflow/scripts/run_coloc.R"

rule generate_cosie_results:
    input:
        data_path = \
            "output/simulation/single_causal_variant/pve_{pve}/ld_{linkage}/gene_{gene}/data",
        genotype_model_path = \
            "output/simulation/single_causal_variant/pve_{pve}/ld_{linkage}/gene_{gene}/model_genotype",
        summary_model_path = \
            "output/simulation/single_causal_variant/pve_{pve}/ld_{linkage}/gene_{gene}/model_summary"
    output:
        genotype_output = \
            "output/simulation/single_causal_variant/pve_{pve}/ld_{linkage}/gene_{gene}/pairs_genotype",
        summary_output = \
            "output/simulation/single_causal_variant/pve_{pve}/ld_{linkage}/gene_{gene}/pairs_genotype",




