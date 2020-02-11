configfile: "config/config.yaml"

rule all:
    input:
        expand(
            "output/simulation/single_causal_variant/pve_{pve}/ld_{linkage}/gene_{gene}/pairs_summary",
            pve=config["pves"], linkage=config["linkages"], gene=config["genes"]
        )

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

rule fit_pairwise_cosie_summary:
    input:
        "output/simulation/single_causal_variant/pve_{pve}/ld_{linkage}/gene_{gene}/data"
    output:
        "output/simulation/single_causal_variant/pve_{pve}/"
        "ld_{linkage}/gene_{gene}/model_summary_pairwise/"
        "t1_{tissue1}_t2_{tissue2}_model_summary"
    scripts:
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

rule make_tissue_pair_components_table:
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
            "output/simulation/single_causal_variant/pve_{pve}/ld_{linkage}/gene_{gene}/pairs_summary"
    script:
        "workflow/scripts/make_tissue_pair_components_table.py"

