import itertools
configfile: "config/config.yaml"

rule all:
    input:
        expand(
            "output/simulation/single_causal_variant/pve_{pve}/ld_{linkage}/gene_{gene}/pairs_summary",
            pve=config["pves"], linkage=config["linkages"], gene=config["genes"]
        ),
        expand(
            "output/simulation/single_causal_variant/pve_{pve}/ld_{linkage}/gene_{gene}/coloc",
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

rule fit_summary_model:
    input:
        "output/simulation/single_causal_variant/pve_{pve}/ld_{linkage}/gene_{gene}/data"
    output:
        "output/simulation/single_causal_variant/pve_{pve}/ld_{linkage}/gene_{gene}/model_summary"
    script:
        "workflow/scripts/fit_cosie_summary.py"

rule fit_genotype_model:
    input:
        "output/simulation/single_causal_variant/pve_{pve}/ld_{linkage}/gene_{gene}/data"
    output:
        "output/simulation/single_causal_variant/pve_{pve}/ld_{linkage}/gene_{gene}/model_genotype"
    script:
        "workflow/scripts/fit_cosie_genotype.py"

rule fit_pairwise_summary_model:
    input:
        "output/simulation/single_causal_variant/pve_{pve}/ld_{linkage}/gene_{gene}/data"
    output:
        ("output/simulation/single_causal_variant/pve_{pve}/"
        "ld_{linkage}/gene_{gene}/pairwise_summary/"
        "t1_{tissue1}_t2_{tissue2}_model_summary")
    script:
        "workflow/scripts/fit_cosie_summary.py"

rule fit_pairwise_genotype_model:
    input:
        "output/simulation/single_causal_variant/pve_{pve}/ld_{linkage}/gene_{gene}/data"
    output:
        ("output/simulation/single_causal_variant/pve_{pve}/"
        "ld_{linkage}/gene_{gene}/pairwise_genotype/"
        "t1_{tissue1}_t2_{tissue2}_model_genotype")
    script:
        "workflow/scripts/fit_cosie_summary.py"

rule all_pairwise_pairs:
    input:
        expand(
            ("output/simulation/single_causal_variant/pve_{pve}/"
            "ld_{linkage}/gene_{gene}/pairwise_summary/pairs_summary"),
            pve=config['pves'], linkage=config['linkages'],
            gene=config['genes']
        )

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

rule make_pairwise_pair_components_table:
    """
    same rule as above except for the pairwise outputs
    """

    input:
        data_path = "output/simulation/single_causal_variant/pve_{pve}/ld_{linkage}/gene_{gene}/data",
        summary_model_paths = expand(
            "output/simulation/single_causal_variant/pve_{pve}/ld_{linkage}/gene_{gene}/pairwise_summary/t1_{tissue1}_t2_{tissue2}_model_summary",
            tissue1=[0, 1, 4], tissue2=[2, 5], pve='{pve}', linkage='{linkage}', gene='{gene}'
        )
    output:
        "output/simulation/single_causal_variant/pve_{pve}/ld_{linkage}/gene_{gene}/pairwise_summary/pairs_table"
    script:
        "workflow/scripts/make_pairwise_pair_components_table.py"

