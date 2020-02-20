import itertools
import pickle
import numpy as np
import pandas as pd
import glob

configfile: "config/config.yaml"

# terminal rules
rule generate_figures:
    input:
        expand(
            "output/GTEx/gene_{gene}/summary.zscores.png", gene=config['chr22_genes']
        ),
        expand(
            "output/GTEx/gene_ENSG00000100078.3/regressed_genotype_cov/summary.zscores.png", gene=config['chr22_genes']
        ),
        expand(
            "output/GTEx/gene_ENSG00000100078.3/tissue_specific_cov/summary.zscores.png", gene=config['chr22_genes']
        )

rule all_tissue_pairs:
    input:
        expand(
            "output/simulation/single_causal_variant/pve_{pve}/ld_{linkage}/gene_{gene}/pairs_summary",
            pve=config["pves"], linkage=config["linkages"], gene=config["genes"]
        )

rule all_pairwise_pairs:
    input:
        expand(
            ("output/simulation/single_causal_variant/pve_{pve}/"
            "ld_{linkage}/gene_{gene}/pairwise_summary/pairs_table"),
            pve=config['pves'], linkage=config['linkages'],
            gene=config['genes']
        )

rule all_coloc:
    input:
        expand(
            "output/simulation/single_causal_variant/pve_{pve}/ld_{linkage}/gene_{gene}/coloc",
            pve=config["pves"], linkage=config["linkages"], gene=config["genes"]
        )

rule run_multiple_causal_variant_simulation:
    input:
        expand(
            "output/simulation/multiple_causal_variant/pve_{pve}/sparsity_0.2/"
            "gene_{gene}/pairs_summary",
            pve=config['pves'], gene=config['genes']
        ),
        expand(
            "output/simulation/multiple_causal_variant/pve_{pve}/sparsity_0.2/"
            "gene_{gene}/coloc",
            pve=config['pves'], gene=config['genes']
        ),
        expand(
            "output/simulation/multiple_causal_variant/pve_{pve}/sparsity_0.2/"
            "gene_{gene}/ecaviar",
            pve=config['pves'], gene=config['genes']
        ),
        expand(
            "output/simulation/multiple_causal_variant/pve_{pve}/sparsity_0.2/"
            "gene_{gene}/max_min_variance_summary",
            pve=config['pves'], gene=config['genes']
        )

rule run_single_causal_variant_simulation:
    input:
        expand(
            "output/simulation/single_causal_variant/"
            "pve_{pve}/ld_{linkage}/gene_{gene}/pairs_summary",
            pve=config["pves"], linkage=config["linkages"], gene=config["genes"]
        ),
        expand(
            "output/simulation/single_causal_variant/"
            "pve_{pve}/ld_{linkage}/gene_{gene}/coloc",
            pve=config["pves"], linkage=config["linkages"], gene=config["genes"]
        ),
        expand(
            "output/simulation/single_causal_variant/"
            "pve_{pve}/ld_{linkage}/gene_{gene}/ecaviar",
            pve=config["pves"], linkage=config["linkages"], gene=config["genes"]
        ),
        expand(
            "output/simulation/single_causal_variant/"
            "pve_{pve}/ld_{linkage}/gene_{gene}/max_min_variance_summary",
            pve=config["pves"], linkage=config["linkages"], gene=config["genes"]
        )

rule run_chr22_cafeh_summary:
    input:
        expand(
            "output/GTEx/gene_{gene}/model_summary", gene=config['chr22_genes']
        )

rule run_chr22_cafeh_summary_ts_cov:
    input:
        expand(
            "output/GTEx/gene_{gene}/tissue_specific_cov/model_summary", gene=config['chr22_genes']
        )

rule run_chr22_cafeh_summary_rg_cov:
    input:
        expand(
            "output/GTEx/gene_{gene}/regressed_genotype_cov/model_summary", gene=config['chr22_genes']
        )

rule run_chr22_cafeh_genotype:
    input:
        expand(
            "output/GTEx/gene_{gene}/model_genotype", gene=config['chr22_genes']
        )

# intermediate rules
rule get_cis_variants:
    output:
        "output/genotypes/{gene}_cis_variants"
    script:
        "workflow/scripts/get_cis_variants.py"

# simulation rules
rule simulate_single_causal_variant:
    input:
        "output/genotypes/{gene}_cis_variants"
    output:
        "output/simulation/single_causal_variant/pve_{pve}/ld_{linkage}/gene_{gene}/data"
    script:
        "workflow/scripts/single_causal_variant.py"

rule simulate_multiple_causal_variant:
    input:
        "output/genotypes/{gene}_cis_variants"
    output:
        "output/simulation/multiple_causal_variant/"
        "pve_{pve}/sparsity_{sparsity}/gene_{gene}/data"
    script:
        "workflow/scripts/multiple_causal_variant.py"

rule get_gtex_data:
    input:
        "output/genotypes/{gene}_cis_variants"
    output:
        "output/GTEx/gene_{gene}/data"
    wildcard_constraints:
        gene = "(?!\/)[^\/]+(?=\/)"
    script:
        "workflow/scripts/get_gtex_data.py"

rule get_tissue_specific_cov:
    input:
        "output/GTEx/gene_{gene}/data"
    output:
        "output/GTEx/gene_{gene}/tissue_specific_cov/data"
    params:
        alpha = 0.001
    script:
        "workflow/scripts/get_tissue_specific_cov.py"

rule get_regressed_genotype_cov:
    input:
        "output/GTEx/gene_{gene}/data"
    output:
        "output/GTEx/gene_{gene}/regressed_genotype_cov/data"
    script:
        "workflow/scripts/get_regressed_genotype_cov.py"

rule fit_summary_model:
    input:
        "output/{path}/data"
    output:
        "output/{path}/model_summary"
    script:
        "workflow/scripts/fit_cafeh_summary.py"

rule fit_genotype_model:
    input:
        "output/{path}/data"
    output:
        "output/{path}/model_genotype"
    script:
        "workflow/scripts/fit_cafeh_genotype.py"

rule fit_pairwise_summary_model:
    input:
        "output/simulation/single_causal_variant/pve_{pve}/ld_{linkage}/gene_{gene}/data"
    output:
        "output/simulation/{simulation}/"
        "{settings}/gene_{gene}/pairwise_summary/"
        "t1_{tissue1}_t2_{tissue2}_model_summary"
    wildcard_constraints:
        simulation = "(?!\/)[^\/]+(?=\/)"
    script:
        "workflow/scripts/fit_cafeh_summary.py"

rule fit_pairwise_genotype_model:
    input:
        "output/simulation/single_causal_variant/pve_{pve}/ld_{linkage}/gene_{gene}/data"
    output:
        ("output/simulation/single_causal_variant/pve_{pve}/"
        "ld_{linkage}/gene_{gene}/pairwise_genotype/"
        "t1_{tissue1}_t2_{tissue2}_model_genotype")
    script:
        "workflow/scripts/fit_cafeh_summary.py"

rule run_coloc:
    input:
        "output/simulation/{simulation}/{settings}/gene_{gene}/data"
    output:
        "output/simulation/{simulation}/{settings}/gene_{gene}/coloc"
    wildcard_constraints:
        simulation = "(?!\/)[^\/]+(?=\/)"
    script:
        "workflow/scripts/run_coloc.R"

rule format_caviar_data_ld:
    input:
        'output/{path}/data'
    output:
        ld_matrix='output/{path}/caviar/data.ld'
    run:
        data = pickle.load(open(input[0], 'rb'))
        np.savetxt(fname=output.ld_matrix, X=data['LD'], delimiter='\t')

rule format_caviar_data_zscore:
    input:
        'output/{path}/data'
    output:
        z_scores='output/{path}/caviar/data.z{tissue}'
    run:
        data = pickle.load(open(input[0], 'rb'))
        zscores = pd.DataFrame(data['zscores'][int(wildcards.tissue)])
        zscores.to_csv(output.z_scores, sep='\t', header=None)

rule run_caviar:
    input:
        ld_matrix='output/{path}/caviar/data.ld',
        z_scores='output/{path}/caviar/data.z{tissue}'
    output:
        'output/{path}/caviar/caviar_t{tissue}.log',
        'output/{path}/caviar/caviar_t{tissue}_post',
        'output/{path}/caviar/caviar_t{tissue}_set'
    shell:
        "workflow/bin/caviar/CAVIAR "
        "-o output/{wildcards.path}/caviar/caviar_t{wildcards.tissue} "
        "-l {input.ld_matrix} "
        "-z {input.z_scores} "
        "-c 2"

rule make_ecaviar_table:
    input:
        data = 'output/{path}/data',
        caviar_posteriors = expand(
            'output/{path}/caviar/caviar_t{tissue}_post',
            path='{path}', tissue=np.arange(7)
        )
        # todo get this to adapt to the number of tissues
    output:
        'output/{path}/ecaviar'
    script:
        'workflow/scripts/make_ecaviar_table.py'

rule make_max_min_variance_table:
    input:
        data = 'output/{path}/data',
        model = 'output/{path}/model_summary'
    output:
        'output/{path}/max_min_variance_summary'
    script:
        'workflow/scripts/make_variance_table.py'

# stat gathering rules
rule make_tissue_pair_components_table:
    input:
        data_path = \
            "output/simulation/{simulation}/{settings}/gene_{gene}/data",
        genotype_model_path = \
            "output/simulation/{simulation}/{settings}/gene_{gene}/model_genotype",
        summary_model_path = \
            "output/simulation/{simulation}/{settings}/gene_{gene}/model_summary"
    output:
        genotype_output = \
            "output/simulation/{simulation}/{settings}/gene_{gene}/pairs_genotype",
        summary_output = \
            "output/simulation/{simulation}/{settings}/gene_{gene}/pairs_summary"
    wildcard_constraints:
        simulation = "(?!\/)[^\/]+(?=\/)"
    script:
        "workflow/scripts/make_tissue_pair_components_table.py"

tissue_pairs = [x for x in itertools.combinations(np.arange(7), 2)]
rule make_pairwise_pair_components_table:
    """
    same rule as above except for the pairwise outputs
    """
    input:
        data_path = "output/simulation/single_causal_variant/pve_{pve}/ld_{linkage}/gene_{gene}/data",
        summary_model_paths = expand(
            "output/simulation/single_causal_variant/pve_{pve}/ld_{linkage}/gene_{gene}/pairwise_summary/t1_{tissue1}_t2_{tissue2}_model_summary",
            tissue1=[x[0] for x in tissue_pairs],
            tissue2=[x[1] for x in tissue_pairs],
            pve='{pve}', linkage='{linkage}', gene='{gene}'
        )
    output:
        "output/simulation/single_causal_variant/pve_{pve}/ld_{linkage}/gene_{gene}/pairwise_summary/pairs_table"
    script:
        "workflow/scripts/make_pairwise_pair_components_table.py"

rule make_cafeh_plots:
    input:
        data_path = 'output/{path}/data',
        model_path = 'output/{path}/model_summary'
    output:
        component_plot_path = report('output/{path}/summary.components.png'),
        zscore_plot_path = report('output/{path}/summary.zscores.png')
    script:
        'workflow/scripts/cafeh_plots.py'


