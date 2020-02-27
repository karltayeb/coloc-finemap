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
        "output/{path}/data"
    output:
        "{path}/pairwise_summary/"
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
            "output/{path}/data",
        summary_model_path = \
            "output/simulation/{path}/model_summary"
    output:
        summary_output = \
            "output/{path}/pairs_summary"
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