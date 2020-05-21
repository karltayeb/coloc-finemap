rule fit_summary_model:
    input:
        "output/{path}/data"
    output:
        "output/{path}/model_summary"
    script:
        "../../workflow/scripts/fit_cafeh_summary.py"

rule fit_cafeh2:
    input:
        "output/{path}/data"
    output:
        "output/{path}/cafeh2.model"
    script:
        "../../workflow/scripts/fit_cafeh2.py"

rule fit_study_model:
    input:
        "output/{path}/data"
    output:
        "output/{path}/model_study"
    script:
        "../../workflow/scripts/fit_susie_summary.py"

rule fit_genotype_model:
    input:
        "output/{path}/genotype_data"
    output:
        "output/{path}/model_genotype"
    script:
        "../../workflow/scripts/fit_cafeh_genotype.py"

rule fit_genotype_model2:
    input:
        "output/{path}/genotype.data"
    output:
        "output/{path}/genotype.model"
    params:
        k=20
    script:
        "../../workflow/scripts/fit_cafeh_genotype.py"

rule fit_genotype_model40:
    input:
        "output/{path}/genotype.data"
    output:
        "output/{path}/genotype.k40.model"
    params:
        k=40
    script:
        "../../workflow/scripts/fit_cafeh_genotype.py"

rule fit_gss20:
    input:
        "output/{path}/genotype.data"
    output:
        "output/{path}/cafeh.k20.gss"
    params:
        k=20
    script:
        "../../workflow/scripts/fit_gss.py"

rule fit_standardized_gss20:
    input:
        expression = 'output/GTEx/{chr}/{gene}/{gene}.expression',
        genotype = 'output/GTEx/{chr}/{gene}/{gene}.raw',
        model="output/GTEx/{chr}/{gene}/genotype.standardized.k20.model"
    output:
        "output/GTEx/{chr}/{gene}/gss.standardized.k20.pi01.model"
    params:
        k=20,
        pi=0.01
    script:
        "../../workflow/scripts/genotype_model_to_gss.py"

rule fit_genotype_model20:
    input:
        "output/{path}/genotype.data"
    output:
        "output/{path}/genotype.k20.model"
    params:
        k=20
    script:
        "../../workflow/scripts/fit_cafeh_genotype.py"


rule fit_standardized_genotype_model20:
    input:
        "output/{path}/genotype.standardized.data"
    output:
        "output/{path}/genotype.standardized.k20.model"
    params:
        k=20
    script:
        "../../workflow/scripts/fit_cafeh_genotype.py"

rule fit_standardized_genotype_model40:
    input:
        "output/{path}/genotype.standardized.data"
    output:
        "output/{path}/genotype.standardized.k40.model"
    params:
        k=40
    script:
        "../../workflow/scripts/fit_cafeh_genotype.py"


common = pd.read_csv('../../output/GTEx/gss_css_common.txt', header=None).values.flatten()
rule run_gss_css_gtex:
    input:
        list(common[:100])

rule gss_css_gtex:
    input:
        expression = 'output/GTEx/{chr}/{gene}/{gene}.expression',
        genotype = 'output/GTEx/{chr}/{gene}/{gene}.raw',
        associations = 'output/GTEx/{chr}/{gene}/{gene}.associations',
    output:
        gss='output/GTEx/{chr}/{gene}/common.gss',
        css='output/GTEx/{chr}/{gene}/common.recompute.css',
        gtex_css='output/GTEx/{chr}/{gene}/common.gtex.css',

    script:
        "../../workflow/scripts/gss_css_gtex.py"


rule fit_pairwise_summary_model:
    input:
        "output/{path}/data"
    output:
        "{path}/pairwise_summary/"
        "t1_{tissue1}_t2_{tissue2}_model_summary"
    wildcard_constraints:
        simulation = "(?!\/)[^\/]+(?=\/)"
    script:
        "../../workflow/scripts/fit_cafeh_summary.py"

rule fit_pairwise_genotype_model:
    input:
        "output/simulation/single_causal_variant/pve_{pve}/ld_{linkage}/gene_{gene}/data"
    output:
        ("output/simulation/single_causal_variant/pve_{pve}/"
        "ld_{linkage}/gene_{gene}/pairwise_genotype/"
        "t1_{tissue1}_t2_{tissue2}_model_genotype")
    script:
        "../../workflow/scripts/fit_cafeh_summary.py"

rule make_max_min_variance_table:
    input:
        data = 'output/{path}/data',
        model = 'output/{path}/model_summary'
    output:
        'output/{path}/max_min_variance_summary'
    script:
        '../../workflow/scripts/make_variance_table.py'

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
        "../../workflow/scripts/make_tissue_pair_components_table.py"

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
        "../../workflow/scripts/make_pairwise_pair_components_table.py"

rule make_cafeh_plots:
    input:
        data_path = 'output/{path}/data',
        model_path = 'output/{path}/model_summary'
    output:
        component_plot_path = report('output/{path}/summary.components.png'),
        zscore_plot_path = report('output/{path}/summary.zscores.png')
    script:
        '../../workflow/scripts/cafeh_plots.py'

rule report_genotype:
    input:
        expression = 'output/GTEx/{chr}/{gene}/{gene}.expression',
        genotype = 'output/GTEx/{chr}/{gene}/{gene}.raw',
        model = 'output/GTEx/{chr}/{gene}/genotype.model'
    output:
        scores = "output/GTEx/{chr}/{gene}/genotype.scores",
        csets = "output/GTEx/{chr}/{gene}/genotype.csets"
    script:
        '../../workflow/scripts/report_genotype.py'


rule report_genotypek20:
    input:
        expression = 'output/GTEx/{chr}/{gene}/{gene}.expression',
        genotype = 'output/GTEx/{chr}/{gene}/{gene}.raw',
        model = 'output/GTEx/{chr}/{gene}/genotype.k20.model'
    output:
        scores = "output/GTEx/{chr}/{gene}/genotype.k20.scores",
        csets = "output/GTEx/{chr}/{gene}/genotype.k20.csets"
    script:
        '../../workflow/scripts/report_genotype.py'

rule report_standardized_genotypek20:
    input:
        expression = 'output/GTEx/{chr}/{gene}/{gene}.expression',
        genotype = 'output/GTEx/{chr}/{gene}/{gene}.raw',
        model = 'output/GTEx/{chr}/{gene}/genotype.standardized.k20.model'
    output:
        scores = "output/GTEx/{chr}/{gene}/genotype.standardized.k20.scores",
        csets = "output/GTEx/{chr}/{gene}/genotype.standardized.k20.csets"
    script:
        '../../workflow/scripts/report_genotype.py'

