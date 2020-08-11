include: 'workflow/snk/data_prep.snk.py'


rule fit_gss20:
    input:
        "output/{path}/genotype.data"
    output:
        "output/{path}/cafeh.k20.gss"
    params:
        k=20
    script:
        "../../workflow/scripts/fit_gss.py"

rule fit_gss_k_pi:
    input:
        # TODO get all inputs
        "output/{path}/{gene}.raw",
        "output/{path}/{gene}.1kG.raw",
        "output/{path}/{gene}.associations"
    output:
        "output/{path}/{gene}.k20.pi01.gss"
    params:
        k=20,
        pi0 = 0.01
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

rule fit_cad_gtex_cafeh:
    input:
        genotype_gtex = 'output/GTEx/{chr}/{gene}/{gene}.raw',
        genotype_1kG = 'output/GTEx/{chr}/{gene}/{gene}.1kG.raw',
        associations = 'output/GTEx/{chr}/{gene}/{gene}.associations',
        v2r = 'output/GTEx/{chr}/{gene}/{gene}.associations'
    output:
        'output/CAD/{chr}/{gene}/{gene}.cad_gtex.css'
    script:
        '../../workflow/scripts/cad_cafeh_ss.py'


rule fit_gtex_no_spike_and_slab:
    input:
        genotype = rule.get_gtex_genotype.output.genotype,
        expression = rule.get_gtex_genotype.output.expression,
        rsid_map = rule.snpid2rsid.rsid_map
        variants = None,
    output:
        "output/{path}/{gene}.k20.pi01.gss"
    params:
        k=20,
        pi0 = 0.01
    script:
        "../../workflow/scripts/fit_gss.py"

rule fit_genotype_model_general:
    input:
        genotype = rule.get_gtex_genotype.output.genotype,
        expression = rule.get_gtex_genotype.output.expression,
        rsid_map = rule.snpid2rsid.rsid_map
    output:
        genotype_model = '',
        gss = ''
        "output/{path}/genotype.standardized.k40.model"
    params:
        k=40
    script:
        "../../workflow/scripts/fit_cafeh_genotype.py"
