rule fit_genotype_model:
    input:
        genotype = 'output/GTEx/{chr}/{gene}/{gene}.raw',
        expression = 'output/GTEx/{chr}/{gene}/{gene}.expression',
        snp2rsid = 'output/GTEx/{chr}/{gene}/{gene}.snp2rsid'
    output:
        model = 'output/GTEx/{chr}/{gene}/{gene}.cafeh_genotype'
    params:
        K = 20,
        p0k = 1.0
    run:
        from cafeh.independent_model_ss import CAFEHG
        from cafeh.fitting import forward_fit_procedure

        from utils.misc import load_gtex_genotype, load_gtex_expression
        genotype = load_gtex_genotype(wildcards.gene)
        X = np.nan_to_num(genotype.values - np.nanmean(genotype.values, 0)).T
        expression = load_gtex_expression(wildcards.gene)
        Y = expression.values
        covariates = pd.read_csv(
            '/work-zfs/abattle4/karl/cosie_analysis/output/GTEx/covariates.csv', sep='\t', index_col=[0, 1])
        covariates = covariates.loc[expression.index.values].loc[:, genotype.index.values]
        model = CAFEHG(
            X=X, Y=Y, K=params.K, covariates=covariates,
            study_ids=expression.index.values, snp_ids=genotype.columns.values, sample_ids=genotype.index.values)
        model.prior_activity = np.ones(params.K) * params.p0k
        forward_fit_procedure(model, verbose=True, update_covariate_weights=True)
        model.save(output.model)




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
