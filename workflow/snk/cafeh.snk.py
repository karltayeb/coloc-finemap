rule fit_genotype_model:
    input:
        genotype = 'output/GTEx/{chr}/{gene}/{gene}.raw',
        expression = 'output/GTEx/{chr}/{gene}/{gene}.expression',
        snp2rsid = 'output/GTEx/{chr}/{gene}/{gene}.snp2rsid'
    output:
        model = 'output/GTEx/{chr}/{gene}/{gene}.cafeh_genotype'
    params:
        K = 20,
        p0k = 1.0,
        tolerance = 1e-4
    group: "g"
    run:
        from cafeh.independent_model_ss import CAFEHG
        from cafeh.fitting import forward_fit_procedure

        from utils.misc import load_gtex_genotype, load_gtex_expression
        genotype = load_gtex_genotype(wildcards.gene)
        expression = load_gtex_expression(wildcards.gene)

        snp_ids = genotype.columns.values
        sample_ids = np.intersect1d(genotype.index.values, expression.columns.values)
        study_ids = expression.index.values

        X = np.nan_to_num(genotype.loc[sample_ids].values - np.nanmean(genotype.loc[sample_ids].values, 0)).T
        Y = expression.loc[:, sample_ids].values

        covariates = pd.read_csv(
            '/work-zfs/abattle4/karl/cosie_analysis/output/GTEx/covariates.csv', sep='\t', index_col=[0, 1])
        covariates = covariates.loc[study_ids].loc[:, sample_ids]
        model = CAFEHG(
            X=X, Y=Y, K=params.K, covariates=covariates,
            study_ids=study_ids, snp_ids=snp_ids, sample_ids=sample_ids)
        model.prior_activity = np.ones(params.K) * params.p0k
        model.tolerance = params.tolerance

        # set weight_var reasonably
        model.fit(
            max_iter=1,
            verbose=True,
            update_weights=True,
            update_pi=True,
            ARD_weights=False,
            update_active=False,
            update_covariate_weights=True,
            update_variance=False
        )

        # set tissue_variance to reflect covariates
        model.fit(
            max_iter=1,
            verbose=True,
            update_weights=False,
            update_pi=False,
            ARD_weights=False,
            update_active=False,
            update_covariate_weights=True,
            update_variance=True
        )

        # fit w/o ARD to get good initialization
        model.fit(
            max_iter=20,
            verbose=True,
            update_weights=True,
            update_pi=True,
            ARD_weights=False,
            update_active=False,
            update_covariate_weights=True,
            update_variance=False
        )

        # fit model
        model.fit(
            max_iter=100,
            verbose=True,
            update_weights=True,
            update_pi=True,
            ARD_weights=True,
            update_active=False,
            update_covariate_weights=True,
            update_variance=True
        )
        model.save(output.model)


rule fit_cafeh_genotype_ss:
    input:
        genotype = 'output/GTEx/{chr}/{gene}/{gene}.raw',
        expression = 'output/GTEx/{chr}/{gene}/{gene}.expression',
        snp2rsid = 'output/GTEx/{chr}/{gene}/{gene}.snp2rsid',
        model = 'output/GTEx/{chr}/{gene}/{gene}.cafeh_genotype'
    output:
        model = 'output/GTEx/{chr}/{gene}/{gene}.cafeh_genotype_ss'
    params:
        K = 20,
        p0k = 0.1,
        tolerance = 1e-5
    group: "g"
    run:
        from cafeh.independent_model_ss import CAFEHG
        from cafeh.fitting import forward_fit_procedure

        from utils.misc import load_gtex_genotype, load_gtex_expression
        genotype = load_gtex_genotype(wildcards.gene)
        expression = load_gtex_expression(wildcards.gene)

        snp_ids = genotype.columns.values
        sample_ids = np.intersect1d(genotype.index.values, expression.columns.values)
        study_ids = expression.index.values

        X = np.nan_to_num(genotype.loc[sample_ids].values - np.nanmean(genotype.loc[sample_ids].values, 0)).T
        Y = expression.loc[:, sample_ids].values

        covariates = pd.read_csv(
            '/work-zfs/abattle4/karl/cosie_analysis/output/GTEx/covariates.csv', sep='\t', index_col=[0, 1])
        covariates = covariates.loc[study_ids].loc[:, sample_ids]
        model = CAFEHG(
            X=X, Y=Y, K=params.K, covariates=covariates,
            study_ids=study_ids, snp_ids=snp_ids, sample_ids=sample_ids)

        # initialize model with cafeh_genotype model params
        init_model = pickle.load(open(input.model, 'rb'))
        init_model._decompress_model()
        init_model.__dict__.pop('precompute')
        model.__dict__.update(init_model.__dict__)

        model.prior_activity = np.ones(model.dims['K']) * params.p0k
        model.tolerance = params.tolerance

        # fit with spike and slab
        model.fit(
            max_iter=100,
            verbose=True,
            update_weights=True,
            update_pi=True,
            ARD_weights=True,
            update_active=True,
            update_covariate_weights=True,
            update_variance=True
        )
        model.save(output.model)


rule generate_snp_report:
    input:
        model = 'output/GTEx/{chr}/{gene}/{gene}.cafeh_genotype_ss',
        snp2rsid = 'output/GTEx/{chr}/{gene}/{gene}.snp2rsid'
    output:
        report = 'output/GTEx/{chr}/{gene}/{gene}.cafeh_genotype_ss.variant_report'
    script:
        '../../workflow/scripts/report_genotype.py'

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


rule fit_gwas_gtex_z_impute_cafeh:
    input:
        genotype_gtex = 'output/GTEx/{chr}/{gene}/{gene}.raw',
        associations = 'output/GTEx/{chr}/{gene}/{gene}.associations',
        v2r = 'output/GTEx/{chr}/{gene}/{gene}.snp2rsid'
    output:
        z_imp_model='output/{study}/{phenotype}/{chr}/{gene}/{gene}.{phenotype}.z_imputed.css'
    params:
        impute=True,
        K=20
    script:
        '../../workflow/scripts/ukbb_gtex_cafeh_ss.py'

rule fit_gwas_gtex_z_cafeh:
    input:
        genotype_gtex = 'output/GTEx/{chr}/{gene}/{gene}.raw',
        associations = 'output/GTEx/{chr}/{gene}/{gene}.associations',
        v2r = 'output/GTEx/{chr}/{gene}/{gene}.snp2rsid'
    output:
        z_imp_model='output/{study}/{phenotype}/{chr}/{gene}/{gene}.{phenotype}.z.css'
    params:
        impute=False,
        K=20
    script:
        '../../workflow/scripts/ukbb_gtex_cafeh_ss.py'
