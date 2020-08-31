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
    run:
        import pickle
        import numpy as np
        import pandas as pd
        from cafeh.cafeh_ss import CAFEH
        from cafeh.misc import plot_components
        from utils.misc import *

        def _get_minalpha(pi):
            """
            report the minimum alpha value to include this snp in cs
            """
            argsort = np.flip(np.argsort(pi))
            resort = np.argsort(argsort)
            cumsum = np.cumsum(pi[argsort])
            minalpha = np.roll(cumsum, 1)
            minalpha[0] = 0
            return minalpha[resort]

        def get_minalpha(model):
            return  pd.DataFrame(
                np.array([_get_minalpha(model.pi[k]) for k in range(model.dims['K'])]).T,
                index=model.snp_ids
            )

        gene = snakemake.wildcards.gene

        # load a model
        model = pickle.load(open(input.model, 'rb'))
        model._decompress_model()

        study_pip = model.get_study_pip().T
        table = study_pip.reset_index().melt(id_vars='index').rename(columns={
            'index': 'variant_id',
            'variable': 'study',
            'value': 'pip' 
        })

        all_pip = pd.DataFrame({'pip': model.get_pip(), 'variant_id': model.snp_ids, 'study': 'all'})
        table = pd.concat([table, all_pip], sort=True)

        v2r = load_var2rsid(gene)
        table.loc[:, 'rsid'] = table.variant_id.apply(lambda x: v2r.get(x, '-'))

        top_component = pd.Series(model.pi.argmax(0), index=model.snp_ids).to_dict()
        table.loc[:, 'top_component'] = table.variant_id.apply(lambda x: top_component.get(x))

        minalpha = get_minalpha(model).to_dict()
        table.loc[:, 'alpha'] = [minalpha.get(k).get(v) for k, v in zip(table.top_component.values, table.variant_id.values)]

        rank = pd.DataFrame({k: np.argsort(np.flip(np.argsort(model.pi[k]))) for k in range(model.dims['K'])}, index=model.snp_ids).to_dict()
        table.loc[:, 'rank'] = [rank.get(k).get(v) for k, v in zip(table.top_component.values, table.variant_id.values)]

        active = pd.DataFrame(model.active, index=model.study_ids)
        active.loc['all'] = (model.active.max(0) > 0.5).astype(int)
        active = active.to_dict()
        table.loc[:, 'p_active'] = [active.get(k).get(s) for k, s in zip(table.top_component.values, table.study.values)]

        pi = pd.Series(model.pi.max(0), index=model.snp_ids).to_dict()
        table.loc[:, 'pi'] = table.variant_id.apply(lambda x: pi.get(x))

        table.loc[:, 'chr'] = get_chr(gene)
        table.loc[:, 'start'] = table.variant_id.apply(lambda x: int(x.split('_')[1]))
        table.loc[:, 'end'] = table.start + 1
        table.loc[:, 'gene'] = gene

        table.loc[:, ['chr', 'start', 'end', 'variant_id', 'rsid', 'study', 'gene', 'pip', 'pi', 'top_component', 'p_active', 'alpha']]

        table = table.loc[:, ['chr', 'start', 'end', 'variant_id', 'rsid', 'study', 'pip', 'top_component', 'p_active', 'pi', 'alpha', 'rank']]
        small_table = table[table.p_active > 0.5].sort_values(by=['chr', 'start'])
        small_table.to_csv(output.report, sep='\t', index=None

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
