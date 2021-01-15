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


rule fit_susie:
    input:
        genotype = 'output/GTEx/{chr}/{gene}/{gene}.raw',
        expression = 'output/GTEx/{chr}/{gene}/{gene}.expression',
        snp2rsid = 'output/GTEx/{chr}/{gene}/{gene}.snp2rsid'
    output:
        'output/GTEx/{chr}/{gene}/{gene}.susie.variant_report'
    params:
        K = 20,
        p0k = 1.0,
        tolerance = 1e-4
    group: "g"
    run:
        from cafeh.independent_model_ss import CAFEHG
        from cafeh.fitting import fit_all
        from cafeh.model_queries import summary_table
        from utils.misc import load_gtex_genotype, load_gtex_expression

        def make_table(model, gene):
            table = summary_table(model)

            # annotate table
            if table.shape[0] > 0:
                v2r = load_var2rsid(gene)
                table.loc[:, 'rsid'] = table.variant_id.apply(lambda x: v2r.get(x, '-'))
                table.loc[:, 'chr'] = table.variant_id.apply(lambda x: (x.split('_')[0]))
                table.loc[:, 'start'] = table.variant_id.apply(lambda x: int(x.split('_')[1]))
                table.loc[:, 'end'] = table.start + 1
                table.loc[:, 'gene'] = gene
            
            table = table.loc[:, ['chr', 'start', 'end', 'variant_id', 'rsid', 'study', 'pip', 'top_component', 'p_active', 'pi', 'alpha', 'rank', 'effect', 'effect_var']]
            return table

        genotype = load_gtex_genotype(wildcards.gene)
        expression = load_gtex_expression(wildcards.gene)

        snp_ids = genotype.columns.values
        sample_ids = np.intersect1d(genotype.index.values, expression.columns.values)
        study_ids = expression.index.values

        genotype = genotype.loc[sample_ids]
        Y = expression.loc[:, sample_ids].values

        covariates = pd.read_csv(
            '/work-zfs/abattle4/karl/cosie_analysis/output/GTEx/covariates.csv',
            sep='\t', index_col=[0, 1]
        )

        tables = []
        for i in range(Y.shape[0]):
            print(study_ids[i])
            y = Y[i]
            sample_mask = ~np.isnan(y)
            y = y[sample_mask]

            X = genotype[sample_mask].values
            X = np.nan_to_num(X) - np.nanmean(X, 0)
            X = X / X.std(0)

            snp_mask = ~np.any(np.isnan(X), 0)
            X = X[:, snp_mask]

            cov = covariates[sample_ids[sample_mask]][
                covariates.index.get_level_values(0) == study_ids[i]]
            cov = cov.T.values
            H = cov @ np.linalg.pinv(cov)

            y = y - H @ y
            X = X - H @ X

            model = CAFEHG(
                X=X.T, Y=y[None], K=5,
                study_ids=study_ids[[[i]]],
                snp_ids=snp_ids[snp_mask], sample_ids=sample_ids[sample_mask])
            model.prior_activity = np.ones(5) * 0.1
            fit_all(model, update_active=False, ARD_weights=False, max_iter=100)
            fit_all(model, update_active=True, ARD_weights=True, max_iter=100)
            tables.append(make_table(model, wildcards.gene))

        table = pd.concat(tables)
        table = table.sort_values(['chr', 'start'])
        table.to_csv(output[0], sep='\t')

rule generate_snp_report:
    input:
        model = 'output/GTEx/{chr}/{gene}/{gene}.cafeh_genotype_ss',
        snp2rsid = 'output/GTEx/{chr}/{gene}/{gene}.snp2rsid'
    output:
        report = 'output/GTEx/{chr}/{gene}/{gene}.cafeh_genotype_ss.variant_report'
    run:
        from cafeh.independent_model_ss import CAFEHG
        from cafeh.fitting import fit_all
        from cafeh.model_queries import summary_table
        from utils.misc import load_gtex_genotype, load_gtex_expression
        print('!!!!!!!!!!!')
        def make_table(model, gene):
            table = summary_table(model)

            # annotate table
            if table.shape[0] > 0:
                v2r = load_var2rsid(gene)
                table.loc[:, 'rsid'] = table.variant_id.apply(lambda x: v2r.get(x, '-'))
                table.loc[:, 'chr'] = table.variant_id.apply(lambda x: (x.split('_')[0]))
                table.loc[:, 'start'] = table.variant_id.apply(lambda x: int(x.split('_')[1]))
                table.loc[:, 'end'] = table.start + 1
                table.loc[:, 'gene'] = gene

            table = table.loc[:, ['chr', 'start', 'end', 'variant_id', 'rsid', 'study', 'pip', 'top_component', 'p_active', 'pi', 'alpha', 'rank', 'effect', 'effect_var']]
            return table
        # load a model
        gene = wildcards.gene
        model = pickle.load(open(input.model, 'rb'))
        model._decompress_model()
        table = make_table(model, gene)
        table.to_csv(output.report, sep='\t', index=None)

rule fit_gwas_gtex_z_cafeh:
    input:
        genotype_gtex = 'output/GTEx/{chr}/{gene}/{gene}.raw',
        associations = 'output/GTEx/{chr}/{gene}/{gene}.associations',
        sumstats='output/{study}/{phenotype}/{phenotype}.tsv.bgz',
        tabix_index='output/{study}/{phenotype}/{phenotype}.tsv.bgz.tbi',
        v2r = 'output/GTEx/{chr}/{gene}/{gene}.snp2rsid'
    output:
        variant_report='output/{study}/{phenotype}/{chr}/{gene}/{gene}.{phenotype}.z.variant_report',
        coloc_report='output/{study}/{phenotype}/{chr}/{gene}/{gene}.{phenotype}.z.coloc_report',
        model='output/{study}/{phenotype}/{chr}/{gene}/{gene}.{phenotype}.z.css'
    params:
        K=20,
        zscore=True
    group: 'report'
    script:
        '../../workflow/scripts/ukbb_gtex_cafeh_ss.py'

rule fit_gwas_gtex_z_cafeh:
    input:
        genotype_gtex = 'output/GTEx/{chr}/{gene}/{gene}.raw',
        associations = 'output/GTEx/{chr}/{gene}/{gene}.associations',
        sumstats='output/{study}/{phenotype}/{phenotype}.tsv.bgz',
        tabix_index='output/{study}/{phenotype}/{phenotype}.tsv.bgz.tbi',
        v2r = 'output/GTEx/{chr}/{gene}/{gene}.snp2rsid'
    output:
        variant_report='output/{study}/{phenotype}/{chr}/{gene}/{gene}.{phenotype}.variant_report',
        coloc_report='output/{study}/{phenotype}/{chr}/{gene}/{gene}.{phenotype}.coloc_report',
        model='output/{study}/{phenotype}/{chr}/{gene}/{gene}.{phenotype}.css'
    params:
        K=20,
        zscore=False
    group: 'report'
    script:
        '../../workflow/scripts/ukbb_gtex_cafeh_ss.py'

rule fit_gwas_gtex_z_pairwise_cafeh:
    input:
        genotype_gtex = 'output/GTEx/{chr}/{gene}/{gene}.raw',
        associations = 'output/GTEx/{chr}/{gene}/{gene}.associations',
        v2r = 'output/GTEx/{chr}/{gene}/{gene}.snp2rsid'
    output:
        variant_report='output/{study}/{phenotype}/{chr}/{gene}/{gene}.{phenotype}.z.pairwise.variant_report',
        associations='output/{study}/{phenotype}/{chr}/{gene}/{gene}.{phenotype}.z.pairwise.coloc_report',
        model='output/{study}/{phenotype}/{chr}/{gene}/{gene}.{phenotype}.z.pairwise.css'
    params:
        impute=False,
        K=5
    group: 'report'
    script:
        '../../workflow/scripts/cafeh_single_tissue_x_gwas.py'

