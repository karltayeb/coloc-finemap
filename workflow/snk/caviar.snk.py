rule format_caviar_data_ld:
    input:
        '{path}/{prefix}.data'
    output:
        ld_matrix = temp('{path}/caviar/{prefix}.data.ld')
        zscores = temp(enumerate('output/{path}/caviar/{prefix}.data.zscores{tissue}', t=range(20)))

    run:
        data = pickle.load(open(input[0], 'rb'))
        LD = np.corrcoef(data['X'])
        LD = pickle.load(open(input[0], 'rb'))['LD']
        np.savetxt(fname=output.ld_matrix, X=LD, delimiter='\t')

        X = data['X'],
        Y = data['Y']
        n = Y.shape[1]
        xx = np.einsum('nm,nm->n', X, X)
        B = Y@X.T / xx
        S2 = np.sum((Y[:, None] - B[..., None] * X)**2, 2) / (xx * (n-2))
        Z = B / np.sqrt(S2)
        for t, z in enumerate(Z):
            zscores = pd.DataFrame(z, index=data['snp_ids'])
            zscores.to_csv(output.z_scores[t], sep='\t', header=None)

rule run_caviar:
    input:
        ld_matrix = '{path}/caviar/{prefix}.data.ld',
        z_scores = '{path}/caviar/{prefix}.data.zscores{tissue}'
    output:
        '{path}/caviar/caviar.t{tissue}.log',
        '{path}/caviar/caviar.t{tissue}.post',
        '{path}/caviar/caviar.t{tissue}.set'
    shell:
        "workflow/bin/caviar/CAVIAR "
        "-o output/{wildcards.path}/caviar/caviar_t{wildcards.tissue} "
        "-l {input.ld_matrix} "
        "-z {input.z_scores} "
        "-c 3"

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