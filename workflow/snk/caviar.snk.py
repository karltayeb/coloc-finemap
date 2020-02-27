rule format_caviar_data_ld:
    input:
        'output/{path}/data'
    output:
        ld_matrix = temp('output/{path}/caviar/data.ld{tissue}')
    run:
        LD = pickle.load(open(input[0], 'rb'))['LD']
        if np.ndim(LD) == 3:
            LD = LD[int(wildcards.tissue)]
        np.savetxt(fname=output.ld_matrix, X=LD, delimiter='\t')

rule format_caviar_data_zscore:
    input:
        'output/{path}/data'
    output:
        z_scores = temp('output/{path}/caviar/data.z{tissue}')
    run:
        data = pickle.load(open(input[0], 'rb'))
        zscores = pd.DataFrame(data['zscores'][int(wildcards.tissue)])
        zscores.to_csv(output.z_scores, sep='\t', header=None)

rule run_caviar:
    input:
        ld_matrix = 'output/{path}/caviar/data.ld{tissue}',
        z_scores = 'output/{path}/caviar/data.z{tissue}'
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