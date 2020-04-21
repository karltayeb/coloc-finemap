import pickle
import numpy as np
import pandas as pd

rule format_caviar_data_ld:
    input:
        '{path}/{prefix}.data'
    output:
        ld_matrix = temp('{path}/caviar/{prefix}.data.ld'),
        zscores = temp(expand('{path}/caviar/{prefix}.data.zscores{t}',
            path='{path}', prefix='{prefix}', t=list(range(20))))
    wildcard_constraints:
        prefix='[^/]+'
    run:
        data = pickle.load(open(input[0], 'rb'))
        LD = np.corrcoef(data['X'])
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
        '{path}/caviar/{prefix}.caviar.t{tissue}.log',
        '{path}/caviar/{prefix}.caviar.t{tissue}.post',
        '{path}/caviar/{prefix}.caviar.t{tissue}.set'
    wildcard_constraints:
        prefix='[^/]+'
    shell:
        "workflow/bin/caviar/CAVIAR "
        "-o output/{wildcards.path}/caviar/{wildcards.prefix}.caviar.t{wildcards.tissue} "
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