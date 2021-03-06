rule format_coloc_data:
    input:
        '{path}/{prefix}.data'
    output:
        coloc_data = temp('{path}/{prefix}.coloc_data'),
    wildcard_constraints:
        prefix='[^/]+'
    run:
        data = pickle.load(open(input[0], 'rb'))
        print('computing summary stats')
        X = data['X']
        Y = data['Y']
        n = Y.shape[1]
        xx = np.einsum('nm,nm->n', X, X)
        B = Y@X.T / xx
        S2 = np.sum((Y[:, None] - B[..., None] * X)**2, 2) / (xx * (n-2))
        data = {
        	'Y': Y,
            'sdY': Y.std(1),
        	'betas': B,
        	's2': S2
        }
        pickle.dump(data, open(output.coloc_data, 'wb'))

rule run_coloc:
    input:
        "output/simulation/{simulation}/{settings}/gene_{gene}/data"
    output:
        "output/simulation/{simulation}/{settings}/gene_{gene}/coloc"
    wildcard_constraints:
        simulation = "(?!\/)[^\/]+(?=\/)"
    script:
        "../../workflow/scripts/run_coloc.R"

paths = glob.glob('output/sim/multiple/chr2/*')
rule run_coloc_sim:
	input:
		expand('{path}/genotype.sim.t20.pve{pve}.coloc',
			path=paths, pve=['01', '05', '10'])

rule run_coloc2:
    input:
        "{path}/{prefix}.coloc_data"
    output:
        "{path}/{prefix}.coloc"
    wildcard_constraints:
        prefix='[^/]+'
    script:
        "../../workflow/scripts/run_coloc.R"

