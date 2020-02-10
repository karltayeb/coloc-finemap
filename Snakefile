configfile: "config/config.yaml"
rule get_cis_variants:
    output:
        "output/genotypes/{gene}_cis_variants"
    script:
        "workflow/scripts/get_cis_variants.py"

rule random_gene_list:
    output:
        "output/gene_list.txt"
    run:
        import pandas as pd
        import numpy as np
        gencode = pd.read_csv(
            '/work-zfs/abattle4/lab_data/genomic_annotation_data/'
            'gencode.v19.genes.v6p.patched_contigs_TSS.bed',
            sep='\t', header=None
        )
        
        gencode = gencode.loc[~gencode.loc[:, 0].isin(['chrMT', 'chrX', 'chrY'])]
        sites = np.random.choice(gencode.shape[0], 100)
        gencode = gencode.iloc[sites]
        gencode.to_csv(output[0])

rule simulate_single_causal_variant:
    input:
        "output/genotypes/{gene}_cis_variants"
    output:
        "output/simulation/single_causal_variant/{gene}/ld_{linkage}_pve_{pve}_data"
    script:
        "scripts/single_causal_variant.py"

rule fit_cosie_summary:
    input:
        "output/simulation/single_causal_variant/{gene}/ld_{linkage}_pve_{pve}_data"
    output:
        "output/models/cosie_summary/gene_{gene}_ld_{linkage}_pve_{pve}_model"
    script:
        "scripts/fit_cosie_summary.py"

rule fit_cosie_genotype:
    input:
        "output/simulation/single_causal_variant/{gene}/ld_{linkage}_pve_{pve}_data"
    output:
        "output/models/cosie_genotype/gene_{gene}_ld_{linkage}_pve_{pve}_model"
    script:
        "scripts/fit_cosie_genotype.py"
        
