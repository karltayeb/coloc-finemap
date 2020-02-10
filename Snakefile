configfile: "config.yaml"
rule get_cis_variants:
    input:
    output:
        "output/genotypes/{gene}_cis_variants"
    run:
        import pandas as pd
        import numpy as np

        def get_cis_variants(genotype, tss, window=500000, size=1000):
            pos = np.array([int(x.split('_')[1]) for x in genotype.index.values])
            #cis_variants = genotype.iloc[np.abs(pos - tss) < window]
            cis_variants = genotype.iloc[
                np.arange(genotype.shape[0])[np.argsort(np.abs(pos - tss))[:size]]
            ]
            return cis_variants
        gencode = pd.read_csv('/work-zfs/abattle4/lab_data/genomic_annotation_data/gencode.v19.genes.v6p.patched_contigs_TSS.bed', sep='\t')
        gene = gencode.loc[gencode.iloc[:, 3] == wildcards.gene]
        chr_num = gene.iloc[0, 0]
        tss = gene.iloc[0, 1]
        gene_name = gene.iloc[0, 3]

        genotype_path = ("/work-zfs/abattle4/lab_data/"
            "GTEx_v8_trans_eqtl_data_processed_by_brian/processed_genotypes/"
            "GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_{}_dosage_MAF_05.txt".format(chr_num))
        genotype = pd.read_csv(genotype_path, sep='\t', index_col=0, na_values='-', low_memory=True)
        genotype = genotype.iloc[~np.any(np.isnan(genotype.values), 1)]
        
        cis_variants = get_cis_variants(genotype, tss)
        cis_variants.to_csv(output[0])

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
        
