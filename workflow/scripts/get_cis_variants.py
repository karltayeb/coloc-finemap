import pandas as pd
import numpy as np

def get_cis_variants(genotype, tss, window=500000, size=1000):
    pos = np.array([int(x.split('_')[1]) for x in genotype.index.values])
    #cis_variants = genotype.iloc[np.abs(pos - tss) < window]
    cis_variants = genotype.iloc[
        np.arange(genotype.shape[0])[np.argsort(np.abs(pos - tss))[:size]]
    ]
    return cis_variants

gencode = pd.read_csv(
    '/work-zfs/abattle4/lab_data/genomic_annotation_data/gencode.v19.genes.v6p.patched_contigs_TSS.bed', sep='\t')
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
cis_variants.to_csv(snakemake.output[0])