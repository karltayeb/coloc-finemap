import json
import glob
import pandas as pd

test_gene = snakemake.wildcards.gene
association_indices = {
    x.split('/')[-1].split('.')[0]: json.load(open(x, 'r'))
    for x in snakemake.input
}
gene_tissue_map = {}
for tissue in association_indices:
    for gene in association_indices[tissue]:
        if gene in gene_tissue_map:
            gene_tissue_map[gene].append(tissue)
        else:
            gene_tissue_map[gene] = [tissue]

gene_expression_df = {}
for tissue in gene_tissue_map[test_gene]:
    expression = pd.read_csv('/work-zfs/abattle4/lab_data/GTEx_v8/ciseQTL/GTEx_Analysis_v8_eQTL_expression_matrices/'
                             + '/{}.v8.normalized_expression.bed.gz'.format(tissue), sep='\t')
    gene_expression_df[tissue] = expression.loc[expression.gene_id == test_gene]

gene_expression = pd.concat(gene_expression_df, sort=False).droplevel(1)
gene_expression = gene_expression.iloc[:, 4:]

gene_expression.to_csv(snakemake.output[0], sep='\t')