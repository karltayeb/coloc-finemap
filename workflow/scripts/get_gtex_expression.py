import json
import glob
import pandas as pd
import numpy as np

test_gene = snakemake.wildcards.gene
association_indices = {
    x.split('/')[-1].split('.')[0]: json.load(open(x, 'r'))
    for x in snakemake.input
}
gene_tissue_map = json.load(open('data/GTEx/gene2tissues'))
gene_expression_df = {}
for tissue in gene_tissue_map[test_gene]:
    expression = pd.read_csv('/work-zfs/abattle4/lab_data/GTEx_v8/ciseQTL/GTEx_Analysis_v8_eQTL_expression_matrices/'
                             + '/{}.v8.normalized_expression.bed.gz'.format(tissue), sep='\t')
    gene_expression_df[tissue] = expression.loc[expression.gene_id == test_gene]

"""
dfs = {}
paths = glob.glob('/work-zfs/abattle4/lab_data/GTEx_v8/ciseQTL/GTEx_Analysis_v8_eQTL_expression_matrices/*.bed.gz')

for name, group in all_expression[all_expression.gene_id.isin(protein_coding.gene)].groupby(['#chr', 'gene_id']):
    chrom, gene = name
    print(gene)
    group.index = group.index.droplevel(1)
    group = group.iloc[:, 1:-3]
    group = group.loc[:, ~np.all(np.isnan(group), axis=0)]
    os.makedirs('output/GTEx/{}/{}'.format(chrom, gene))
    group.to_csv('output/GTEx/{}/{}/{}.expression'.format(chrom, gene, gene), sep='\t')
"""
gene_expression = pd.concat(gene_expression_df, sort=False).droplevel(1)
gene_expression = gene_expression.iloc[:, 4:]
gene_expression = gene_expression.loc[:, ~np.all(np.isnan(gene_expression), 0)]
gene_expression.to_csv(snakemake.output[0], sep='\t')