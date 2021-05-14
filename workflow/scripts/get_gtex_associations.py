import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from itertools import islice
import glob
import json

association_indices = {
        x.split('/')[-1].split('.')[0]: json.load(open(x, 'r'))
        for x in glob.glob('/work-zfs/abattle4/karl/cosie_analysis/output/GTEx/index/*.index')
    }

gene_tissue_map = {}
for tissue in association_indices:
    for gene in association_indices[tissue]:
        if gene in gene_tissue_map:
            gene_tissue_map[gene].append(tissue)
        else:
            gene_tissue_map[gene] = [tissue]

def get_associations(test_gene):
    results = {}
    for test_tissue in gene_tissue_map[test_gene]:
        lines = []
        with open(
            '/work-zfs/abattle4/lab_data/GTEx_v8/ciseQTL/'
            'GTEx_Analysis_v8_eQTL_all_associations/{}.allpairs.txt'.format(
                test_tissue), 'r') as f:
            f.seek(association_indices[test_tissue][test_gene])
            last_gene = ''
            for line in f:
                gene = line.split('\t')[0]
                if gene != test_gene and len(lines) > 1:
                    break
                else:
                    lines.append(line)
                    last_gene = gene
        results[test_tissue] = pd.DataFrame([x.strip().split('\t') for x in lines])

    results = pd.concat(results)
    results = results.rename({
        0: 'gene_id',
        1: 'variant_id',
        2: 'tss_distance',
        3: 'ma_samples',
        4: 'ma_count',
        5: 'maf',
        6: 'pval_nominal',
        7: 'slope',
        8: 'slope_se'
    }, axis='columns')

    results.slope = results.slope.astype(float)
    results.slope_se = results.slope_se.astype(float)

    results = results[results.maf.astype(float) > 0.01]
    results.index.levels[0].name = 'tissue'
    results.reset_index(inplace=True)
    results.drop('level_1', axis=1, inplace=True)
    return results

test_gene = snakemake.wildcards.gene
results = get_associations(test_gene)
results.to_csv(snakemake.output.associations)

"""
a = pd.pivot_table(results, index='tissue', columns='variant_id', values='zscore')
a.to_csv(snakemake.output[0], sep='\t')

a = pd.pivot_table(results, index='tissue', columns='variant_id', values='slope')
a.to_csv(snakemake.output[1], sep='\t')

a = pd.pivot_table(results, index='tissue', columns='variant_id', values='slope_se')
a.to_csv(snakemake.output[2], sep='\t')
"""
