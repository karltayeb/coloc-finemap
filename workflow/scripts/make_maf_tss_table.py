import numpy as np
import pandas as pd
def find(arr, target, low, high):
    if high - low <= 1:
        return high
    
    mid = int((low + high) / 2)
    if arr[mid] > target:
        return find(arr, target, low, mid)
    else:
        return find(arr, target, mid, high)

gencode = pd.read_csv(('/work-zfs/abattle4/lab_data/GTEx_v8/references/'
                      'gencode.v26.GRCh38.genes.gtf'), sep='\t', skiprows=6, header=None)

g = gencode[gencode.iloc[:, 2] =='gene']
tss = g.apply(lambda x: x.values[3] if x.values[6] is '+' else x.values[4], axis=1)
gene_id = g.apply(lambda x: x.values[8].split(';')[0].split('"')[1], axis=1)
gene_info = pd.concat([g.iloc[:, 0], tss, gene_id], keys=['chromosome', 'tss', 'gene'], axis=1)
gene_info_dict = {'chr{}'.format(k): gene_info.loc[gene_info.chromosome == 'chr{}'.format(k)].sort_values('tss') for k in range(1, 23)}

i=0
with open(snakemake.output[0], 'w') as g:
    with open(snakemake.input[0], 'r') as f:
        f.readline()
        for line in f:
            try:
                _, variant_id, ref, alt, alt_freq, obs_ct = line.split('\t')
                chrom, pos = variant_id.split('_')[0], int(variant_id.split('_')[1])
                alt_freq = float(alt_freq)
                maf = np.min([1 - alt_freq, alt_freq])
                low = find(gene_info_dict[chrom].tss.values - pos, -500000, 0,
                           gene_info_dict[chrom].tss.values.size)
                high = find(gene_info_dict[chrom].tss.values - pos, 500000, 0,
                            gene_info_dict[chrom].tss.values.size)
                gene_info_dict[chrom].iloc[low:high].apply(
                    lambda x: print('{}\t{}\t{}\t{:.5f}\t{}'.format(
                        chrom, variant_id, x.gene, maf, pos - x.tss), file=g), axis=1
                )
            except Exception:
                break
            i += 1
            if i % 5000 == 0:
                print(i)


"""
import os
from collections import defaultdict
import pickle
import numpy as np
import pandas as pd
import yaml

maf_bins = np.linspace(0, 1, 51)
tss_bins = np.linspace(-500000, 500000, 41)
bins = {}

for m in maf_bins:
    for t in tss_bins:
        bins[float(m) ,float(t)] = []

for f in os.listdir('output/enrichment/GTEx_maf_tss/'):
    print(f)
    df = pd.read_csv('output/enrichment/GTEx_maf_tss/{}'.format(f), sep='\t', header=None)
    maf_binning = np.digitize(df.iloc[:, 3].astype(float), maf_bins)
    tss_binning = np.digitize(df.iloc[:, 4].astype(float), tss_bins)
    variants = df.iloc[:, 1].values
    for m, t, variant in zip(maf_binning, tss_binning, variants):
        bins[maf_bins[m], tss_bins[t]].append(variant)
        
    with open('output/enrichment/bins', 'wb') as f:
        pickle.dump(bins, f)
"""