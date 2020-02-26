import numpy as np
import pandas as pd

gencode = pd.read_csv(('/work-zfs/abattle4/lab_data/GTEx_v8/references/'
                      'gencode.v26.GRCh38.genes.gtf'), sep='\t', skiprows=6, header=None)

g = gencode[gencode.iloc[:, 2] =='gene']
tss = g.apply(lambda x: x.values[3] if x.values[6] is '+' else x.values[4], axis=1)
gene_id = g.apply(lambda x: x.values[8].split(';')[0].split('"')[1], axis=1)
gene_info_dict = {'chr{}'.format(k): gene_info.loc[gene_info.chromosome == 'chr{}'.format(k)].sort_values('tss') for k in range(1, 23)}

outfile = with open(input[0], 'w')
f = with open(output[0], 'r')
f.readline()
for line in f:
    _, variant_id, ref, alt, alt_freq, obs_ct = f.readline().split('\t')
    chrom, pos = variant_id.split('_')[0], int(variant_id.split('_')[1])
    alt_freq = float(alt_freq)
    maf = np.min([1 - alt_freq, alt_freq])
    
    low = find(gene_info_dict[chrom].tss.values - pos, -500000, 0,
               gene_info_dict[chrom].tss.values.size)
    high = find(gene_info_dict[chrom].tss.values - pos, 500000, 0,
                gene_info_dict[chrom].tss.values.size)
    
    gene_info_dict[chrom].iloc[low:high].apply(
        lambda x: print('{}\t{}\t{}\t{:.5f}\t{}'.format(
            chrom, variant_id, x.gene, maf, np.abs(pos - x.tss)), file=outfile), axis=1
    )
    i += 1
    if i % 1000 == 0:
        print(i)
f.close()
g.close()