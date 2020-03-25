import numpy as np
import pandas as pd
import json
"""
given a bed file listing variants
construct a background set matched for MAF and TSS distance
"""
def get_matched_set(input_bed, output_bed):
    maf_bins = np.linspace(0, 1, 51)
    tss_bins = np.linspace(-500000, 500000, 51)

    # put variants into bins
    bins = []
    results = pd.read_csv(input_bed, sep='\t', header=None)
    for chrom, group in results.groupby(0):
        print('\n', chrom, group.shape[0])
        print('\t', end='')
        df = pd.read_csv('maf/{}.afreq'.format(chrom))
        df.loc[:, 'pos'] = df.ID.apply(lambda x: int(x.split('_')[1]))
        
        for i, record in group.iterrows():
            try:
                dtss = record[1] - gencode.loc[record[3]].tss

                pos = record[1]
                maf = df[df.pos == pos].ALT_FREQS.values[0]
                maf = np.min([maf, 1-maf])

                maf_bin = np.digitize(maf, maf_bins)
                tss_bin = np.digitize(dtss, tss_bins)

                bins.append((maf_bin, tss_bin))
            except:
                print('!', end='')

    # select variants from bins
    matched_snps = []
    binning = json.load(open(input[1], 'r'))
    for match in bins:
        try:
            snps = np.array(list(binning[str(match[0])][str(match[1])].keys()))
            snps = snps[np.random.choice(snps.size, 100)]
            matched_snps.append(list(snps))
        except Exception:
            '!'

    with open(output[0], 'w') as f:
        for snp in np.unique(np.concatenate(matched_snps)):
            chromosome, pos = snp.split('_')[:2]
            pos = int(pos)
            print('{}\t{}\t{}\t{}\t{}'.format(chromosome, pos, pos+1, 'matched', snp.strip()), file=f)