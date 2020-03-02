import pandas as pd
import numpy as np

rule create_matched_variant_set:
    input:
        'output/enrichment/{prefix}.bed'
        'output/enrichment/GTEx_maf_tss_binned/bins.{suffix}'
    output:
        'output/enrichment/{prefix}.{suffix}.matched.bed'
    wildcard_constraints:
        suffix='^[x][a-z][a-z]$'
    run:
        # put variants into bins
        bins = []
        results = pd.read_csv(input[0], sep='\t', header=None)
        for chrom, group in results.groupby(0):
            print(chrom)
            df = pd.read_csv('maf/{}.afreq'.format(chrom))
            df.loc[:, 'pos'] = df.ID.apply(lambda x: int(x.split('_')[1]))
            
            for i, record in group.iterrows():
                try:
                    dtss = record[1] - gencode.loc[record[3]].tss

                    pos = int(record[4].split('_')[1])
                    maf = df[df.pos == pos].ALT_FREQS.values[0]
                    maf = np.min([maf, 1-maf])

                    maf_bin = np.digitize(maf, maf_bins)
                    tss_bin = np.digitize(dtss, tss_bins)

                    bins.append((maf_bin, tss_bin))
                except:
                    print('!')

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

suffixes=[
    'xab', 'xac', 'xad', 'xae', 'xaf', 'xag', 'xah', 'xai' ,'xaj', 'xak', 'xal', 'xam',
    'xan', 'xao', 'xap', 'xaq', 'xar', 'xas', 'xat', 'xau', 'xav', 'xaw', 'xax', 'xay', 'xaz',
    'xba', 'xbb', 'xbc', 'xbd', 'xbe', 'xbf', 'xbg', 'xbh', 'xbi' ,'xbj', 'xbk', 'xbl', 'xbm',
    'xbn', 'xbo', 'xbp', 'xbq', 'xbr', 'xbs', 'xbt'
]
rule merge_variant_sets:
    input:
        expand('output/enrichment/{prefix}.{suffix}.matched.bed',
            prefix='{prefix}', suffix=suffixes)
    wildcard_constraints:
        suffix='^[x][a-z][a-z]$',
        prefix= '[^.]+'
    output:
        sorted_bed='output/enrichment/{prefix}.sorted.bed',
        merged_bed='output/enrichment/{prefix}.matched.bed'
    shell:
        'cat {input} | sort -k1,1, -k2,2n > {output.sorted_bed}'
        'merge -i {output.sorted_bed} > {output.merged_bed}'
