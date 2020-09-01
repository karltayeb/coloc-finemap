import pandas as pd
import numpy as np
from coloc.cafeh import CAFEH

suffixes=[
    'xab', 'xac', 'xad', 'xae', 'xaf', 'xag', 'xah', 'xai' ,'xaj', 'xak', 'xal', 'xam',
    'xan', 'xao', 'xap', 'xaq', 'xar', 'xas', 'xat', 'xau', 'xav', 'xaw', 'xax', 'xay', 'xaz',
    'xba', 'xbb', 'xbc', 'xbd', 'xbe', 'xbf', 'xbg', 'xbh', 'xbi' ,'xbj', 'xbk', 'xbl', 'xbm',
    'xbn', 'xbo', 'xbp', 'xbq', 'xbr', 'xbs', 'xbt'
]

rule create_matched_variant_set:
    """
    """
    input:
        'output/{path}/{prefix}.bed',
        'output/enrichment/GTEx_maf_tss_binned/bins.{suffix}'
    output:
        temp('output/enrichment/{prefix}.{suffix}.matched.bed')
    wildcard_constraints:
        prefix= '(\w.)+'
    run:
        maf_bins = np.linspace(0, 1, 51)
        tss_bins = np.linspace(-500000, 500000, 51)

        # put variants into bins
        bins = []
        results = pd.read_csv(input[0], sep='\t', header=None)
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
                print('!', end='')

        with open(output[0], 'w') as f:
            for snp in np.unique(np.concatenate(matched_snps)):
                chromosome, pos = snp.split('_')[:2]
                pos = int(pos)
                print('{}\t{}\t{}\t{}\t{}'.format(chromosome, pos, pos+1, 'matched', snp.strip()), file=f)

rule merge_variant_sets:
    input:
        expand('output/enrichment/{prefix}.{suffix}.matched.bed',
            path='{path}', prefix='{prefix}', suffix=suffixes)
    wildcard_constraints:
        prefix= '[^.]+'
    output:
        matched='output/enrichment/{prefix}.matched.bed',
        sorted='output/enrichment/{prefix}.matched.sorted.bed',
        merged='output/enrichment/{prefix}.matched.merged.bed'
    shell:
        'cat {input} > {output.matched}'
        '\nsort -k1,1 -k2,2n {output.matched} > {output.sorted}'
        '\nbedtools merge -i {output.sorted} > {output.merged}'


rule gtex_get_variant_sets:
    input:
        'output/GTEx/variant_reports/{tissue}.all_genes.variant_report'
    output:
        test_set = 'output/GTEx/enrichment/{tissue}.test.bed',
        bacgkround_set = 'output/GTEx/enrichment/{tissue}.background.bed'
    run:
        import numpy as np
        import pandas as pd
        from tqdm import tqdm

        def pair2bin(dtss, maf):
            """
            put (dtss, maf) pair into bin
            """
            dtss_bin = int((dtss + 1e6) / 20000)
            maf_bin = '{:.2f}'.format(maf)
            return '{}/{}.{}'.format(dtss_bin, dtss_bin, maf_bin)

        df = pd.read_csv('../../output/GTEx/variant_reports/Liver.all_genes.variant_report', sep='\t')

        df = df[
            (df.p_active > 0.9)
            & (df.pip > 0.2)
            & (df.alpha < 0.95)
        ]

        # put variant, gene pair into bins
        df.loc[:, 'dtss'] = df.start - df.tss
        df.loc[:, 'bin'] = [pair2bin(dtss, maf) for dtss, maf in zip(df.dtss, df.maf)]

        # count number of variants in each bin
        bins = df.bin.value_counts().to_dict()

        background = []
        for b, count in tqdm(list(bins.items())):
            bin_path = '../../output/GTEx/maf_dtss_binned/{}.bed'.format(b)
            bin_df = pd.read_csv(bin_path, sep='\t', header=None)
            background.append(bin_df.iloc[np.random.choice(bin_df.shape[0], count*5, replace=True)])

        background_df = pd.concat(background)

        outfile = 'test_set'
        df.sort_values(by=['chr', 'start']).drop_duplicates(['chr', 'start'])\
            .to_csv(outfile, sep='\t', header=False, index=False)

        outfile = 'test_background'
        background_df.sort_values(by=[0, 1]).drop_duplicates([0, 1])\
            .to_csv(outfile, sep='\t', header=False, index=False)
