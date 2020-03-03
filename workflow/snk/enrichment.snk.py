import pandas as pd
import numpy as np
from coloc.cafeh import CAFEH

suffixes=[
    'xab', 'xac', 'xad', 'xae', 'xaf', 'xag', 'xah', 'xai' ,'xaj', 'xak', 'xal', 'xam',
    'xan', 'xao', 'xap', 'xaq', 'xar', 'xas', 'xat', 'xau', 'xav', 'xaw', 'xax', 'xay', 'xaz',
    'xba', 'xbb', 'xbc', 'xbd', 'xbe', 'xbf', 'xbg', 'xbh', 'xbi' ,'xbj', 'xbk', 'xbl', 'xbm',
    'xbn', 'xbo', 'xbp', 'xbq', 'xbr', 'xbs', 'xbt'
]

rule make_bed_summary_genotype:
    input:
        data='output/{path}/data',
        model='output/{path}/model_genotype'
    output:
        'ouput/{path}/genotype.credible_sets.bed'
    script:
        'workflow/scripts/make_credible_set_bed_genotype.py'

rule create_matched_variant_set:
    input:
        'output/{path}/{prefix}.bed',
        'output/enrichment/GTEx_maf_tss_binned/bins.{suffix}'
    output:
        temp('output/enrichment/{prefix}.{suffix}.matched.bed')
    wildcard_constraints:
        prefix= '[^.]+'
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
                '!'

        with open(output[0], 'w') as f:
            for snp in np.unique(np.concatenate(matched_snps)):
                chromosome, pos = snp.split('_')[:2]
                pos = int(pos)
                print('{}\t{}\t{}\t{}\t{}'.format(chromosome, pos, pos+1, 'matched', snp.strip()), file=f)

rule merge_variant_sets:
    input:
        expand('output/{path}/{prefix}.{suffix}.matched.bed',
            prefix='{prefix}', suffix=suffixes)
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
