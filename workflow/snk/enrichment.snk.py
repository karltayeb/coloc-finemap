rule get_tissue_expressed_genes_bed:
    output:
        'output/GTEx/tissue_expressed_genes/{tissue}.genes.bed'
    run:
        import pandas as pd
        import numpy as np

        tissue = wildcards.tissue
        expression_path = '/work-zfs/abattle4/lab_data/GTEx_v8/ciseQTL/'\
            'GTEx_Analysis_v8_eQTL_expression_matrices/{}.v8.normalized_expression.bed.gz'.format(tissue)

        gencode_path = '/work-zfs/abattle4/lab_data/annotation/gencode.v26/gencode.v26.annotation.gene.txt'

        print('load data')
        gencode = pd.read_csv(gencode_path, sep='\t')
        tissue_gene = pd.read_csv(expression_path, sep='\t', usecols=np.arange(4))

        print('filtering')
        pseudo_gene_types = [x for x in gencode.gene_type.unique() if 'pseudo' in x]
        gencode_tissue = gencode[gencode.gene_id.isin(tissue_gene.gene_id) & ~(gencode.gene_type.isin(pseudo_gene_types))]

        lines = []
        for _, row in gencode_tissue.iterrows():
            if row.strand == '-':
                start = row.end_pos
            else:
                start = row.start_pos
            line = {
                'chr': row.chr, 'start': start, 'end': start + 1,
                'gene_id': row.gene_id, 'strand': 1 if (row.strand == '+') else -1
            }
            lines.append(line)
            
        tissue_gene_bed = pd.DataFrame(lines)

        print('save')
        tissue_gene_bed.sort_values(['chr', 'start']).to_csv(output[0], sep='\t', index=None, header=False)

rule get_tissue_specific_variant_gene_pairs:
    input:
        'output/GTEx/tissue_expressed_genes/{tissue}.genes.bed'
    output:
        'output/GTEx/tissue_specific_variant_gene_pairs/{tissue}.variant_gene_pairs.bed'
    run:
        import subprocess
        cmd = "bedtools closest -a output/GTEx/GTEx.afreq.ldscore.bed -b {input} -d "\
            "| awk \'{{print $1, $2, $3, $4, $5, $10, $12, $13, $14, $15,$15* $16}}\' > {output}"
        cmd = cmd.format(input=input[0], output=output[0])
        print(cmd)
        subprocess.run(cmd, shell=True)


rule bin_tissue_specific_variant_gene_pairs:
    input:
        'output/GTEx/tissue_specific_variant_gene_pairs/{tissue}.variant_gene_pairs.bed'
    output:
        expand('output/GTEx/tissue_specific_variant_gene_pairs/{tissue}/{tissue}.variant_gene_pairs.maf_bin_{maf_bin}.bed',
            maf_bin=np.arange(1, 26), allow_missing=True)
    run:
        import pandas as pd
        import numpy as np

        def digitize_column(df, col, n_bins):
            """
            add {col}_bin and {col}_bin_range to df
            discritize numeric valued col in n_bins quantile bins
            """
            bins = np.quantile(df.loc[:, col], np.linspace(0, 1, n_bins+1))
            bins[-1] = bins[-1] + 1 # so that we dont have a singleton bin
            print(bins)
            
            bin_col = '{}_bin'.format(col)
            bin_range_col = '{}_bin_range'.format(col)
            bin2range = {x: '[{}, {})'.format(bins[x-1], bins[x]) for x in range(1, n_bins + 1)}
            
            df.loc[:, bin_col] = np.digitize(df.loc[:, col], bins)
            df.loc[:, bin_range_col] = df.loc[:, bin_col].apply(lambda x: bin2range.get(x))

        print('loading')
        df = pd.read_csv(input[0], sep=' ', usecols=[0, 1, 2, 3, 4, 5, 10], header=None)
        df.columns = np.array(['chr', 'start', 'end', 'variant_id', 'maf', 'ldscore', 'tss'])
        df = df[(df.tss > -1e6) & (df.tss < 1e6)]

        print('binning')
        digitize_column(df, 'maf', 25)
        digitize_column(df, 'ldscore', 10)
        digitize_column(df, 'tss', 40)

        for maf_bin, group in df.groupby('maf_bin'):
            save_path = 'output/GTEx/tissue_specific_variant_gene_pairs/{tissue}/{tissue}.variant_gene_pairs.maf_bin_{maf_bin}.bed'.format(
                tissue=wildcards.tissue, maf_bin=maf_bin)
            print(save_path)
            group.to_csv(save_path, sep='\t', index=False)

rule gtex_get_variant_sets:
    input:
        'output/GTEx/variant_reports/{tissue}.all_genes.variant_report'
    output:
        test = 'output/GTEx/enrichment/{tissue}.test.bed',
        background = 'output/GTEx/enrichment/{tissue}.background.bed'
    run:
        tissue = wildcards.tissue

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

        df = pd.read_csv('output/GTEx/variant_reports/{}.all_genes.variant_report'.format(tissue), sep='\t')

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
            try:
                bin_path = 'output/GTEx/maf_dtss_binned/{}.bed'.format(b)
                bin_df = pd.read_csv(bin_path, sep='\t', header=None)
                background.append(bin_df.iloc[np.random.choice(bin_df.shape[0], count*5, replace=True)])
            except Exception:
                continue
        background_df = pd.concat(background)

        df.loc[:, 'chr_num'] = df.loc[:, 'chr'].apply(lambda x: int(x.replace('chr', '')))
        df.sort_values(by=['chr_num', 'start']).drop_duplicates(['chr', 'start'])\
            .to_csv(output.test, sep='\t', header=False, index=False)

        background_df.loc[:, 'chr_num'] = background_df.loc[:, 0].apply(lambda x: int(x.replace('chr', '')))
        background_df.sort_values(by=['chr_num', 1]).drop_duplicates([0, 1])\
            .to_csv(output.background, sep='\t', header=False, index=False)


rule gtex_filtered_variants_and_background:
    input:
        'output/GTEx/variant_reports/{tissue}.all_genes.variant_report'
    output:
        test = 'output/GTEx/enrichment/{analysis_id}/{tissue}.test.bed',
        background = 'output/GTEx/enrichment/{analysis_id}/{tissue}.background.bed'
    params:
        filters = lambda wildcards: config['enrichment_filters'][wildcards.analysis_id]
    run:
        tissue = wildcards.tissue
        analysis_id = wildcards.analysis_id
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

        if analysis_id == 'eqtl':
            genes = pd.read_csv('output/GTEx/protein_coding_autosomal_egenes.txt', sep='\t').gene.values
            eqtls = pd.read_csv(
                '/work-zfs/abattle4/lab_data/GTEx_v8/ciseQTL/GTEx_Analysis_v8_eQTL/'
                '{}.v8.signif_variant_gene_pairs.txt'.format(tissue),
                sep='\t', index_col=None)
            eqtls = eqtls[eqtls.gene_id.isin(genes)]
            idx = eqtls.groupby('gene_id').pval_nominal.idxmin().values
            eqtls = eqtls.loc[idx]

            eqtls.loc[:, 'chr'] = eqtls.variant_id.apply(lambda x: x.split('_')[0])
            eqtls.loc[:, 'start'] = eqtls.variant_id.apply(lambda x: int(x.split('_')[1]))
            eqtls.loc[:, 'end'] = eqtls.loc[:, 'start'] + 1
            eqtls.loc[:, 'study'] = tissue
            df = eqtls.loc[:, ['chr', 'start', 'end', 'variant_id', 'tss_distance', 'maf']]

        else:
            print('using filter: {}'.format(params.filters))
            df = pd.read_csv('output/GTEx/variant_reports/{}.all_genes.variant_report'.format(tissue), sep='\t')
            print('\t {} total records'.format(df.shape[0]))
            df = df[eval(params.filters)]
            print('\t {} remaining records'.format(df.shape[0]))

            # add tss_distance
            df.loc[:, 'tss_distance'] = df.start - df.tss

        """
        if 'eqtltop' in analysis_id:
            print('fetching top eqtls in GTEx')
            tissue_significant = pd.read_csv(
                'output/GTEx/nominally_significant_associations/{}.nominally_significant.txt'.format(tissue),
                sep='\t', index_col=None)
            gene2count = df.gene.value_counts()

            new_df = pd.concat([group.nsmallest(gene2count.get(gene, 0), 'pval_nominal')
                for gene, group in tqdm(tissue_significant.groupby('gene_id')) if gene in gene2count])

            new_df.loc[:, 'chr'] = new_df.variant_id.apply(lambda x: x.split('_')[0])
            new_df.loc[:, 'start'] = new_df.variant_id.apply(lambda x: int(x.split('_')[1]))
            new_df.loc[:, 'end'] = new_df.loc[:, 'start'] + 1
            new_df.loc[:, 'study'] = tissue

            df = new_df.loc[:, ['chr', 'start', 'end', 'variant_id', 'tss_distance', 'maf']]
        """
        # save test set of unique chr pos
        print('save test set')
        df.loc[:, 'chr_num'] = df.loc[:, 'chr'].apply(lambda x: int(x.replace('chr', '')))
        df.sort_values(by=['chr_num', 'start']).drop_duplicates(['chr', 'start'])\
            .to_csv(output.test, sep='\t', header=False, index=False)

        # put variant, gene pair into bins
        print('binning')
        df.loc[:, 'bin'] = [pair2bin(dtss, maf) for dtss, maf in zip(df.tss_distance, df.maf)]

        # count number of variants in each bin
        bins = df.bin.value_counts().to_dict()

        print('constructing background set')
        background = []
        for b, count in tqdm(list(bins.items())):
            try:
                bin_path = 'output/GTEx/maf_dtss_binned/{}.bed'.format(b)
                bin_df = pd.read_csv(bin_path, sep='\t', header=None)
                background.append(bin_df.iloc[np.random.choice(bin_df.shape[0], count*5, replace=True)])
            except Exception:
                continue
        # save background set of unique chr pos
        background_df = pd.concat(background)
        f = ~background_df.iloc[:, 4].isin(df.iloc[:, 3])
        background_df = background_df[f]
        print('removed {} variants from background'.format((~f).sum()))

        print('save background set')
        background_df.loc[:, 'chr_num'] = background_df.loc[:, 0].apply(lambda x: int(x.replace('chr', '')))
        background_df.sort_values(by=['chr_num', 1]).drop_duplicates([0, 1])\
            .to_csv(output.background, sep='\t', header=False, index=False)


rule roadmap_enrichment:
    input:
        test = 'output/GTEx/enrichment/{analysis_id}/{tissue}.test.bed',
        background = 'output/GTEx/enrichment/{analysis_id}/{tissue}.background.bed'
    output:
        'output/GTEx/enrichment/roadmap/{analysis_id}/{tissue}.roadmap.enrichments'
    run:
        import os
        import pybedtools
        import glob
        import numpy as np
        import pandas as pd
        from tqdm import tqdm
        from scipy.stats import fisher_exact
        from subprocess import check_output

        intersect_template = 'bedtools intersect -a {} -b {} -sorted| wc -l'

        def count_intersection(a, b):
            cmd = intersect_template.format(a, b)
            print(cmd)
            return int(check_output(cmd, shell=True))

        def count_lines(a):
            cmd = 'cat {} | wc -l'.format(a)
            print(cmd)
            return int(check_output(cmd, shell=True))

        def contingency_table(test_path, background_path, annot_path):
            """
            test and background DO NOT intersect
            return [[test/annotation, tesn n annotation],
                    [background/annotation, [test n annotation]]
            """
            test_in_annot = count_intersection(test_path, annot_path)
            test_not_in_annot = count_lines(test_path) - test_in_annot

            background_in_annot = count_intersection(background_path, annot_path)
            background_not_in_annot = count_lines(background_path) - background_in_annot

            return np.array([[test_in_annot, test_not_in_annot],
                    [background_in_annot, background_not_in_annot]])

        def get_record(analysis_id, tissue, eid, annotation_type):
            eid2celltype = pd.read_csv('output/annotations/roadmap/EIDlegend.txt',
                    sep='\t', index_col=0, header=None).iloc[:, 0].to_dict()
            contingency_entry_labels = np.array([['test_in_annot', 'test_not_in_annot'],
                        ['background_in_annot', 'background_not_in_annot']])

            test = 'output/GTEx/enrichment/{}/{}.test.bed'.format(analysis_id, tissue)
            background = 'output/GTEx/enrichment/{}/{}.background.bed'.format(analysis_id, tissue)
            annot_file = '{}.{}.bed'.format(eid, annotation_type)
            annotation = 'output/annotations/roadmap/{}'.format(annot_file)

            ct = np.array(contingency_table(test, background, annotation))
            odds, p = fisher_exact(ct)

            record = {a: b for a, b in zip(contingency_entry_labels.flatten(), ct.flatten())}
            record.update({
                'p': p,
                'odds_ratio': odds,
                'EID': annot_file.split('.')[0],
                'cell_type': eid2celltype.get(annot_file.split('.')[0]),
                'annotation_type': annot_file.split('.')[1],
                'tissue': tissue,
                'analysis_id': analysis_id
            })
            return record

        tissue = wildcards.tissue
        analysis_id = wildcards.analysis_id

        # get annotation file
        annotation_files = os.listdir('output/annotations/roadmap/')
        annotation_files = [x for x in annotation_files if ('enhancer' in x or 'promoter' in x)]

        records = []
        for annot_file in tqdm(annotation_files):
            eid, annotation_type, _ = annot_file.split('.')
            record = get_record(analysis_id, tissue, eid, annotation_type)
            records.append(record)

        pd.DataFrame(records).to_csv(output[0], sep='\t', index=None)
