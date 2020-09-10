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

        print('using filter: {}'.format(params.filters))
        df = pd.read_csv('output/GTEx/variant_reports/{}.all_genes.variant_report'.format(tissue), sep='\t')
        print('\t {} total records'.format(df.shape[0]))
        df = df[eval(params.filters)]
        print('\t {} remaining records'.format(df.shape[0]))
        # put variant, gene pair into bins
        df.loc[:, 'dtss'] = df.start - df.tss
        df.loc[:, 'bin'] = [pair2bin(dtss, maf) for dtss, maf in zip(df.dtss, df.maf)]

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

        # save test set of unique chr pos
        df.loc[:, 'chr_num'] = df.loc[:, 'chr'].apply(lambda x: int(x.replace('chr', '')))
        df.sort_values(by=['chr_num', 'start']).drop_duplicates(['chr', 'start'])\
            .to_csv(output.test, sep='\t', header=False, index=False)

        # save background set of unique chr pos
        background_df = pd.concat(background)
        f = ~background_df.iloc[:, 4].isin(df.iloc[:, 3])
        print('removing {} variants from background'.format((~f).sum()))
        background_df = background_df[f]
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

        def contingency_table(test_path, background_path, annot_path):
            """
            test and background DO NOT intersect
            return [[test/annotation, tesn n annotation],
                    [background/annotation, [test n annotation]]
            """
            intersect_template = 'bedtools intersect -a {} -b {} -sorted| wc -l'
            count_intersection = lambda a, b: int(check_output(intersect_template.format(a, b), shell=True))
            count_lines = lambda a: int(check_output('cat {} | wc -l'.format(a), shell=True))

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

            test = pybedtools.BedTool('output/GTEx/enrichment/{}/{}.test.bed'.format(analysis_id, tissue))
            background = pybedtools.BedTool('output/GTEx/enrichment/{}/{}.background.bed'.format(analysis_id, tissue))
            annot_file = '{}.{}.bed'.format(eid, annotation_type)
            annotation_path = 'output/annotations/roadmap/{}'.format(annot_file)
            annotation = pybedtools.BedTool(annotation_path)

            ct = np.array(contingency_table(test, background-test, annotation))
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
