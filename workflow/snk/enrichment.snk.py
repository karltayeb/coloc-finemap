rule get_tissue_expressed_genes_bed:
    output:
        temp('output/GTEx/enrichment/{tissue}.genes.bed')
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

rule get_tissue_variant_gene_bank:
    input:
        'output/GTEx/enrichment/{tissue}.genes.bed'
    output:
        'output/GTEx/enrichment/bank/{tissue}.bank.bed'
    run:
        import subprocess
        cmd = "bedtools closest -a output/GTEx/GTEx.afreq.ldscore.bed -b {input} -d "\
            "| awk \' BEGIN {{FS = \"\\t\"}}; {{print $1, $2, $3, $4, $5, $10, $12, $13, $14, $15,$15* $16}}\' > {output}"
        cmd = cmd.format(input=input[0], output=output[0])
        print(cmd)
        subprocess.run(cmd, shell=True)

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
        df = pd.read_csv(output[0], sep=' ', usecols=[0, 1, 2, 3, 4, 5, 10], header=None)
        df.columns = np.array(['chr', 'start', 'end', 'variant_id', 'maf', 'ldscore', 'tss'])
        df = df[(df.tss > -1e6) & (df.tss < 1e6)]

        print('binning')
        digitize_column(df, 'maf', 25)
        digitize_column(df, 'ldscore', 10)
        digitize_column(df, 'tss', 40)
        df.to_csv(output[0], sep='\t', index=False)

rule gtex_make_test_set:
    input:
        'output/GTEx/variant_reports/{method}/{tissue}.all_genes.variant_report',
        bank = 'output/GTEx/enrichment/bank/{tissue}.bank.bed'
    output:
        test = temp('output/GTEx/enrichment/{method}/{analysis_id}/{tissue}.test_temp.bed'),
        test_binned = 'output/GTEx/enrichment/{method}/{analysis_id}/{tissue}.test.bed',
    params:
        filters = lambda wildcards: config['enrichment_filters'][wildcards.analysis_id]
    group: "tissue_analysis"
    run:
        tissue = wildcards.tissue
        analysis_id = wildcards.analysis_id
        method = wildcards.method
        import pandas as pd
        import subprocess

        if method != 'eqtl':
            df = pd.read_csv(input[0], sep='\t')
            df = df[eval(params.filters)]
        else:
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
            df = eqtls.loc[:, ['chr', 'start', 'end', 'variant_id', 'variant_id', 'study', 'gene_id']]

        # save temp test file
        df.iloc[:, :7].sort_values(['chr', 'start']).to_csv(output.test, sep='\t', index=False, header=False)
        
        # use bedtools to make test file with bin info
        cmd = 'bedtools intersect -a {} -b {} -wa -wb -sorted > {}'.format(output.test, input.bank, output.test_binned)
        print(cmd)
        subprocess.run(cmd, shell=True)


rule gtex_make_background_set:
    input:
        test = 'output/GTEx/enrichment/{method}/{analysis_id}/{tissue}.test.bed',
        bank = 'output/GTEx/enrichment/bank/{tissue}.bank.bed'
    output:
        background = 'output/GTEx/enrichment/{method}/{analysis_id}/{tissue}.background.bed'
    group: "tissue_analysis"
    run:
        import pandas as pd
        import tqdm

        test = pd.read_csv(input.test, sep='\t', usecols=[0,1,2,3,14,16,18], header=None)
        test.columns = ['chr', 'start', 'end', 'variant_id', 'maf_bin', 'ldscore_bin', 'dtss_bin']

        bank = pd.read_csv(input.bank, usecols=[0,1,2,3,7,9,11],  sep='\t', header=None)
        bank.columns = ['chr', 'start', 'end', 'variant_id', 'maf_bin', 'ldscore_bin', 'dtss_bin']

        bin2count = test.groupby(['maf_bin', 'ldscore_bin', 'dtss_bin']).agg(['count']).iloc[:, 0].to_dict()

        background = []
        for key, grp in tqdm.tqdm(bank.groupby(['maf_bin', 'ldscore_bin', 'dtss_bin'])):
            background.append(grp.sample(5*bin2count.get(key, 0)))

        background = pd.concat(background)
        background = background.sort_values(['chr', 'start'])
        background = background[~background.variant_id.isin(test.variant_id)]
        background.to_csv(output.background, sep='\t', index=False, header=False)


rule roadmap_enrichment:
    input:
        test = 'output/GTEx/enrichment/{method}/{analysis_id}/{tissue}.test.bed',
        background = 'output/GTEx/enrichment/{method}/{analysis_id}/{tissue}.background.bed'
    output:
        'output/GTEx/enrichment/{method}/roadmap/{analysis_id}/{tissue}.roadmap.enrichments'
    group: "tissue_analysis"
    run:
        import os
        import pybedtools
        import glob
        import numpy as np
        import pandas as pd
        from tqdm import tqdm
        from scipy.stats import fisher_exact
        from subprocess import check_output

        intersect_template = 'bedtools intersect -a {} -b {} -sorted | wc -l'

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

        def get_record(test, background, tissue, eid, annotation_type):
            eid2celltype = pd.read_csv('output/annotations/roadmap/EIDlegend.txt',
                    sep='\t', index_col=0, header=None).iloc[:, 0].to_dict()
            contingency_entry_labels = np.array([['test_in_annot', 'test_not_in_annot'],
                        ['background_in_annot', 'background_not_in_annot']])
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
                'analysis_id': analysis_id,
                'method': method
            })
            return record

        tissue = wildcards.tissue
        analysis_id = wildcards.analysis_id
        method = wildcards.method
        
        # get annotation file
        annotation_files = os.listdir('output/annotations/roadmap/')
        annotation_files = [x for x in annotation_files if ('enhancer' in x or 'promoter' in x)]

        records = []
        for annot_file in tqdm(annotation_files):
            eid, annotation_type, _ = annot_file.split('.')
            record = get_record(input[0], input[1], tissue, eid, annotation_type)
            records.append(record)

        pd.DataFrame(records).to_csv(output[0], sep='\t', index=None)
