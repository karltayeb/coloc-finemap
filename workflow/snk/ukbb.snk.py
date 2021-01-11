rule download_ukbb_pheno:
    input:
        manifest='output/UKBB/manifest.txt',
        pheno2manifest='output/UKBB/pheno2manifest'
    output:
        temp_save_file='output/UKBB/{phenotype}/_{phenotype}.tsv.bgz'
    run:
        import pandas as pd
        import subprocess
        import json

        manifest = pd.read_csv(input.manifest, sep='\t').set_index(['Phenotype Code', 'Sex'])
        pheno2manifest = json.load(open(input.pheno2manifest, 'r'))

        print('downloading summary stats')
        phenotype = wildcards.phenotype
        source_file = manifest.loc[:, 'wget command']\
            .loc[pheno2manifest[phenotype]]\
            .values[0]\
            .split(' ')[1]

        cmd = 'wget {} -O {}'.format(source_file, output.temp_save_file)
        print(cmd)
        subprocess.run(cmd, shell=True)

rule ukbb_build_index:
    input:
        var='output/UKBB/variants.tsv.bgz',
        sumstats='output/UKBB/{phenotype}/_{phenotype}.tsv.bgz'
    output:
        save_file='output/UKBB/{phenotype}/{phenotype}.tsv.bgz',
        tabix_index='output/UKBB/{phenotype}/{phenotype}.tsv.bgz.tbi'
    run:
        shell("paste <(zcat {input.var}) <(zcat {input.sumstats}) | bgzip > {output.save_file}")
        shell("tabix -s 2 -b 3 -e 3 -S 1 {output.save_file}")

