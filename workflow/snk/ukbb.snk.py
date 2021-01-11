rule download_ukbb_pheno:
    input:
        manifest='output/UKBB/manifest.txt',
        pheno2manifest='output/UKBB/pheno2manifest'
    output:
        temp_sumstats='output/UKBB/{phenotype}/_{phenotype}.tsv.bgz'
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

        cmd = 'wget {} -O {}'.format(source_file, output.temp_sumstats)
        print(cmd)
        subprocess.run(cmd, shell=True)

rule ukbb_build_index:
    input:
        var='output/UKBB/variants.tsv.bgz',
        sumstats_raw='output/UKBB/{phenotype}/_{phenotype}.tsv.bgz'
    output:
        sumstats='output/UKBB/{phenotype}/{phenotype}.tsv.bgz',
        tabix_index='output/UKBB/{phenotype}/{phenotype}.tsv.bgz.tbi'
    run:
        shell("paste <(zcat {input.var}) <(zcat {input.sumstats_raw}) | bgzip > {output.sumstats}")
        shell("tabix -s 2 -b 3 -e 3 -S 1 {output.sumstats}")


rule ukbb_get_hits:
    input:
        sumstats='output/UKBB/{phenotype}/{phenotype}.tsv.bgz'
    output:
        hits='output/UKBB/{phenotype}/{phenotype}.hits.txt'
    run:
        shell("zcat {input.sumstats} | awk '{{if($NF < 1e-5){{print}}}}' > {output.hits}")

