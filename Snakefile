import itertools
import pickle
import numpy as np
import pandas as pd
import glob
from coloc.misc import repair_model

configfile: "config/config.yaml"

include: 'workflow/snk/cafeh.snk.py'
include: 'workflow/snk/caviar.snk.py'
include: 'workflow/snk/coloc.snk.py'
include: 'workflow/snk/data_prep.snk.py'
include: 'workflow/snk/simulations.snk.py'
include: 'workflow/snk/enrichment.snk.py'

gtex_genes = np.loadtxt('/work-zfs/abattle4/karl/cosie_analysis/config/random_gene_list.txt', dtype=str)
# terminal rules
rule generate_figures:
    input:
        expand(
            "output/GTEx/gene_{gene}/tissue_specific_cov/summary.zscores.png", gene=gtex_genes
        )

k20_genes = ['/'.join(x.split('/')[:-1]) for x in glob.glob('output/GTEx/chr[2-8]/*/genotype.standardized.k20.model')]
rule repair_k20_models:
    input:
        expand(
            "{path}/genotype.standardized.k20.repaired.log", path=k20_genes
        )

g_to_gss = pd.read_csv('output/GTEx/g_to_gss.txt').values.flatten()
rule fit_gss_from_g:
    input:
        list(g_to_gss[:1000])

rule run_component_cluster_enrichment:
    input:
        test='{path}/{group}.merged.bed',
        background='{path}/{group}.background.merged.bed'
    output:
        "{path}/{group}.enrichment"
    script:
        "workflow/scripts/run_enrichment.py"

tissues = [x.split('/')[-1].split('.')[0] for x in
    glob.glob('output/enrichment/tissue_specific_components/*.background.merged.bed')]

rule run_tissue_component_enrichment:
    input:
        expand('output/enrichment/tissue_components2/{tissue}.enrichment', tissue=tissues),
        expand('output/enrichment/eqtltop/{tissue}.eqtltop2.enrichment', tissue=tissues),

rule run_tissue_component_enrichment2:
    input:
        expand('output/enrichment/tissue_components/{tissue}.enrichment', tissue=tissues)

rule repair_k20_model:
    input:
        '{path}/genotype.standardized.k20.model'
    output:
        '{path}/genotype.standardized.k20.repaired.log'
    run:
        repair_model(input[0])
        print('model repaired', file=open(output[0], 'w'))

#run_genes = [x.split('/')[-2] for x in glob.glob('output/GTEx/*/*/genotype.model')]
run_genes = np.loadtxt('output/GTEx/run_genes.txt', dtype=str)


#k20_genes = ['/'.join(x.split('/')[:-1]) for x in glob.glob('output/GTEx/*/*/genotype.k20.model')]
rule generate_genotype_reportsk20:
    input:
        expand(
            "{path}/genotype.k20.csets", path=k20_genes
        )

#k20_genes = ['/'.join(x.split('/')[:-1]) for x in glob.glob('output/GTEx/*/*/genotype.standardized.k20.model')]
rule generate_standardized_genotype_reportsk20:
    input:
        expand(
            "{path}/genotype.standardized.k20.csets", path=k20_genes
        )


rule generate_genotype_reports:
    input:
        expand(
            "{path}/genotype.csets", path=run_genes
        )

rule run_gtex:
    input:
        expand("output/GTEx/gene_{gene}/tissue_specific_cov/model_summary", gene=gtex_genes)

rule run_gtex_genotype:
    input:
        expand("output/GTEx/gene_{gene}/tissue_specific_cov/model_genotype", gene=gtex_genes)

rule all_tissue_pairs:
    input:
        expand(
            "output/simulation/single_causal_variant/pve_{pve}/ld_{linkage}/gene_{gene}/pairs_summary",
            pve=config["pves"], linkage=config["linkages"], gene=config["genes"]
        )

rule all_pairwise_pairs:
    input:
        expand(
            ("output/simulation/single_causal_variant/pve_{pve}/"
            "ld_{linkage}/gene_{gene}/pairwise_summary/pairs_table"),
            pve=config['pves'], linkage=config['linkages'],
            gene=config['genes']
        )

rule all_coloc:
    input:
        expand(
            "output/simulation/single_causal_variant/pve_{pve}/ld_{linkage}/gene_{gene}/coloc",
            pve=config["pves"], linkage=config["linkages"], gene=config["genes"]
        )

rule run_multiple_causal_variant_simulation:
    input:
        expand(
            "output/simulation/multiple_causal_variant/pve_{pve}/sparsity_0.2/"
            "gene_{gene}/pairs_summary",
            pve=config['pves'], gene=config['genes']
        ),
        expand(
            "output/simulation/multiple_causal_variant/pve_{pve}/sparsity_0.2/"
            "gene_{gene}/coloc",
            pve=config['pves'], gene=config['genes']
        ),
        expand(
            "output/simulation/multiple_causal_variant/pve_{pve}/sparsity_0.2/"
            "gene_{gene}/ecaviar",
            pve=config['pves'], gene=config['genes']
        ),
        expand(
            "output/simulation/multiple_causal_variant/pve_{pve}/sparsity_0.2/"
            "gene_{gene}/max_min_variance_summary",
            pve=config['pves'], gene=config['genes']
        )

rule run_multiple_causal_variant_tissue_specific_cov:
    input:
        expand(
            "output/simulation/multiple_causal_variant/pve_{pve}/sparsity_0.2/"
            "gene_{gene}/tissue_specific_cov/ecaviar",
            pve=config['pves'], gene=config['genes']
        ),
        expand(
            "output/simulation/multiple_causal_variant/pve_{pve}/sparsity_0.2/"
            "gene_{gene}/tissue_specific_cov/max_min_variance_summary",
            pve=config['pves'], gene=config['genes']
        )

rule run_single_causal_variant_simulation:
    input:
        expand(
            "output/simulation/single_causal_variant/"
            "pve_{pve}/ld_{linkage}/gene_{gene}/pairs_summary",
            pve=config["pves"], linkage=config["linkages"], gene=config["genes"]
        ),
        expand(
            "output/simulation/single_causal_variant/"
            "pve_{pve}/ld_{linkage}/gene_{gene}/coloc",
            pve=config["pves"], linkage=config["linkages"], gene=config["genes"]
        ),
        expand(
            "output/simulation/single_causal_variant/"
            "pve_{pve}/ld_{linkage}/gene_{gene}/ecaviar",
            pve=config["pves"], linkage=config["linkages"], gene=config["genes"]
        ),
        expand(
            "output/simulation/single_causal_variant/"
            "pve_{pve}/ld_{linkage}/gene_{gene}/max_min_variance_summary",
            pve=config["pves"], linkage=config["linkages"], gene=config["genes"]
        )

rule run_chr22_cafeh_summary:
    input:
        expand(
            "output/GTEx/gene_{gene}/model_summary", gene=config['chr22_genes']
        )

rule run_chr22_cafeh_summary_ts_cov:
    input:
        expand(
            "output/GTEx/gene_{gene}/tissue_specific_cov/model_summary", gene=config['chr22_genes']
        )

rule run_chr22_cafeh_summary_rg_cov:
    input:
        expand(
            "output/GTEx/gene_{gene}/regressed_genotype_cov/model_summary", gene=config['chr22_genes']
        )

rule run_chr22_cafeh_genotype:
    input:
        expand(
            "output/GTEx/gene_{gene}/model_genotype", gene=config['chr22_genes']
        )

