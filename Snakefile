import itertools
import pickle
import numpy as np
import pandas as pd
import glob
from coloc.misc import *

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

k20_genes = ['/'.join(x.split('/')[:-1]) for x in glob.glob('output/GTEx/chr[2-8]/*/genotype.k20.model')]
rule repair_k20_models:
    input:
        expand(
            "{path}/genotype.standardized.k20.repaired.log", path=k20_genes
        )

rule repair_k20_model:
    input:
        '{path}/genotype.standardized.k20.model'
    output:
        '{path}/genotype.standardized.k20.repaired.log'
    script:
        repair_model(input[0])
        print('reapired_model', f=open(output[0], 'w'))

#run_genes = [x.split('/')[-2] for x in glob.glob('output/GTEx/*/*/genotype.model')]
run_genes = np.loadtxt('output/GTEx/run_genes.txt', dtype=str)


k20_genes = ['/'.join(x.split('/')[:-1]) for x in glob.glob('output/GTEx/*/*/genotype.k20.model')]
rule generate_genotype_reportsk20:
    input:
        expand(
            "{path}/genotype.k20.csets", path=k20_genes
        )

k20_genes = ['/'.join(x.split('/')[:-1]) for x in glob.glob('output/GTEx/*/*/genotype.standardized.k20.model')]
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

