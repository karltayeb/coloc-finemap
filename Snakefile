import itertools
import pickle
import numpy as np
import pandas as pd
import glob

configfile: "config/config.yaml"

include: 'workflow/snk/cafeh.snk.py'
include: 'workflow/snk/caviar.snk.py'
include: 'workflow/snk/coloc.snk.py'
include: 'workflow/snk/data_prep.snk.py'
include: 'workflow/snk/simulations.snk.py'

# terminal rules
rule generate_figures:
    input:
        expand(
            "output/GTEx/gene_{gene}/tissue_specific_cov/summary.zscores.png", gene=config['chr22_genes']
        )

gtex_genes = np.loadtxt('/work-zfs/abattle4/karl/cosie_analysis/config/random_gene_list.txt', dtype=str)
rule run_gtex:
    input:
        expand("output/GTEx/gene_{gene}/tissue_specific_cov/model_summary", gene=gtex_genes)

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

