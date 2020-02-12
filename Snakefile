import itertools
import pickle
import numpy as np
import pandas as pd

configfile: "config/config.yaml"

# terminal rules
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
        )
# intermediate rules
rule get_cis_variants:
    output:
        "output/genotypes/{gene}_cis_variants"
    script:
        "workflow/scripts/get_cis_variants.py"

# simulation rules
rule simulate_single_causal_variant:
    input:
        "output/genotypes/{gene}_cis_variants"
    output:
        "output/simulation/single_causal_variant/pve_{pve}/ld_{linkage}/gene_{gene}/data"
    script:
        "workflow/scripts/single_causal_variant.py"

rule simulate_multiple_causal_variant:
    input:
        "output/genotypes/{gene}_cis_variants"
    output:
        "output/simulation/multiple_causal_variant/"
        "pve_{pve}/sparsity_{sparsity}/gene_{gene}/data"
    script:
        "workflow/scripts/multiple_causal_variant.py"

# model fitting rules
rule fit_summary_model:
    input:
        "output/simulation/{simulation}/{settings}/gene_{gene}/data"
    output:
        "output/simulation/{simulation}/{settings}/gene_{gene}/model_summary"
    wildcard_constraints:
        simulation = "(?!\/)[^\/]+(?=\/)"
    script:
        "workflow/scripts/fit_cosie_summary.py"

rule fit_genotype_model:
    input:
        "output/simulation/{simulation}/{settings}/gene_{gene}/data"
    output:
        "output/simulation/{simulation}/{settings}/gene_{gene}/model_genotype"
    wildcard_constraints:
        simulation = "(?!\/)[^\/]+(?=\/)"
    script:
        "workflow/scripts/fit_cosie_genotype.py"

rule fit_pairwise_summary_model:
    input:
        "output/simulation/single_causal_variant/pve_{pve}/ld_{linkage}/gene_{gene}/data"
    output:
        "output/simulation/{simulation}/"
        "{settings}/gene_{gene}/pairwise_summary/"
        "t1_{tissue1}_t2_{tissue2}_model_summary"
    wildcard_constraints:
        simulation = "(?!\/)[^\/]+(?=\/)"
    script:
        "workflow/scripts/fit_cosie_summary.py"

rule fit_pairwise_genotype_model:
    input:
        "output/simulation/single_causal_variant/pve_{pve}/ld_{linkage}/gene_{gene}/data"
    output:
        ("output/simulation/single_causal_variant/pve_{pve}/"
        "ld_{linkage}/gene_{gene}/pairwise_genotype/"
        "t1_{tissue1}_t2_{tissue2}_model_genotype")
    script:
        "workflow/scripts/fit_cosie_summary.py"

rule run_coloc:
    input:
        "output/simulation/{simulation}/{settings}/gene_{gene}/data"
    output:
        "output/simulation/{simulation}/{settings}/gene_{gene}/coloc"
    wildcard_constraints:
        simulation = "(?!\/)[^\/]+(?=\/)"
    script:
        "workflow/scripts/run_coloc.R"

rule format_caviar_data_ld:
    input:
        'output/{path}/data'
    output:
        ld_matrix='output/{path}/data.ld'
    run:
        data = pickle.load(open(input[0], 'rb'))
        import pdb; pdb.set_trace()
        np.savetxt(fname=output.ld_matrix, X=data['LD'], delimiter='\t')

rule format_caviar_data_zscore:
    input:
        'ouput/{path}/data'
    output:
        z_scores='output/{path}/data.z{tissue}'
    run:
        data = pickle.load(open(input[0], 'rb'))
        zscores = pd.DataFrame(data['zscores'][int(wildcards.tissue)])
        zscores.to_csv(ouput.z_scores, sep='\t', index=Falseq)

# stat gathering rules
rule make_tissue_pair_components_table:
    input:
        data_path = \
            "output/simulation/{simulation}/{settings}/gene_{gene}/data",
        genotype_model_path = \
            "output/simulation/{simulation}/{settings}/gene_{gene}/model_genotype",
        summary_model_path = \
            "output/simulation/{simulation}/{settings}/gene_{gene}/model_summary"
    output:
        genotype_output = \
            "output/simulation/{simulation}/{settings}/gene_{gene}/pairs_genotype",
        summary_output = \
            "output/simulation/{simulation}/{settings}/gene_{gene}/pairs_summary"
    wildcard_constraints:
        simulation = "(?!\/)[^\/]+(?=\/)"
    script:
        "workflow/scripts/make_tissue_pair_components_table.py"

tissue_pairs = [x for x in itertools.combinations(np.arange(7), 2)]
rule make_pairwise_pair_components_table:
    """
    same rule as above except for the pairwise outputs
    """
    input:
        data_path = "output/simulation/single_causal_variant/pve_{pve}/ld_{linkage}/gene_{gene}/data",
        summary_model_paths = expand(
            "output/simulation/single_causal_variant/pve_{pve}/ld_{linkage}/gene_{gene}/pairwise_summary/t1_{tissue1}_t2_{tissue2}_model_summary",
            tissue1=[x[0] for x in tissue_pairs],
            tissue2=[x[1] for x in tissue_pairs],
            pve='{pve}', linkage='{linkage}', gene='{gene}'
        )
    output:
        "output/simulation/single_causal_variant/pve_{pve}/ld_{linkage}/gene_{gene}/pairwise_summary/pairs_table"
    script:
        "workflow/scripts/make_pairwise_pair_components_table.py"

