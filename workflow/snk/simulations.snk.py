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
    wildcard_constraints:
        gene = "[^\/]+(?=\/)"
    script:
        "workflow/scripts/multiple_causal_variant.py"


rule simulate_multiple_causal_variant:
    input:
        genotype="output/GTEx/{chr}/{gene}/{gene}.raw"
    output:
        model="output/simulation/multiple_causal/{chr}/{gene}/genotype.model",
        info="output/simulation/multiple_causal/{chr}/{gene}/sim.info",
    wildcard_constraints:
        gene = "[^\/]+(?=\/)"
    script:
        "workflow/scripts/sim_multiple_causal_variants.py"