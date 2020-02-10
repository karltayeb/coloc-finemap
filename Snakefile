configfile: "config/config.yaml"

rule get_cis_variants:
    output:
        "output/genotypes/{gene}_cis_variants"
    script:
        "workflow/scripts/get_cis_variants.py"

rule simulate_single_causal_variant:
    input:
        "output/genotypes/{gene}_cis_variants"
    output:
        "output/simulation/single_causal_variant/{gene}/ld_{linkage}_pve_{pve}_data"
    script:
        "workflow/scripts/single_causal_variant.py"

rule fit_cosie_summary:
    input:
        "output/simulation/single_causal_variant/{gene}/ld_{linkage}_pve_{pve}_data"
    output:
        "output/models/cosie_summary/gene_{gene}_ld_{linkage}_pve_{pve}_model"
    script:
        "workflow/scripts/fit_cosie_summary.py"

rule fit_cosie_genotype:
    input:
        "output/simulation/single_causal_variant/{gene}/ld_{linkage}_pve_{pve}_data"
    output:
        "output/models/cosie_genotype/gene_{gene}_ld_{linkage}_pve_{pve}_model"
    script:
        "workflow/scripts/fit_cosie_genotype.py"
        
