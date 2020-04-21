import glob

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

chr_gene = glob.glob('/work-zfs/abattle4/karl/cosie_analysis/output/GTEx/chr2/*/*.standardized.k20.model')
chr_gene = ['/'.join(x.split('/')[7:9]) for x in chr_gene][:500]
rule simulate:
    input:
        expand("output/sim/multiple/{chr_gene}/genotype.sim.t20.pve01.model", chr_gene=chr_gene),
        expand("output/sim/multiple/{chr_gene}/genotype.sim.t20.pve05.model", chr_gene=chr_gene),
        expand("output/sim/multiple/{chr_gene}/genotype.sim.t20.pve10.model", chr_gene=chr_gene)

rule simulate_multiple_causal_variant2:
    input:
        genotype="output/GTEx/{chr}/{gene}/{gene}.raw"
    output:
        model="output/sim/multiple/{chr}/{gene}/genotype.sim.t{t}.pve{pve}.model",
        susie="output/sim/multiple/{chr}/{gene}/genotype.sim.t{t}.pve{pve}.susie",
        info="output/sim/multiple/{chr}/{gene}/sim.t{t}.pve{pve}.info",
        data="output/sim/multiple/{chr}/{gene}/genotype.sim.t{t}.pve{pve}.data"

    wildcard_constraints:
        gene = "[^\/]+(?=\/)"
    script:
        "../../workflow/scripts/sim_multiple_causal_variants.py"