import glob

chr_gene = glob.glob('/work-zfs/abattle4/karl/cosie_analysis/output/GTEx/chr2/*/*.standardized.k20.model')
chr_gene = ['/'.join(x.split('/')[7:9]) for x in chr_gene][:200]
rule simulate:
    input:
        expand("output/sim/multiple/{chr_gene}/genotype.sim.t20.pve01.model", chr_gene=chr_gene),
        expand("output/sim/multiple/{chr_gene}/genotype.sim.t20.pve05.model", chr_gene=chr_gene),
        expand("output/sim/multiple/{chr_gene}/genotype.sim.t20.pve10.model", chr_gene=chr_gene),
        expand("output/sim/multiple/{chr_gene}/genotype.sim.t20.pve20.model", chr_gene=chr_gene)


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


rule simulate_single:
    input:
        expand("output/sim/single/{chr_gene}/genotype.sim.t10.pve01.model", chr_gene=chr_gene),
        expand("output/sim/single/{chr_gene}/genotype.sim.t10.pve05.model", chr_gene=chr_gene),
        expand("output/sim/single/{chr_gene}/genotype.sim.t10.pve10.model", chr_gene=chr_gene),
        expand("output/sim/single/{chr_gene}/genotype.sim.t10.pve20.model", chr_gene=chr_gene)

rule simulate_single_causal_variant2:
    input:
        genotype="output/GTEx/{chr}/{gene}/{gene}.raw"
    output:
        model="output/sim/single/{chr}/{gene}/genotype.sim.t{t}.pve{pve}.model",
        info="output/sim/single/{chr}/{gene}/sim.t{t}.pve{pve}.info",
        data="output/sim/single/{chr}/{gene}/genotype.sim.t{t}.pve{pve}.data"

    wildcard_constraints:
        gene = "[^\/]+(?=\/)"
    script:
        "../../workflow/scripts/sim_single_causal_variants.py"