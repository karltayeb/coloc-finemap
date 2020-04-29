import glob

chr_gene = glob.glob('/work-zfs/abattle4/karl/cosie_analysis/output/GTEx/chr2/*/*.standardized.k20.model')
chr_gene = ['/'.join(x.split('/')[7:9]) for x in chr_gene][:200]
rule simulate_single:
    input:
        expand("output/sim/single_random/{chr_gene}/genotype.sim.t10.pve01.model", chr_gene=chr_gene),
        expand("output/sim/single_random/{chr_gene}/genotype.sim.t10.pve01.coloc", chr_gene=chr_gene),
        expand("output/sim/single_random/{chr_gene}/genotype.sim.t10.pve01.ecaviar", chr_gene=chr_gene),

        expand("output/sim/single_random/{chr_gene}/genotype.sim.t10.pve05.model", chr_gene=chr_gene),
        expand("output/sim/single_random/{chr_gene}/genotype.sim.t10.pve05.coloc", chr_gene=chr_gene),
        expand("output/sim/single_random/{chr_gene}/genotype.sim.t10.pve05.ecaviar", chr_gene=chr_gene),

        expand("output/sim/single_random/{chr_gene}/genotype.sim.t10.pve10.model", chr_gene=chr_gene),
        expand("output/sim/single_random/{chr_gene}/genotype.sim.t10.pve10.coloc", chr_gene=chr_gene),
        expand("output/sim/single_random/{chr_gene}/genotype.sim.t10.pve10.ecaviar", chr_gene=chr_gene),

        expand("output/sim/single_random/{chr_gene}/genotype.sim.t10.pve20.model", chr_gene=chr_gene),
        expand("output/sim/single_random/{chr_gene}/genotype.sim.t10.pve20.coloc", chr_gene=chr_gene),
        expand("output/sim/single_random/{chr_gene}/genotype.sim.t10.pve20.ecaviar", chr_gene=chr_gene),


        expand("output/sim/multiple_random/{chr_gene}/genotype.sim.t10.pve01.model", chr_gene=chr_gene),
        expand("output/sim/multiple_random/{chr_gene}/genotype.sim.t10.pve01.coloc", chr_gene=chr_gene),
        expand("output/sim/multiple_random/{chr_gene}/genotype.sim.t10.pve01.ecaviar", chr_gene=chr_gene),

        expand("output/sim/multiple_random/{chr_gene}/genotype.sim.t10.pve05.model", chr_gene=chr_gene),
        expand("output/sim/multiple_random/{chr_gene}/genotype.sim.t10.pve05.coloc", chr_gene=chr_gene),
        expand("output/sim/multiple_random/{chr_gene}/genotype.sim.t10.pve05.ecaviar", chr_gene=chr_gene),

        expand("output/sim/multiple_random/{chr_gene}/genotype.sim.t10.pve10.model", chr_gene=chr_gene),
        expand("output/sim/multiple_random/{chr_gene}/genotype.sim.t10.pve10.coloc", chr_gene=chr_gene),
        expand("output/sim/multiple_random/{chr_gene}/genotype.sim.t10.pve10.ecaviar", chr_gene=chr_gene),

        expand("output/sim/multiple_random/{chr_gene}/genotype.sim.t10.pve20.model", chr_gene=chr_gene),
        expand("output/sim/multiple_random/{chr_gene}/genotype.sim.t10.pve20.coloc", chr_gene=chr_gene),
        expand("output/sim/multiple_random/{chr_gene}/genotype.sim.t10.pve20.ecaviar", chr_gene=chr_gene)


rule run_susie_on_simulations:
    input:
        data="output/sim/multiple/{chr}/{gene}/genotype.sim.t{t}.pve{pve}.data"
    output:
        susie="output/sim/multiple/{chr}/{gene}/genotype.sim.t{t}.pve{pve}.susie",
    wildcard_constraints:
        gene = "[^\/]+(?=\/)"
    params:
        sample_effects=False
    script:
        "../../workflow/scripts/sim_susie.py"

rule simulate_multiple_causal_variant_fixed_effect_size:
    input:
        genotype="output/GTEx/{chr}/{gene}/{gene}.raw"
    output:
        model="output/sim/multiple/{chr}/{gene}/genotype.sim.t{t}.pve{pve}.model",
        info="output/sim/multiple/{chr}/{gene}/sim.t{t}.pve{pve}.info",
        data="output/sim/multiple/{chr}/{gene}/genotype.sim.t{t}.pve{pve}.data"
    wildcard_constraints:
        gene = "[^\/]+(?=\/)"
    params:
        sample_effects=False
    script:
        "../../workflow/scripts/sim_multiple_causal_variants.py"

rule simulate_multiple_causal_variant_random_effect_size:
    input:
        genotype="output/GTEx/{chr}/{gene}/{gene}.raw"
    output:
        model="output/sim/multiple_random/{chr}/{gene}/genotype.sim.t{t}.pve{pve}.model",
        info="output/sim/multiple_random/{chr}/{gene}/sim.t{t}.pve{pve}.info",
        data="output/sim/multiple_random/{chr}/{gene}/genotype.sim.t{t}.pve{pve}.data"
    wildcard_constraints:
        gene = "[^\/]+(?=\/)"
    params:
        sample_effects=True
    script:
        "../../workflow/scripts/sim_multiple_causal_variants.py"

rule simulate_single_causal_variant_fixed_effect_size:
    input:
        genotype="output/GTEx/{chr}/{gene}/{gene}.raw"
    output:
        model="output/sim/single/{chr}/{gene}/genotype.sim.t{t}.pve{pve}.model",
        info="output/sim/single/{chr}/{gene}/sim.t{t}.pve{pve}.info",
        data="output/sim/single/{chr}/{gene}/genotype.sim.t{t}.pve{pve}.data"
    wildcard_constraints:
        gene = "[^\/]+(?=\/)"
    params:
        sample_effects=False
    script:
        "../../workflow/scripts/sim_single_causal_variant.py"


rule simulate_single_causal_variant_random_effect_size:
    input:
        genotype="output/GTEx/{chr}/{gene}/{gene}.raw"
    output:
        model="output/sim/single_random/{chr}/{gene}/genotype.sim.t{t}.pve{pve}.model",
        info="output/sim/single_random/{chr}/{gene}/sim.t{t}.pve{pve}.info",
        data="output/sim/single_random/{chr}/{gene}/genotype.sim.t{t}.pve{pve}.data"
    wildcard_constraints:
        gene = "[^\/]+(?=\/)"
    params:
        sample_effects=True
    script:
        "../../workflow/scripts/sim_single_causal_variant.py"