import glob

chr_gene = glob.glob('/work-zfs/abattle4/karl/cosie_analysis/output/GTEx/chr2/*/*.standardized.k20.model')
chr_gene = ['/'.join(x.split('/')[7:9]) for x in chr_gene][:100]

rule simulate_n:
    input:
        expand("output/sim/n_causal_variants/{chr_gene}/sim.t{t}.n{snps_per_tissue}.pve{pve}.susie",
            chr_gene=chr_gene, t=10, snps_per_tissue=[1,2,3,4,5], pve=['05', '10', '20'])

rule ecaviar_simulate_n:
    input:
        expand("output/sim/n_causal_variants/{chr_gene}/sim.t{t}.n{snps_per_tissue}.pve{pve}.ecaviar",
            chr_gene=chr_gene[:20], t=10, snps_per_tissue=[1,2,3,4,5], pve=['05', '10', '20'])

rule coloc_simulate_n:
    input:
        expand("output/sim/n_causal_variants/{chr_gene}/sim.t{t}.n{snps_per_tissue}.pve{pve}.coloc",
            chr_gene=chr_gene[:20], t=10, snps_per_tissue=[1,2,3,4,5], pve=['05', '10', '20'])

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
        data="output/sim/{sim}/{chr}/{gene}/genotype.sim.t{t}.pve{pve}.data"
    output:
        susie="output/sim/{sim}/{chr}/{gene}/genotype.sim.t{t}.pve{pve}.susie",
    wildcard_constraints:
        gene = "[^\/]+(?=\/)"
    params:
        sample_effects=False
    script:
        "../../workflow/scripts/sim_susie.py"

rule run_susie_on_n_causal_variant_simulations:
    input:
        data="output/sim/n_causal_variants/{chr}/{gene}/sim.t{t}.n{snps_per_tissue}.pve{pve}.data"
    output:
        susie="output/sim/n_causal_variants/{chr}/{gene}/sim.t{t}.n{snps_per_tissue}.pve{pve}.susie",
    wildcard_constraints:
        snps_per_tissue = "\d+",
        gene = "[^\/]+(?=\/)"
    params:
        sample_effects=False
    script:
        "../../workflow/scripts/sim_susie.py"

rule simulate_n_causal_variant_random_effect_size:
    input:
        genotype="output/GTEx/{chr}/{gene}/{gene}.raw"
    output:
        model="output/sim/n_causal_variants/{chr}/{gene}/sim.t{t}.n{snps_per_tissue}.pve{pve}.genotype.model",
        info="output/sim/n_causal_variants/{chr}/{gene}/sim.t{t}.n{snps_per_tissue}.pve{pve}.info",
        data="output/sim/n_causal_variants/{chr}/{gene}/sim.t{t}.n{snps_per_tissue}.pve{pve}.data"
    wildcard_constraints:
        snps_per_tissue = "\d+",
        gene = "[^\/]+(?=\/)"
    params:
        sample_effects=True
    script:
        "../../workflow/scripts/sim_n_causal_variants.py"

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