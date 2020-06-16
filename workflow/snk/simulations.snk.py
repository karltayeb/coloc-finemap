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
            chr_gene=chr_gene, t=10, snps_per_tissue=[1,2,3,4,5], pve=['05', '10', '20'])

rule coloc_simulate_n:
    input:
        expand("output/sim/n_causal_variants/{chr_gene}/sim.t{t}.n{snps_per_tissue}.pve{pve}.coloc",
            chr_gene=chr_gene, t=10, snps_per_tissue=[1,2,3,4,5], pve=['05', '10', '20'])

rule gss_simulate_n:
    input:
        expand("output/sim/n_causal_variants/{chr_gene}/sim.t{t}.n{snps_per_tissue}.pve{pve}.gss",
            chr_gene=chr_gene, t=10, snps_per_tissue=[1,2,3,4,5], pve=['05', '10', '20'])

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

rule run_gss_on_n_causal_variant_simulations:
    input:
        data="output/sim/n_causal_variants/{chr}/{gene}/sim.t{t}.n{snps_per_tissue}.pve{pve}.data",
        info="output/sim/n_causal_variants/{chr}/{gene}/sim.t{t}.n{snps_per_tissue}.pve{pve}.info"
    output:
        model="output/sim/n_causal_variants/{chr}/{gene}/sim.t{t}.n{snps_per_tissue}.pve{pve}.gss",
    wildcard_constraints:
        snps_per_tissue = "\d+",
        gene = "[^\/]+(?=\/)"
    params:
        sample_effects=False
    script:
        "../../workflow/scripts/sim_gss.py"

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


sim_spec = pd.read_csv('output/sim/ld/sim_spec.txt', sep='\t')
def get_model(wildcards):
    return sim_spec[sim_spec.sim_id == wildcards.sim_id].source_model_path.values[0]

rule sim_from_model:
    input:
        get_model
    output:
        'output/sim/ld/{gene}/{sim_id}.sim',
    script:
        '../../workflow/scripts/sim_from_model.py'


sim_files = []
for i, row in sim_spec.iterrows():
    sim_files.append('output/sim/ld/{}/{}.sim'.format(row.gene, row.sim_id))
rule run_sim_from_model:
    input:
        expand({path}, path=sim_files)