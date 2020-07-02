import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import pysam
import json
import ast

def autocast(s):
    try:
        return ast.literal_eval(s)
    except Exception:
        return s

def make_variant2rsid(snps):
    """
    take a list of snps
    use GTEx variantid 2 rsid map to create a dictionary
    {variant_id: rsid}
    """
    var2rs = pysam.TabixFile('../../output/GTEx/variantid2rsid.tsv.gz')
    pos = np.array([int(snp.split('_')[1]) for snp in snps])
    chromosome = snps[0].split('_')[0]
    records = [[
        autocast(x) for x in record.strip().split('\t')]
        for record in var2rs.fetch(chromosome, pos.min()-1, pos.max()+1)]
    records = pd.DataFrame(records)
    records.columns = var2rs.header[0].split('\t')
    snp2rsid = records[records.variant_pos.isin(pos)].set_index('variant_id').loc[:, 'rs_id_dbSNP151_GRCh38p7']
    snp2rsid = snp2rsid[snp2rsid!='.']
    return snp2rsid.to_dict()


def make_variant2rsid_old(snps):
    is_biallelic = lambda rec: len(rec.alleles) == 2
    vcf = pysam.VariantFile("/work-zfs/abattle4/lab_data/dbSNP/GCF_000001405.38.gz")

    # get chromosome and linear position
    pos = np.array([int(snp.split('_')[1]) for snp in snps])
    chr_idx = int(snps[0].split('_')[0][3:]) - 1
    chromosome = list(vcf.header.contigs)[chr_idx]

    varid2posrefalt = {}
    for snp in snps:
        varid2posrefalt[snp] = tuple(snp.split('_')[1:4])

    posrefalt2rs = {}
    pos2rs = defaultdict(list)
    rs2rec = {}

    for rec in vcf.fetch(chromosome, pos.min(), pos.max()):
        if is_biallelic(rec) and (rec.pos in pos):
            posrefalt2rs[(str(rec.pos), *rec.alleles[:2])] = rec.id
            rs2rec[rec.id] = rec
            pos2rs[rec.pos].append(rec.id)

    print('{} biallelic variants'.format(len(posrefalt2rs)))

    varid2rs = {}
    missed = []
    for snp in snps:
        try:
            varid2rs[snp] = posrefalt2rs[varid2posrefalt[snp]]
        except Exception:
            missed.append(tuple(snp.split('_')[1:4]))
    print('{} biallelic variants mapped to rsid'.format(len(varid2rs)))


# input snplist
snps = pd.read_csv(snakemake.input[0])
snps = snps.values.flatten()
varid2rs = make_variant2rsid(snps)

np.savetxt(snakemake.output[0], np.array(list(varid2rs.values())), fmt='%s')
json.dump(varid2rs, open(snakemake.output[1], 'w'))
