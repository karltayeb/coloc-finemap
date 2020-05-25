import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
from pysam import VariantFile
import json

is_biallelic = lambda rec: len(rec.alleles) == 2
# input snplist
snps = pd.read_csv(snakemake.input[0])
snps = snps.values.flatten()

vcf = VariantFile("/work-zfs/abattle4/lab_data/dbSNP/GCF_000001405.38.gz")

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

np.savetxt(snakemake.output[0], np.array(list(varid2rs.values())), fmt='%s')
json.dump(varid2rs, open(snakemake.output[1], 'w'))