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
    var2rs = pysam.TabixFile('output/GTEx/variantid2rsid.tsv.gz')
    pos = np.array([int(snp.split('_')[1]) for snp in snps])
    chromosome = snps[0].split('_')[0]
    records = [[
        autocast(x) for x in record.strip().split('\t')]
               for record in var2rs.fetch(chromosome, pos.min()-1, pos.max()+1)]
    records = pd.DataFrame(records)
    records.columns = var2rs.header[0].split('\t')
    snp2rsid = records[records.variant_pos.isin(pos)].set_index('variant_id').loc[:, 'rs_id_dbSNP151_GRCh38p7']
    snp2rsid = snp2rsid[snp2rsid !='.']
    return snp2rsid.to_dict()

# input snplist
snps = pd.read_csv(snakemake.input[0])
snps = snps.values.flatten()
varid2rs = make_variant2rsid(snps)

np.savetxt(snakemake.output[0], np.array(list(varid2rs.values())), fmt='%s')
json.dump(varid2rs, open(snakemake.output[1], 'w'))
