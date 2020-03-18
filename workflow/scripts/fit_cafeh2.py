import numpy as np
import pandas as pd
import glob
import pickle
from coloc.cafeh2 import CAFEH2
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib

def make_mask(ld):
    excise = []
    for i, row in enumerate(ld):
        r = row[i+1:]
        excise.append(list((np.arange(r.size) + i + 1)[np.abs(r) > 0.95]))
    excise = np.unique(np.concatenate(excise))
    mask = ~np.isin(np.arange(ld.shape[0]), excise)
    return mask

def make_bed(cafeh, output_dir, gene):
    cs, p = cafeh.get_credible_sets()
    bed = []
    for k in range(20):
        if p[k] > 0.9:
            for variant in cs[k]:
                chom, pos = variant.split('_')[:2]
                pos = int(pos)
                line = '{}\t{}\t{}\t{}'.format(chom, pos, pos+1, k)
                bed.append(line)
    with open(output_dir + '{}.credible_sets'.format(gene), 'w') as f:
        f.write('\n'.join(bed))

def assign(obj, dictionary):
    for key in dictionary.keys():
        obj.__dict__[key] = dictionary[key]

data = pickle.load(open(snakemake.input[0], 'rb'))
ld = data['X']
cov = np.cov(data['Y'].T)
gene = path.split('/')[-2].split('_')[1]
dat2 = {
    'R1': ld,
    'R2': ld + np.diag(np.diag(cov)),
    'Y': data['Y'],
    'snp_ids': data['snp_ids'],
    'tissue_ids': data['tissue_ids']
}

cafeh2 = CAFEH2(**dat2, K=20)
cafeh2.fit(max_iter=500, update_active=False, update_weights=True,
          update_pi=True, ARD_weights=True, verbose=True)

path = '/'.join(snakemake.output[0].split('/')[:-1])
name = snakemake.output[0].split('/')[-1]
# save model
cafeh2.save(path, name, save_data=False)

# make credible set bed
make_bed(cafeh2, path, gene)
print('gene: {}, converged: {}, iters: {}, ELBO: {}'.format(
    gene, len(cafeh2.elbos) < 500, len(cafeh2.elbos), cafeh2.elbos[-1]
))