import numpy as np
import pandas as pd
import pybedtools
from collections import namedtuple
from scipy.stats import fisher_exact
import json
import glob

import matplotlib.pyplot as plt

def contingency_table(test, background, annotation):
    """
    test and background do not intersect
    return [[test/annotation, tesn n annotation],
            [background/annotation, [test n annotation]]
    """
    test_in_annot = test.intersect(annotation).count()
    background_in_annot = background.intersect(annotation).count()
    test_size = test.count()
    background_size = background.count()
    return np.array([[test_in_annot, test_size - test_in_annot],
            [background_in_annot, background_size - background_in_annot]])

# open beds
k = snakemake.wildcards.cluster

test = pybedtools.BedTool(snakemake.input.test)

background = pybedtools.BedTool(snakemake.input.background)
background = background - test

"""
all_cset = pybedtools.BedTool('output/enrichment/component_clusters/all.sorted.merged.bed')
all_cset = all_cset - test

all_tested = pybedtools.BedTool('output/enrichment/component_clusters/all_tested_variants.merged.bed')
all_tested = all_tested - test
"""

background_sets = {
    'dtss_maf_matched_background': background,
}


annotations = {}
for f in glob.glob('/work-zfs/abattle4/marios/annotations/ct_annot/*.bed'):
    key = f.split('/')[-1].split('_')[0]
    annotations[key] = pybedtools.BedTool(f)

annotations['enhancer'] = pybedtools.BedTool(
    '/work-zfs/abattle4/marios/GTEx_v8/coloc/Enhancer_regions_cross_tissue_Roadmap_25state.bed'
)

annotations['promoter'] = pybedtools.BedTool(
    '/work-zfs/abattle4/marios/GTEx_v8/coloc/Promoter_regions_cross_tissue_Roadmap_25state.bed'
)

results = []
for annot in annotations:
    print(annot)
    annot_bed = annotations[annot]
    for background in background_sets:
        print('\t' + background)
        background_bed = background_sets[background]
        ct = contingency_table(test, background_bed, annot_bed)
        odds, p = fisher_exact(ct)
        enrichment = {
            'test_set': k,
            'background_set': background,
            'annotation': annot,
            'test_in_annot': ct[0, 0],
            'test_not_annot': ct[0, 1],
            'background_in_annot': ct[1, 0],
            'background_not_annot': ct[1, 1],
            'odds_ratio': odds,
            'p': p
        }
        results.append(enrichment)
    print('\tsaving...')
    pd.DataFrame(results).to_csv(snakemake.output[0], sep='\t')
