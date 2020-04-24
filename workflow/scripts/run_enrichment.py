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
    test_in_annot = (test + annotation).count()
    background_in_annot = (background + annotation).count()
    test_size = test.count()
    background_size = background.count()
    return np.array(
        [[test_in_annot, test_size - test_in_annot],
        [background_in_annot, background_size - background_in_annot]]
    )

# open beds
k = snakemake.wildcards.group

test = pybedtools.BedTool(snakemake.input.test)

background = pybedtools.BedTool(snakemake.input.background)
background = background - test

eqtltop = pybedtools.BedTool(snakemake.input.eqtltop)
eqtltop = eqtltop

results = []
promoter_paths = glob.glob('output/enrichment/annotations/ct_annot/*.enhancer.bed')
for annot_path in promoter_paths:
    annot_label = annot_path.split('/')[-1][:-4]
    print(annot_label)

    annot = pybedtools.BedTool(annot_path)
    ct = contingency_table(test, background, annot)
    odds, p = fisher_exact(ct)
    enrichment = {
        'test_set': k,
        'background_set': 'dtss_maf_matched_background',
        'annotation': annot_label,
        'test_in_annot': ct[0, 0],
        'test_not_annot': ct[0, 1],
        'background_in_annot': ct[1, 0],
        'background_not_annot': ct[1, 1],
        'odds_ratio': odds,
        'p': p
    }
    results.append(enrichment)

promoter_paths = glob.glob('output/enrichment/annotations/ct_annot/all.promoter.bed')
for annot_path in promoter_paths:
    annot_label = annot_path.split('/')[-1][:-4]
    print(annot_label)

    annot = pybedtools.BedTool(annot_path)
    ct = contingency_table(test, background, annot)
    odds, p = fisher_exact(ct)
    enrichment = {
        'test_set': k,
        'background_set': 'dtss_maf_matched_background',
        'annotation': annot_label,
        'test_in_annot': ct[0, 0],
        'test_not_annot': ct[0, 1],
        'background_in_annot': ct[1, 0],
        'background_not_annot': ct[1, 1],
        'odds_ratio': odds,
        'p': p
    }
    results.append(enrichment)

pd.DataFrame(results).to_csv(snakemake.output[0], sep='\t')

