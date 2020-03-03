from coloc.cafeh import CAFEH
import pickle
from utils import assign

# load model
cafeh = CAFEH(**pickle.load(open(snakemake.input.model, 'rb')), K=20)
# load data
assign(cafeh, pickle.load(open(snakemake.input.data, 'rb')))

# generate credible sets
cs, p = cafeh.get_credible_sets()

