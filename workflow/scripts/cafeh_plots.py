import pickle
import numpy as np
from coloc.cafeh import CAFEH

def assign(obj, dictionary):
    for key in dictionary.keys():
        obj.__dict__[key] = dictionary[key]

model_path = snakemake.input.model_path
data_path = snakemake.input.data_path
component_plot_path = snakemake.output.component_plot_path
zscore_plot_path = snakemake.output.zscore_plot_path

cafeh = CAFEH(np.zeros((1, 1)), np.zeros((1, 1)), 1)


assign(cafeh, pickle.load(open(model_path, 'rb')))
path = '/'.join(model_path.split('/')[:-1])
data = pickle.load(open(data_path, 'rb'))
cafeh.X = data['LD']
cafeh.Y = data['zscores']

cafeh.plot_components(save_path=component_plot_path)
cafeh.plot_decomposed_zscores(save_path=zscore_plot_path)