import pickle
import numpy as np
from coloc.ard_ser import MVNFactorSER

def assign(obj, dictionary):
    for key in dictionary.keys():
        obj.__dict__[key] = dictionary[key]

model_path = snakemake.input[0]
component_plot_path = snakemake.output.component_plot_path
zscore_plot_path = snakemake.output.zscore_plot_path

cafeh = MVNFactorSER(np.zeros((1, 1)), np.zeros((1, 1)), 1)


assign(cafeh, pickle.load(open(model_path, 'rb')))
path = '/'.join(model_path.split('/')[:-1])
data = pickle.load(open('{}/data'.format(path), 'rb'))
cafeh.X = data['LD']
cafeh.Y = data['zscores']

cafeh.plot_components(save_path=component_plot_path)
cafeh.plot_decomposed_zscores(save_path=zscore_plot_path)