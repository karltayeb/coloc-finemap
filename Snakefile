import itertools
import pickle
import numpy as np
import pandas as pd
import glob

configfile: "config/config.yaml"

include: 'workflow/snk/cafeh.snk.py'
include: 'workflow/snk/caviar.snk.py'
include: 'workflow/snk/coloc.snk.py'
include: 'workflow/snk/data_prep.snk.py'
include: 'workflow/snk/simulations.snk.py'
include: 'workflow/snk/enrichment.snk.py'
include: 'workflow/snk/terminal_rules.snk.py'
