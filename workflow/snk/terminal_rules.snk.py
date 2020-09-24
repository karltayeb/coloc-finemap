import numpy as np
import json
import glob
import pandas as pd
from collections import defaultdict
from os.path import isfile

def get_paths(request):
    """
    get list of tile from request file
    """
    paths = open(request, 'r').read().strip().split('\n')
    paths = [x for x in paths if not isfile(x)]
    return paths

rule terminal_rule:
    input:
        expand('{path}', path=get_paths(config['request']))
