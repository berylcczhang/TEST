#
#This script displays the output of the test problem: example
#See section 10 of the slug manual 
#

import numpy as np
import matplotlib.pyplot as plt
import os
import os.path as osp
import copy

try:
    from slugpy import *    # If slugpy is already in our path
except ImportError:
    # If import failed, try to find slugpy in $SLUG_DIR
    if 'SLUG_DIR' in os.environ:
        import sys
        cur_path = copy.deepcopy(sys.path)
        sys.path.append(os.environ['SLUG_DIR'])
        from slugpy import *
        sys.path = cur_path
    else:
        raise ImportError("No module named slugpy")

# Get slug dir
if 'SLUG_DIR' in os.environ:
    slugdir = os.environ['SLUG_DIR']
else:
    slugdir = '.'


# Read runs
