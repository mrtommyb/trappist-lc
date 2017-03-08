from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
import sys

%matplotlib inline

from astropy.io import fits as pyfits
from astropy.stats.funcs import median_absolute_deviation as MAD
from scipy.ndimage import label

sys.path.append('/Users/tom/Projects/MontetStar/')
from martinsff import martinsff
import extract_lc

import yash_bls
from blssearch import doSearch, plotSearch, get_qf

sys.path.append('/Users/tom/gitcode/')
import dave.fileio.mastio
import dave.tpf2lc.tpf2lc as tpf2lc