# /// script
# requires-python = ">=3.13"
# dependencies = [
#     "joblib",
#     "matplotlib",
#     "numpy",
#     "scipy",
# ]
# ///

from tabnanny import verbose

import numpy as np

import test_module.repre_sample_2D as rstd

MYGEOM = rstd.GeomReduction(500, 3, 10, 10, 1, 1, weighted=False, pdfcomp = "KLdiv", intweights=0, verbose=False)

# Set the seed for reproducibility
np.random.seed(42)
# Generate random float values around 2 with some noise
excit_e = 2 + (np.random.rand(3, 500) - 0.5) * 4
excit_e = np.clip(excit_e, 0, 4)
trans_dip_mom = 10 + np.random.rand(3, 500, 3) * 90

MYGEOM.read_data_direct(excit_e, trans_dip_mom)
MYGEOM.reduce_geoms()
