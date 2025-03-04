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

import rep_sampler_2d as rstd #import the GeomReduction class from the repre_sample_2D module within the test_module package

MYGEOM = rstd.GeomReduction(500, 3, 20, 100, 1, 1, weighted=False, pdfcomp = "KLdiv", intweights=0, verbose=False) #create an instance of GeomReduction with specified parameters: 500 samples, 3 states, 20 representative molecules, 100 cycles, 1 core, 1 job, without weighting, using KL divergence for PDF comparison, no integer weights, and verbose off

np.random.seed(42) #set the seed for reproducibility

excit_e = 2 + (np.random.rand(500, 3) - 0.5) * 4 #generate random excitation energies around 2 eV with noise, and clip them between 0 and 4 eV
excit_e = np.clip(excit_e, 0, 4)

trans_dip_mom_x = 10 + np.random.rand(500, 3) * 90 #generate random transition dipole moment components for x, y, and z with values between 10 and 100
trans_dip_mom_y = 10 + np.random.rand(500, 3) * 90
trans_dip_mom_z = 10 + np.random.rand(500, 3) * 90

MYGEOM.read_data_direct(excit_e, trans_dip_mom_x, trans_dip_mom_y, trans_dip_mom_z) #directly feed the generated data into the GeomReduction instance
MYGEOM.reduce_geoms() #start the geometry reduction process to select representative geometries
