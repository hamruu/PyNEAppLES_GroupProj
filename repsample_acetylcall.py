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
import json

from src.pyneapples import rep_sampler_2d as rstd #import the GeomReduction class from the repre_sample_2D module within the test_module package

file_path = 'acetaldehyde/acetaldehyde_source.json'

#Read the JSON file
with open(file_path, 'r') as json_file:
    acetalydehyde = json.load(json_file)

oscillator_strengths = []
excitation_energies_cm = []

for entry in acetalydehyde:
    for item in entry:
        oscillator_strengths.append(item['oscillator_strengths'])
        excitation_energies_cm.append(item['excitation_energies_cm'])

oscillator_strengths_array = np.array(oscillator_strengths)
excitation_energies_cm_array = np.array(excitation_energies_cm)

MYGEOM = rstd.GeomReduction(2000, 3, 5, 1000, 8, 16, weighted=True, pdfcomp = "KLdiv", intweights=False, verbose=False, dim1=False) #create an instance of GeomReduction with specified parameters: 500 samples, 3 states, 20 representative molecules, 100 cycles, 1 core, 1 job, without weighting, using KL divergence for PDF comparison, no integer weights, and verbose off

np.random.seed(42) #set the seed for reproducibility

MYGEOM.read_data_direct_osc(excitation_energies_cm_array, oscillator_strengths_array) #directly feed the generated data into the GeomReduction instance
MYGEOM.reduce_geoms() #start the geometry reduction process to select representative geometries
