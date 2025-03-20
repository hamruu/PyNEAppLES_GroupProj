# Acetaldehyde Module: Example Usage

This folder contains all the necessary data to reproduce an exemplar usage of the module for the photoabsorption cross-section calculation of acetaldehyde.

## Files and Their Descriptions

### `optfreq_mp2_vtz.out`
- Output file from an ORCA 5 optimization frequency calculation of acetaldehyde.
- This file was inputted into the [**Harmonwig**](https://github.com/ispg-group/harmonwig) program to generate a Wigner sample of geometries.

### `harmonix_samples.xyz`
- Contains 2000 Wigner-sampled acetaldehyde geometries.
- This file serves as an input for the **AtmoSpec** program.

### `acetaldehyde_course.json`
- A `.json` file containing oscillator strengths and excitation energies for each geometry.
- Extracted from **AtmoSpec**.

### `subsamples/subsamples_1state`
- Directory containing `.xyz` files of subset geometries.
- These geometries were prepared using the **PyNEAppLES** module in `repsample_acetylcall.py` and are ready for input back into **AtmoSpec**.

### `Results_1state` and `Results_3states`
- Contain `.tsv` files with the resultant photoabsorption cross-sections.
  - **High Theory Level:** EOM-CCSD.
  - **Low Theory Level:** TDA/TDDFT.
- These results are provided for each tested geometry count.

---

## Workflow Overview

1. **Geometry Optimization:** Start with the `optfreq_mp2_vtz.out` output file from ORCA 5.
2. **Wigner Sampling:** Use the **Harmonwig** program to create a Wigner sample (`harmonix_samples.xyz`).
3. **Spectral Analysis:** Input Wigner-sampled geometries into **AtmoSpec** to generate oscillator strength data (`acetaldehyde_course.json`).
4. **Subset Selection:** Run `repsample_acetylcall.py` to select representative subsets of geometries.
5. **Cross-Section Calculation:** Analyze the subsets using **AtmoSpec** at both EOM-CCSD and TDA/TDDFT theory levels, saving results into `.tsv` files.
