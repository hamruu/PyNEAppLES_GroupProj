import numpy as np
import pytest
from src.pyneapples.rep_sampler_2d import GeomReduction

@pytest.fixture
def geom_reduction_instance():
    np.random.seed(42)
    gr = GeomReduction(500, 3, 20, 100, 1, 1, weighted=False, pdfcomp="KLdiv", intweights=0, verbose=False)
    
    excit_e = 2 + (np.random.rand(500, 3) - 0.5) * 4
    excit_e = np.clip(excit_e, 0, 4)
    
    trans_dip_mom_x = 10 + np.random.rand(500, 3) * 90
    trans_dip_mom_y = 10 + np.random.rand(500, 3) * 90
    trans_dip_mom_z = 10 + np.random.rand(500, 3) * 90

    gr._trans_dip_mom_x = trans_dip_mom_x
    gr._trans_dip_mom_y = trans_dip_mom_y
    gr._trans_dip_mom_z = trans_dip_mom_z

    gr.read_data_direct(excit_e, trans_dip_mom_x, trans_dip_mom_y, trans_dip_mom_z)
    return gr

def test_reduction_count(geom_reduction_instance):
    gr = geom_reduction_instance
    gr.reduce_geoms()  
    reduced = gr.subsamples
    expected_count = 20
    assert len(reduced) == expected_count, f"Expected {expected_count} reduced geometries, got {len(reduced)}"

def test_transition_dipole_components_range(geom_reduction_instance):
    gr = geom_reduction_instance
    for comp, label in zip(
        [gr._trans_dip_mom_x, gr._trans_dip_mom_y, gr._trans_dip_mom_z],
        ['x', 'y', 'z']
    ):
        min_val = comp.min()
        max_val = comp.max()
        assert min_val >= 10, f"Transition dipole moment component '{label}' has a minimum value {min_val} which is less than 10."
        assert max_val <= 100, f"Transition dipole moment component '{label}' has a maximum value {max_val} which is greater than 100."

