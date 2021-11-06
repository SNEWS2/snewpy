#!/bin/env python
from pathlib import Path
import numpy as np
import os
import pytest


basedir = Path.cwd()
snowglobes_dir = os.environ['SNOWGLOBES']
outdir = Path(snowglobes_dir)/'out'

#simulation setup
model_file = basedir/'models/Kuroda_2020/LnuR10B12.dat'
model_type =model_file.parent.stem
xform_type = 'AdiabaticMSW_NMO'

fluence_name = f'test_{model_type}_{xform_type}'
fluence_file = model_file.parent/f'{fluence_name}.tar.bz2'
detectors= 'all'

#checks setup
ref_dir = basedir/'out_ref1' 

from snewpy import snowglobes

#main steps defined here
def generate():
    f = snowglobes.generate_fluence(model_file, model_type, xform_type, d=10, 
                                    output_filename=fluence_name)
    return f

def simulate():
    #cleanup output directory first
    for f in outdir.glob('*.*'):
        f.unlink()
    snowglobes.simulate(snowglobes_dir,fluence_file, detector_input=detectors)

def collate():
    #recreate the file to avoid crashes
    Path(f'{snowglobes_dir}/fluxes/{fluence_name}.dat').touch()
    try:
        result = snowglobes.collate(snowglobes_dir,fluence_file, detector_input=detectors, 
                                            skip_plots=True, remove_generated_files=False)
    finally:
        #return back to original directory
        os.chdir(basedir)
    return result

def read_data(fname):
    return np.loadtxt(fname,comments=['---','Total','#','Energy'], dtype='f8')

def _check_file_in_directory(ref_file, check_dir):
    """ check that *check_dir* has a file 
    with the same name and contents as *ref_file* """
    file=Path(ref_file)
    check_file = Path(check_dir)/file.name
    assert check_file!=file
    assert check_file.exists()
    d0 = read_data(file)
    if d0.size>0:
        d1 = read_data(check_file)
        print(d0-d1)
        assert np.allclose(d0,d1)

#-------------------------------------------
@pytest.mark.timing
def test_benchmark_generate(benchmark):
    benchmark(generate)

@pytest.mark.timing
def test_benchmark_simulate(benchmark):
    benchmark(simulate)

@pytest.mark.timing
def test_benchmark_collate(benchmark):
    benchmark(collate)
#-------------------------------------------
def test_generate():
    #cleanup fluence file
    fluence_file.unlink(missing_ok=True)
    #generate a new fluence file
    result = generate()
    assert result == str(fluence_file)
    assert fluence_file.exists()

def test_collate():
    tables = collate()
    assert tables

def test_output_files_match():
    #check that each reference file has 
    ref_files = list(Path(ref_dir).glob('*.dat'))
    assert len(ref_files)>0
    for ref_file in ref_files:
        if ref_file.stat().st_size>0:
            _check_file_in_directory(ref_file, outdir)

