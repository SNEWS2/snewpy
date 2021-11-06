#!/bin/env python
from pathlib import Path
import numpy as np

import os
import snewpy.snowglobes

import pytest

#def pytest_addoption(parser):
#    parser.addoption("--snmodel", action="store_true", help="run all combinations")

detectors= 'all'
basedir = Path.cwd()
model_file = basedir/'models/OConnor_2015/M1_neutrinos.dat'
model_file = basedir/'models/Kuroda_2020/LnuR10B12.dat'
ref_dir =    basedir/'out_ref1'
snowglobes_dir = os.environ['SNOWGLOBES']
outdir = Path(snowglobes_dir)/'out'

fluence_file = model_file.parent/'fluence.tar.bz2'


def generate():
    return snewpy.snowglobes.generate_fluence(
                                    model_path=model_file, 
                                    model_type=model_file.parent.stem,
                                    transformation_type='NoTransformation', 
                                    d=5, 
                                    output_filename='fluence'
                                    )

def simulate():
    snewpy.snowglobes.simulate(snowglobes_dir,fluence_file, detector_input=detectors)

def collate():
    #recreate the file to avoid crashes
    Path(snowglobes_dir+'/fluxes/fluence.dat').touch()
    snewpy.snowglobes.collate(snowglobes_dir,fluence_file, detector_input=detectors, 
                                        skip_plots=True, remove_generated_files=False)
    #return back to original directory
    os.chdir(basedir)

def read_data(fname):
    return np.loadtxt(fname,comments=['---','Total','#','Energy'], dtype='f8')

def check_file_in_directory(ref_file, check_dir):
    """ check that *check_dir* has a file 
    with the same name and contents as *ref_file* """
    file=Path(ref_file)
    check_file = Path(check_dir)/file.name
    assert check_file.exists()
    d0 = read_data(file)
    if d0.size>0:
        d1 = read_data(check_file)
        print(d0-d1)
        assert np.allclose(d0,d1)

#-------------------------------------------
def notest_generate():
    #cleanup fluence file
    fluence_file.unlink(missing_ok=True)
    #generate a new fluence file
    result = generate()
    assert result == str(fluence_file)
    assert fluence_file.exists()

#def test_benchmark_generate(benchmark):
#    benchmark(generate)

def test_benchmark_simulate(benchmark):
    benchmark(simulate)

def test_benchmark_collate(benchmark):
    benchmark(collate)

def test_output_files_match():
    #check that each reference file has 
    ref_files = list(Path(ref_dir).glob('*.dat'))
    assert len(ref_files)>0
    for ref_file in ref_files:
        if ref_file.stat().st_size>0:
            check_file_in_directory(ref_file, outdir)


if __name__ == '__main__':
    simulate()
    collate()
