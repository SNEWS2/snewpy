#!/bin/env python
from pathlib import Path
import numpy as np
import os
import pytest
import tarfile
#mark all tests in this file

snowglobes_dir = os.environ.get('SNOWGLOBES','')
pytestmark = pytest.mark.snowglobes 

basedir = Path.cwd()
outdir = Path(snowglobes_dir)/'out'

#simulation setup
model_file = basedir/'models/Kuroda_2020/LnuR10B12.dat'
model_type = model_file.parent.stem
xform_type = 'AdiabaticMSW_NMO'
detectors= 'all'

fluence_name = f'test_{model_type}_{xform_type}'
fluence_file = model_file.parent/f'{fluence_name}.tar.bz2'
output_archive = model_file.parent/f'{fluence_name}_SNOprocessed.tar.gz'

#checks setup
reference_archive = basedir/'python/snewpy/test/reference_data'/output_archive.name

from snewpy import snowglobes

#-------------------------------------------------
def _cleanup(path, pattern='*.*'):
    "Remove all the files from path/ matching pattern"
    if path.exists():
        for f in path.glob(pattern):
            f.unlink()
    
def _read_data(fname):
    return np.loadtxt(fname,comments=['---','Total','#','Energy'], dtype='f8')

#--------------------------------------------
#main steps defined here
def generate():
    f = snowglobes.generate_fluence(model_file, model_type, xform_type, d=10, 
                                    output_filename=fluence_name)
    return f

def simulate():
    #cleanup output directory first
    _cleanup(outdir)
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

#-------------------------------------------
@pytest.mark.timing
class TestExecutionTime:
    def test_time_to_generate(self,benchmark):
        benchmark(generate)

    def test_time_to_simulate(self,benchmark):
        benchmark(simulate)

    def test_time_to_collate(self,benchmark):
        benchmark(collate)

#-------------------------------------------
class Test_SNOwGLoBES_steps:
    def test_generate(self):
        #cleanup fluence file
        fluence_file.unlink(missing_ok=True)
        #generate a new fluence file
        result = generate()
        assert result == str(fluence_file)
        assert fluence_file.exists()

    def test_simulate(self):
        simulate()
        #get files in output_dir
        files = list(outdir.glob(f'{fluence_name}*.dat'))
        #make sure we got many files
        assert len(files)>100
        for f in files:
            #and they're big enough
            assert f.stat().size_t>100

    def test_collate(self):
        tables = collate()
        assert tables
        files = list(outdir.glob(f'Collate_{fluence_name}*.dat'))
        #make sure we got many files
        assert len(files)>100
        for f in files:
            #and they're big enough
            assert f.stat().st_size>1024 #1KB
#---------------------------------------------
@pytest.mark.crosscheck
def test_run_snoglobes():
    fluence_file.unlink(missing_ok=True)
    generate()
    simulate()
    collate()
    assert output_archive.exists()
    assert output_archive.stat().st_size>102400 #100KB


#@pytest.mark.skipif(Path(reference_archive).exists() == False, reason=f'No reference file {reference_archive} for check')
@pytest.mark.crosscheck
def test_output_files_matches_reference(tmpdir):
    #check that each reference file has the same value
    new_dir = tmpdir/'new'
    ref_dir = tmpdir/'ref'
    #first extract all files
    with tarfile.open(output_archive) as t:
        t.extractall(new_dir)
    with tarfile.open(reference_archive) as t:
        t.extractall(ref_dir)
    #archive files are in a directory
    ref_dir = ref_dir/reference_archive.stem
    new_dir = new_dir/reference_archive.stem
    #now check that each file in the ref_dir is contained in the new_dir
    ref_files = list(Path(ref_dir).glob(f'*.dat'))
    for ref_file in ref_files:
        if ref_file.stat().st_size>0:
            check_file=Path(check_dir)/ref_file.relative_to(ref_dir)
            assert check_file!=file
            assert check_file.exists()
            d0 = read_data(file)
            if d0.size>0:
                d1 = read_data(check_file)
                print(d0-d1)
                assert np.allclose(d0,d1)



