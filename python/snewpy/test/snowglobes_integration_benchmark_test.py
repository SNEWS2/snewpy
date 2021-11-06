#!/bin/env python
from pathlib import Path
import numpy as np
import os
import pytest
#mark all tests in this file

snowglobes_dir = os.environ.get('SNOWGLOBES','')
pytestmark = [pytest.mark.snowglobes, 
              pytest.mark.skipif(not snowglobes_dir, reason='Missing $SNOWGLOBES env')]

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
reference_archive = basedir/f'python/snewpy/test/reference_data/{fluence_name}.tar.gz' 

from snewpy import snowglobes

#-------------------------------------------------
def _cleanup(path, pattern='*.*'):
    "Remove all the files from path/ matching pattern"
    if path.exists():
        for f in path.glob(pattern):
            f.unlink()
    
def _read_data(fname):
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
        files = list(outdir.glob('{fluence_name}*.dat'))
        #make sure we got many files
        assert len(files)>1000
        for f in files:
            #and they're big enough
            assert f.stat().size_t>100

    def test_collate(self):
        tables = collate()
        assert tables
        files = list(outdir.glob('Collate_{fluence_name}*.dat'))
        #make sure we got many files
        assert len(files)>1000
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


#Prepare output directories for comparing archive files
@pytest.fixture
def temp_dirs(tmpdir):
    dirs = [tmpdir/'ref', tmpdir/'new']
    for d in dirs:
        d.mkdir()
    return dirs

@pytest.mark.crosscheck
@pytest.mark.skipif(Path(reference_archive).exists() == False, reason='No reference file for check')
def test_output_files_matches_reference(temp_dirs):
    #check that each reference file has the same value
    new_dir,ref_dir = temp_dirsz
    #first extract all files
    with tarfile.open(output_archive) as t:
        t.extractall(new_dir)
    with tarfile.open(reference_archive) as t:
        t.extractall(ref_dir)
    #now check that each file in the ref_dir is contained in the new_dir
    ref_files = list(ref_dir.glob('*.dat'))
    assert len(ref_files)>0
    for ref_file in ref_files:
        if ref_file.stat().st_size>0:
            _check_file_in_directory(ref_file, new_dir)


