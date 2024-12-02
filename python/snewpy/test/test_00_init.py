import snewpy
import importlib
import subprocess

def test_version_exists():
    assert hasattr(snewpy, '__version__')

def test_version_is_consistent_with_variable():
    assert importlib.metadata.version('snewpy') == snewpy.__version__

def test_version_is_consistent_with_git_tag():
    git_tag = subprocess.check_output('git describe --abbrev=0'.split())
    git_tag = git_tag.decode('ascii').strip().lstrip('v')
    assert git_tag == snewpy.__version__