import snewpy
  
def test_version_exists():
    assert hasattr(snewpy, '__version__')
