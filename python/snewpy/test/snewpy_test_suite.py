import unittest

def snewpy_test_suite():
    """Returns unittest.TestSuite of desiutil tests.

    This is factored out separately from runtests() so that it can be used by
    ``python setup.py test``.
    """
    from os.path import dirname
    pydir = dirname(dirname(__file__))
    tests = unittest.defaultTestLoader.discover(pydir,
                                                top_level_dir=dirname(pydir))
    return tests

def runtests():
    """Run all tests in snewpy.test.test_*.
    """
    # Load all TestCase classes from snewpy/test/test_*.py
    tests = snewpy_test_suite()
    # Run them
    unittest.TextTestRunner(verbosity=2).run(tests)

if __name__ == "__main__":
    runtests()
