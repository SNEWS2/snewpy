import snewpy
import unittest


class TestInit(unittest.TestCase):
    def test_version_exists(self):
        self.assertTrue(hasattr(snewpy, '__version__'))