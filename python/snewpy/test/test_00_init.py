import snewpy
import unittest

import pytest
pytestmark=pytest.mark.base

class TestInit(unittest.TestCase):
    def test_version_exists(self):
        self.assertTrue(hasattr(snewpy, '__version__'))