import unittest

from .. import gis_tools

class TestGisTools(unittest.TestCase):
    def setUp(self):
        pass

class TestDecomposeWind(TestGisTools):
    def setUp(self):
        pass

    def test_conversions(self):
        """Make sure conversions work.
        """
        c1 = 'ms_kt'
        c2 = 'ms_mph'
        c3 = 'kt_ms'

        for conv in (c1,c2,c3):
            with self.subTest(conv=conv):
                self.assertEqual()
