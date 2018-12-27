import unittest
import pdb
import nose
from unittest.mock import Mock

from evac.utils.grid import Grid

class TestNewGrid(unittest.TestCase):
    def setUp(self):
        opts = dict(urcrnrlat=41.0,
                    urcrnrlon=-4.0,
                    llcrnrlat=40.0,
                    llcrnrlon=-5.0,
                    dx_km=5)
        self.G = Grid(opts)

    def test_newgrid_latlons(self):
        delta = 0.001
        latmax = self.G.lats.max()
        latmin = self.G.lats.min()
        lonmax = self.G.lons.max()
        lonmin = self.G.lons.min()
        self.assertAlmostEqual(latmax,41.0,delta=delta)
        self.assertAlmostEqual(latmin,40.0,delta=delta)
        self.assertAlmostEqual(lonmax,-4.0,delta=delta)
        self.assertAlmostEqual(lonmin,-5.0,delta=delta)

        latmaxloc = N.where(self.G.lats == latmax)
        latminloc = N.where(self.G.lats == latmin)
        lonmaxloc = N.where(self.G.lons == lonmax)
        lonminloc = N.where(self.G.lons == lonmin)
        # self.assertEqual(latmaxloc,

        latshape = self.G.lats.shape
        self.assertEqual(

"""
class TestGridWRFOut(unittest.TestCase):
    def setUp(self):
        W = Mock()
        # more here
        self.G = Grid(W)

    def test_wrfout_latlons(self):
        pass

class TestGridEnsemble(unittest.TestCase):
    def setUp(self):
        E = Mock()
        # more here
        self.G = Grid(E)

    def test_ensemble_latlons(self):
        for dom in (1,2,3):
            pass


class TestGridStageIV(unittest.TestCase):
    def setUp(self):
        ST4 = Mock()
        # more here
        self.G = Grid(ST4)

    def test_st4_latlons(self):
        pass

class TestGridReprojections(unittest.TestCase):
    def setUp(self):
        pass

    def test_newgrid_reproject(self):
        # Interpolating to new domain.
        pass
"""

if __name__ == '__main__':
    unittest.main()
