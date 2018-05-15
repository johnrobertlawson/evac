import unittest
import pdb

from mock import Mock

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
        expect_lats = 0
        expect_lons = 0


class TestGridWRFOut(unittest.TestCase):
    def setUp(self):
        W = Mock()
        # more here
        self.G = Grid(inst=W)

    def test_wrfout_latlons(self):
        pass

class TestGridEnsemble(unittest.TestCase):
    def setUp(self):
        E = Mock()
        # more here
        self.G = Grid(inst=E)

    def test_ensemble_latlons(self):
        for dom in (1,2,3):
            pass


class TestGridStageIV(unittest.TestCase):
    def setUp(self):
        ST4 = Mock()
        # more here
        self.G = Grid(inst=ST4)

    def test_st4_latlons(self):
        pass

class TestGridReprojections(self):
    def setUp(self):
        pass

    def test_newgrid_reproject(self):
        # Interpolating to new domain.

if __name__ == '__main__':
    unittest.main()
