import unittest
import pdb

from evac.utils.grid import Grid

class TestGrid(unittest.TestCase):
    def setUp(self):
        pass

class TestNewGrid(TestGrid):
    def setUp(self):
        opts = dict(urcrnrlat=41.0,
                    urcrnrlon=-4.0,
                    llcrnrlat=40.0,
                    llcrnrlon=-5.0,
                    dx_km=5)
        self.G = Grid(opts)

    def test_newgrid(self):
        pdb.set_trace()

if __name__ == '__main__':
    unittest.main()
