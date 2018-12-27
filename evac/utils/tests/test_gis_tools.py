import unittest

import evac.utils.gis_tools as gis_tools
from evac.utils.evac_numbers as evac_numbers


class TestGisTools(unittest.TestCase):
    def setUp(self):
        pass

class TestDecomposeWind(TestGisTools):
    def setUp(self):
        self.func = gis_tools.decompose_wind
        self.wspd = 5.0
        self.wdir = 189.0
        # c[0] = {'conv':None,'out_u':0.782,'out_v':4.938}

    def test_int_input(self):
        u,v = self.func(int(self.wspd,int(self.wdir))
        u_is_int = evac_numbers.is_integer(u)
        v_is_int = evac_numbers.is_integer(v)
        self.assertTrue(u_is_int)
        self.assertTrue(v_is_int)

    def test_mixed_input(self):
        u,v = self.func(int(wspd),float(wdir))
        u_is_int = evac_numbers.is_int(u)
        v_is_float = evac_numbers.is_float(v)
        self.assertTrue(u_is_int)
        self.assertTrue(v_is_float)

        u,v = self.func(float(wspd),int(wdir))
        u_is_float = evac_numbers.is_float(u)
        v_is_int = evac_numbers.is_int(v)
        self.assertTrue(u_is_float)
        self.assertTrue(v_is_int)

    def test_float_input(self):
        u,v = self.func(float(wspd),float(wdir))
        u_is_float = evac_numbers.is_float(u)
        v_is_float = evac_numbers.is_float(v)
        self.assertTrue(u_is_float)
        self.assertTrue(v_is_float)

class TestConvertSpeed(TestGisTools):
    def setUp(self):
        self.func = gis_tools.convert_velocity_units
        self.V = 0.782

    def test_conversions(self):
        """Make sure conversions work, to two decimal places.
        """
        c = {}

        c[0] = {'conv':'ms_kt','out_V':1.520}
        c[1] = {'conv':'ms_mph','out_V':1.749}
        c[2] = {'conv':'kt_mph','out_V':0.900}
        c[3] = {'conv':'kt_ms','out_V':0.402}


        for idx in range(len(c.keys()):
            with self.subTest(conv=c[idx]['conv']):
                u,v = self.func(self.wspd,self.wdir,convert=conv)
                self.assertAlmostEqual(u,c[idx]['out_u'])
                self.assertAlmostEqual(v,c[idx]['out_v'])

    
