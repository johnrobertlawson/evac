"""Colour tables for plotting.

Adapted from code by Luke Madaus and David-John Gagne II

Todo:
    * Probably should use American English for the title...
"""
import pdb

from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt

# def default(*args):
    # jet_cdict = plt.cm.jet
    # jet_tbl = LinearSegmentedColormap('jet_TBL',jet_cdict)
    # pdb.set_trace()
    # return jet_tbl

def RdBufloat(*args):
    valrange = args
    # Will define a colortable that keeps zero as white
    length = max(valrange)-min(valrange)
    distance = 0 - min(valrange)
    pct = float(distance)/float(length)
    
    RdBufloat_cdict = {'red':((0.00, 0.00, 0.00),
                 (pct, 1.00, 1.00),
                (1.00, 1.00, 1.00)),
        'green':    ((0.00, 0.00, 0.00),
                 (pct, 1.00, 1.00),
                (1.00, 0.00, 0.00)),
        'blue':        ((0.00, 1.00, 1.00),
                 (pct, 1.00, 1.00),
                (1.00, 0.00, 0.00))}
    RdBufloat_coltbl = LinearSegmentedColormap('RdBufloat_COLTBL',RdBufloat_cdict)
    return RdBufloat_coltbl



def rain1(*args):
    rain_cdict ={'red':    ((0.00, 0.00, 0.00),
                (1.00, 0.00, 0.00)),
        'green':    ((0.00, 1.00, 1.00),
                (1.00, 0.20, 0.20)),
        'blue':        ((0.00, 0.00, 0.00),
                (1.00, 0.00, 0.00))}
    rain_coltbl = LinearSegmentedColormap('RAIN_COLTBL',rain_cdict)
    return rain_coltbl

def snow1(*args):
    snow_cdict ={'red':    ((0.00, 0.00, 0.00),
                (1.00, 0.00, 0.00)),
        'green':    ((0.00, 0.00, 0.00),
                (1.00, 0.00, 0.00)),
        'blue':        ((0.00, 1.00, 1.00),
                (1.00, 0.20, 0.20))}
    snow_coltbl = LinearSegmentedColormap('SNOW_COLTBL',snow_cdict)
    return snow_coltbl

def mixprecip1(*args):
    mix_cdict ={'red':    ((0.00, 0.20, 0.20),
                (1.00, 1.00, 1.00)),
        'green':    ((0.00, 0.00, 0.00),
                (1.00, 0.00, 0.00)),
        'blue':        ((0.00, 0.00, 0.00),
                (1.00, 0.00, 0.00))}
    mix_coltbl = LinearSegmentedColormap('MIX_COLTBL',mix_cdict)
    return mix_coltbl

def grays(*args):
    grays_cdict ={'red':    ((0.00, 1.00, 1.00),
                (1.00, 0.05, 0.05)),
        'green':    ((0.00, 1.00, 1.00),
                (1.00, 0.05, 0.05)),
        'blue':        ((0.00, 1.00, 1.00),
                (1.00, 0.05, 0.05))}
    grays_coltbl = LinearSegmentedColormap('GRAYS_COLTBL',grays_cdict)
    return grays_coltbl

def snow2(*args):
    snowf_cdict ={'red':    ((0.00, 0.91, 0.91),
                (0.06, 0.81, 0.81),
                (0.12, 0.51, 0.51),
                (0.18, 0.23, 0.23),
                (0.24, 0.11, 0.11),
                (0.30, 0.00, 0.00),
                (0.36, 0.02, 0.02),
                (0.42, 0.02, 0.02),
                (0.48, 0.03, 0.03),
                (0.54, 0.52, 0.52),
                (0.60, 1.00, 1.00),
                (0.66, 1.00, 1.00),
                (0.72, 1.00, 1.00),
                (0.78, 1.00, 1.00),
                (0.84, 0.70, 0.70),
                (0.90, 0.40, 0.40),
                (1.00, 0.20, 0.20)),
        'green':    ((0.00, 0.80, 0.80),
                (0.06, 0.50, 0.50),
                (0.12, 0.20, 0.20),
                (0.18, 0.00, 0.00),
                (0.24, 0.00, 0.00),
                (0.30, 0.00, 0.00),
                (0.36, 0.24, 0.24),
                (0.42, 0.47, 0.47),
                (0.48, 0.70, 0.70),
                (0.54, 0.85, 0.85),
                (0.60, 1.00, 1.00),
                (0.66, 0.67, 0.67),
                (0.72, 0.33, 0.33),
                (0.78, 0.00, 0.00),
                (0.84, 0.00, 0.00),    
                (0.90, 0.00, 0.00),
                (1.00, 0.00, 0.00)),
        'blue':        ((0.00, 0.98, 0.98),
                (0.06, 0.87, 0.87),
                (0.12, 0.58, 0.58),
                (0.18, 0.69, 0.69),
                (0.24, 0.84, 0.84),
                (0.30, 1.00, 1.00),
                (0.36, 0.69, 0.69),
                (0.42, 0.37, 0.37),
                (0.48, 0.06, 0.06),
                (0.54, 0.03, 0.03),
                (0.60, 0.00, 0.00),
                (0.66, 0.00, 0.00),
                (0.72, 0.00, 0.00),
                (0.78, 0.00, 0.00),
                (0.84, 0.00, 0.00),
                (0.90, 0.00, 0.00),
                (1.00, 0.00, 0.00))}
    snowf_coltbl = LinearSegmentedColormap('SNOWF_COLTBL',snowf_cdict)
    return snowf_coltbl

def reflect(*args):
    reflect_cdict ={'red':    ((0.000, 0.40, 0.40),
                (0.067, 0.20, 0.20),
                (0.133, 0.00, 0.00),
                (0.200, 0.00, 0.00),
                (0.267, 0.00, 0.00),
                (0.333, 0.00, 0.00),
                (0.400, 1.00, 1.00),
                (0.467, 1.00, 1.00),
                (0.533, 1.00, 1.00),
                (0.600, 1.00, 1.00),
                (0.667, 0.80, 0.80),
                (0.733, 0.60, 0.60),
                (0.800, 1.00, 1.00),
                (0.867, 0.60, 0.60),
                (0.933, 1.00, 1.00),
                (1.000, 0.00, 0.00)),
        'green':    ((0.000, 1.00, 1.00),
                (0.067, 0.60, 0.60),
                (0.133, 0.00, 0.00),
                (0.200, 1.00, 1.00),
                (0.267, 0.80, 0.80),
                (0.333, 0.60, 0.60),
                (0.400, 1.00, 1.00),
                (0.467, 0.80, 0.80),
                (0.533, 0.40, 0.40),
                (0.600, 0.00, 0.00),
                (0.667, 0.20, 0.20),
                (0.733, 0.00, 0.00),
                (0.800, 0.00, 0.00),
                (0.867, 0.20, 0.20),
                (0.933, 1.00, 1.00),
                (1.000, 1.00, 1.00)),
        'blue':        ((0.000, 1.00, 1.00),
                (0.067, 1.00, 1.00),
                (0.133, 1.00, 1.00),
                (0.200, 0.00, 0.00),
                (0.267, 0.00, 0.00),
                (0.333, 0.00, 0.00),
                (0.400, 0.00, 0.00),
                (0.467, 0.00, 0.00),
                (0.533, 0.00, 0.00),
                (0.600, 0.00, 0.00),
                (0.667, 0.00, 0.00),
                (0.733, 0.00, 0.00),
                (0.800, 1.00, 1.00),
                (0.867, 0.80, 0.80),
                (0.933, 1.00, 1.00),
                (1.000, 1.00, 1.00))}
    reflect_coltbl = LinearSegmentedColormap('REFLECT_COLTBL',reflect_cdict)
    return reflect_coltbl

def universal_reflect():
    """
    Map a number to a colour.
    """
    pass

def ncdc_modified_ISU(nlvs):
    cmap_dict = {'red':((0.0000, 0.000, 0.000),
                (0.0000, 0.000, 0.000),
                (0.0000, 0.000, 0.000),
                (0.0000, 0.000, 0.000),
                (0.0714, 0.000, 0.000),
                (0.1429, 0.000, 0.000),
                (0.2143, 0.000, 0.000),
                (0.2857, 0.000, 0.000),
                (0.3571, 0.000, 0.000),
                (0.4286, 1.000, 1.000),
                (0.5000, 0.906, 0.906),
                (0.5714, 1.000, 1.000),
                (0.6429, 1.000, 1.000),
                (0.7143, 0.839, 0.839),
                (0.7857, 0.753, 0.753),
                (0.8571, 1.000, 1.000)),
        'green':    ((0.0000, 0.925, 0.925),
                (0.0000, 0.925, 0.925),
                (0.0000, 0.925, 0.925),
                (0.0000, 0.925, 0.925),
                (0.0714, 0.627, 0.627),
                (0.1429, 0.000, 0.000),
                (0.2143, 1.000, 1.000),
                (0.2857, 0.784, 0.784),
                (0.3571, 0.565, 0.565),
                (0.4286, 1.000, 1.000),
                (0.5000, 0.753, 0.753),
                (0.5714, 0.565, 0.565),
                (0.6429, 0.000, 0.000),
                (0.7143, 0.000, 0.000),
                (0.7857, 0.000, 0.000),
                (0.8571, 0.000, 0.000)),

        'blue':        ((0.0000, 0.925, 0.925),
                (0.0000, 0.925, 0.925),
                (0.0000, 0.925, 0.925),
                (0.0000, 0.925, 0.925),
                (0.0714, 0.965, 0.965),
                (0.1429, 0.965, 0.965),
                (0.2143, 0.000, 0.000),
                (0.2857, 0.000, 0.000),
                (0.3571, 0.000, 0.000),
                (0.4286, 0.000, 0.000),
                (0.5000, 0.000, 0.000),
                (0.5714, 0.000, 0.000),
                (0.6429, 0.000, 0.000),
                (0.7143, 0.000, 0.000),
                (0.7857, 0.000, 0.000),
                (0.8571, 1.000, 1.000))}
    mod_coltbl = LinearSegmentedColormap('mod_cref',cmap_dict,nlvs)
    return mod_coltbl

def reflect_ncdc(*args):



    reflect_ncdc_cdict ={'red':((0.0000, 0.000, 0.000),
                (0.0714, 0.000, 0.000),
                (0.1429, 0.000, 0.000),
                (0.2143, 0.000, 0.000),
                (0.2857, 0.000, 0.000),
                (0.3571, 0.000, 0.000),
                (0.4286, 1.000, 1.000),
                (0.5000, 0.906, 0.906),
                (0.5714, 1.000, 1.000),
                (0.6429, 1.000, 1.000),
                (0.7143, 0.839, 0.839),
                (0.7857, 0.753, 0.753),
                (0.8571, 1.000, 1.000),
                (0.9286, 0.600, 0.600),
                (1.000, 0.923, 0.923)),
        'green':    ((0.0000, 0.925, 0.925),
                (0.0714, 0.627, 0.627),
                (0.1429, 0.000, 0.000),
                (0.2143, 1.000, 1.000),
                (0.2857, 0.784, 0.784),
                (0.3571, 0.565, 0.565),
                (0.4286, 1.000, 1.000),
                (0.5000, 0.753, 0.753),
                (0.5714, 0.565, 0.565),
                (0.6429, 0.000, 0.000),
                (0.7143, 0.000, 0.000),
                (0.7857, 0.000, 0.000),
                (0.8571, 0.000, 0.000),
                (0.9286, 0.333, 0.333),
                (1.000, 0.923, 0.923)),

        'blue':        ((0.0000, 0.925, 0.925),
                (0.0714, 0.965, 0.965),
                (0.1429, 0.965, 0.965),
                (0.2143, 0.000, 0.000),
                (0.2857, 0.000, 0.000),
                (0.3571, 0.000, 0.000),
                (0.4286, 0.000, 0.000),
                (0.5000, 0.000, 0.000),
                (0.5714, 0.000, 0.000),
                (0.6429, 0.000, 0.000),
                (0.7143, 0.000, 0.000),
                (0.7857, 0.000, 0.000),
                (0.8571, 1.000, 1.000),
                (0.9286, 0.788, 0.788),
                (1.000, 0.923, 0.923))}
    reflect_ncdc_coltbl = LinearSegmentedColormap('REFLECT_NCDC_COLTBL',reflect_ncdc_cdict)
    return reflect_ncdc_coltbl

def irsat(*args):
    irsat_cdict ={'red':    ((0.000, 1.000, 0.294),
                (0.067, 1.000, 1.000),
                (0.133, 0.804, 0.804),
                (0.200, 0.369, 0.369),
                (0.267, 0.627, 0.627),
                (0.333, 0.804, 0.804),
                (0.400, 1.000, 1.000),
                (0.567, 0.000, 0.000),
                (0.667, 0.400, 0.400),
                (0.700, 0.596, 0.596),
                (0.800, 0.000, 0.000),
                (0.867, 0.416, 0.416),
                (0.933, 0.804, 0.804),
                (1.000, 0.294, 0.294)),
        'green':    ((0.000, 1.000, 0.000),
                (0.067, 0.000, 0.000),
                (0.133, 0.361, 0.361),
                (0.200, 0.149, 0.149),
                (0.267, 0.322, 0.322),
                (0.333, 0.584, 0.584),
                (0.400, 0.757, 0.757),
                (0.567, 0.392, 0.392),
                (0.667, 0.804, 0.804),
                (0.700, 0.961, 0.961),
                (0.800, 0.000, 0.000),
                (0.867, 0.353, 0.353),
                (0.933, 0.000, 0.000),
                (1.000, 0.000, 0.000)),
        'blue':        ((0.000, 1.000, 1.000),
                (0.067, 0.000, 0.000),
                (0.133, 0.360, 0.360),
                (0.200, 0.070, 0.070),
                (0.267, 0.176, 0.176),
                (0.333, 0.047, 0.047),
                (0.400, 0.145, 0.145),
                (0.567, 0.000, 0.000),
                (0.667, 0.667, 0.667),
                (0.700, 1.000, 1.000),
                (0.800, 0.502, 0.502),
                (0.867, 0.804, 0.804),
                (0.933, 0.804, 0.804),
                (1.000, 0.510, 0.510))}
    irsat_coltbl = LinearSegmentedColormap('IRSAT_COLTBL',irsat_cdict)
    return irsat_coltbl

def bw_irsat(*args):
    bw_irsat_cdict ={'red':    ((0.000, 1.000, 1.000),
                (1.000, 0.000, 0.000)),
        'green':    ((0.000, 1.000, 1.000),
                (1.000, 0.000, 0.000)),
        'blue':        ((0.000, 1.000, 1.000),
                (1.000, 0.000, 0.000))}
    bw_irsat_coltbl = LinearSegmentedColormap('BW_IRSAT_COLTBL',bw_irsat_cdict)
    return bw_irsat_coltbl


def precip1(*args):
    precip_cdict ={'red':    ((0.000, 1.000, 1.000),
                (0.004, 0.914, 0.914),
                (0.012, 0.812, 0.812),
                (0.020, 0.514, 0.514),
                (0.040, 0.227, 0.227),
                (0.060, 0.114, 0.114),
                (0.080, 0.000, 0.000),
                (0.100, 0.012, 0.012),
                (0.120, 0.020, 0.020),
                (0.160, 0.031, 0.031),
                (0.200, 0.518, 0.518),
                (0.240, 1.000, 1.000),
                (0.280, 1.000, 1.000),
                (0.320, 1.000, 1.000),
                (0.360, 1.000, 1.000),
                (0.400, 0.702, 0.702),
                (0.500, 0.490, 0.490),
                (0.600, 0.294, 0.294),
                (0.700, 0.196, 0.196),
                (0.800, 0.980, 0.980),
                (1.000, 1.000, 1.000)),
        'green':    ((0.000, 1.000, 0.000),
                (0.004, 0.800, 0.800),
                (0.012, 0.502, 0.502),
                (0.020, 0.200, 0.200),
                (0.040, 0.000, 0.000),
                (0.060, 0.000, 0.000),
                (0.080, 0.000, 0.000),
                (0.100, 0.235, 0.235),
                (0.120, 0.467, 0.467),
                (0.160, 0.702, 0.702),
                (0.200, 0.851, 0.851),
                (0.240, 1.000, 1.000),
                (0.280, 0.667, 0.667),
                (0.320, 0.227, 0.227),
                (0.360, 0.000, 0.000),
                (0.400, 0.000, 0.000),
                (0.500, 0.000, 0.000),
                (0.600, 0.000, 0.000),
                (0.700, 0.000, 0.000),
                (0.800, 0.773, 0.773),
                (1.000, 1.000, 1.000)),
        'blue':        ((0.000, 1.000, 1.000),
                (0.004, 0.976, 0.976),
                (0.012, 0.875, 0.875),
                (0.020, 0.576, 0.576),
                (0.040, 0.690, 0.690),
                (0.060, 0.843, 0.843),
                (0.080, 1.000, 1.000),
                (0.100, 0.686, 0.686),
                (0.120, 0.372, 0.372),
                (0.160, 0.059, 0.059),
                (0.200, 0.031, 0.031),
                (0.240, 0.000, 0.000),
                (0.280, 0.000, 0.000),
                (0.320, 0.000, 0.000),
                (0.360, 0.000, 0.000),
                (0.400, 0.000, 0.000),
                (0.500, 0.000, 0.000),
                (0.600, 0.000, 0.000),
                (0.700, 0.000, 0.000),
                (0.800, 0.980, 0.980),
                (1.000, 1.000, 1.000))}
    precip_coltbl = LinearSegmentedColormap('PRECIP_COLTBL',precip_cdict)
    return precip_coltbl

def dewpoint1(*args):
    dwp_cdict ={'red':    ((0.00, 0.60, 0.60),
                (0.35, 0.70, 0.70),
                (0.40, 0.80, 0.80),
                (0.45, 0.90, 0.90),
                (0.50, 1.00, 1.00),
                (0.55, 0.90, 0.90),
                (0.60, 0.76, 0.76),
                (0.70, 0.64, 0.64),
                (0.75, 0.52, 0.52),
                (0.85, 0.42, 0.42),
                (1.00, 0.32, 0.32)),
        'green':    ((0.00, 0.33, 0.33),
                (0.35, 0.44, 0.44),
                (0.40, 0.56, 0.56),
                (0.45, 0.69, 0.69),
                (0.50, 0.85, 0.85),
                (0.55, 1.00, 1.00),
                (0.60, 0.90, 0.90),
                (0.70, 0.80, 0.80),
                (0.75, 0.70, 0.70),
                (0.85, 0.60, 0.60),
                (1.00, 0.50, 0.50)),
        'blue':        ((0.00, 0.06, 0.06),
                (0.35, 0.17, 0.17),
                (0.40, 0.32, 0.32),
                (0.45, 0.49, 0.49),
                (0.50, 0.70, 0.70),
                (0.55, 0.70, 0.70),
                (0.60, 0.49, 0.49),
                (0.70, 0.32, 0.32),
                (0.75, 0.17, 0.17),
                (0.85, 0.06, 0.06),
                (1.00, 0.05, 0.05))}
     
    dwp_coltbl = LinearSegmentedColormap('DWP_COLTBL',dwp_cdict)
    return dwp_coltbl         

def sftemp(*args):
    sfc_cdict ={'red':    ((0.00, 0.20, 0.20),
                (0.08, 0.40, 0.40),
                (0.17, 0.27, 0.27),
                (0.25, 0.80, 0.80),
                (0.33, 0.20, 0.20),
                (0.42, 0.20, 0.20),
                (0.50, 0.00, 0.00),
                (0.58, 0.99, 0.99),
                (0.67, 1.00, 1.00),
                (0.75, 0.82, 0.82),
                (0.83, 0.53, 0.53),
                (0.92, 0.95, 0.95),
                (1.00, 1.00, 1.00)),

        'green':    ((0.00, 0.20, 0.20),
                (0.08, 0.40, 0.40),
                (0.17, 0.00, 0.00),
                (0.25, 0.60, 0.60),
                (0.33, 0.40, 0.40),
                (0.42, 0.60, 0.60),
                (0.50, 0.39, 0.39),
                (0.58, 0.76, 0.76),
                (0.67, 0.36, 0.36),
                (0.75, 0.02, 0.02),
                (0.83, 0.00, 0.00),
                (0.92, 0.03, 0.03),
                (1.00, 0.60, 0.60)),

        'blue':        ((0.00, 0.60, 0.60),
                (0.08, 0.60, 0.60),
                (0.17, 0.65, 0.65),
                (0.25, 1.00, 1.00),
                (0.33, 1.00, 1.00),
                (0.42, 0.40, 0.40),
                (0.50, 0.07, 0.07),
                (0.58, 0.02, 0.02),
                (0.67, 0.00, 0.00),
                (0.75, 0.01, 0.01),
                (0.83, 0.00, 0.00),
                (0.92, 0.52, 0.52),
                (1.00, 0.80, 0.80))}

     
    sfc_coltbl = LinearSegmentedColormap('SFC_COLTBL',sfc_cdict)
    return sfc_coltbl

def thetae(*args):
    thte_cdict ={'red':    ((0.00, 0.20, 0.20),
                (0.08, 0.40, 0.40),
                (0.17, 0.27, 0.27),
                (0.25, 0.80, 0.80),
                (0.33, 0.20, 0.20),
                (0.42, 0.20, 0.20),
                (0.50, 0.00, 0.00),
                (0.58, 0.99, 0.99),
                (0.67, 1.00, 1.00),
                (0.75, 0.82, 0.82),
                (0.83, 0.53, 0.53),
                (0.92, 0.95, 0.95),
                (1.00, 1.00, 1.00)),

        'green':    ((0.00, 0.20, 0.20),
                (0.08, 0.40, 0.40),
                (0.17, 0.00, 0.00),
                (0.25, 0.60, 0.60),
                (0.33, 0.40, 0.40),
                (0.42, 0.60, 0.60),
                (0.50, 0.39, 0.39),
                (0.58, 0.76, 0.76),
                (0.67, 0.36, 0.36),
                (0.75, 0.02, 0.02),
                (0.83, 0.00, 0.00),
                (0.92, 0.03, 0.03),
                (1.00, 0.60, 0.60)),

        'blue':        ((0.00, 0.60, 0.60),
                (0.08, 0.60, 0.60),
                (0.17, 0.65, 0.65),
                (0.25, 1.00, 1.00),
                (0.33, 1.00, 1.00),
                (0.42, 0.40, 0.40),
                (0.50, 0.07, 0.07),
                (0.58, 0.02, 0.02),
                (0.67, 0.00, 0.00),
                (0.75, 0.01, 0.01),
                (0.83, 0.00, 0.00),
                (0.92, 0.52, 0.52),
                (1.00, 0.80, 0.80))}

    thte_coltbl = LinearSegmentedColormap('THTE_COLTBL',thte_cdict)
    return thte_coltbl

def PkBlfloat(datarange):
   distance = max(datarange)-min(datarange)
   zeroloc = 0-min(datarange)
   pct = float(zeroloc)/float(distance)

   # Now, rescale the colormap so each value
   # is scaled to equal that

   PkBl_data = {'red':   [(0.*pct, 0.1178, 0.1178),(0.015873*pct, 0.195857, 0.195857),
                        (0.031746*pct, 0.250661, 0.250661),(0.047619*pct, 0.295468, 0.295468),
                        (0.063492*pct, 0.334324, 0.334324),(0.079365*pct, 0.369112, 0.369112),
                        (0.095238*pct, 0.400892, 0.400892),(0.111111*pct, 0.430331, 0.430331),
                        (0.126984*pct, 0.457882, 0.457882),(0.142857*pct, 0.483867, 0.483867),
                        (0.158730*pct, 0.508525, 0.508525),(0.174603*pct, 0.532042, 0.532042),
                        (0.190476*pct, 0.554563, 0.554563),(0.206349*pct, 0.576204, 0.576204),
                        (0.222222*pct, 0.597061, 0.597061),(0.238095*pct, 0.617213, 0.617213),
                        (0.253968*pct, 0.636729, 0.636729),(0.269841*pct, 0.655663, 0.655663),
                        (0.285714*pct, 0.674066, 0.674066),(0.301587*pct, 0.691980, 0.691980),
                        (0.317460*pct, 0.709441, 0.709441),(0.333333*pct, 0.726483, 0.726483),
                        (0.349206*pct, 0.743134, 0.743134),(0.365079*pct, 0.759421, 0.759421),
                        (0.380952*pct, 0.766356, 0.766356),(0.396825*pct, 0.773229, 0.773229),
                        (0.412698*pct, 0.780042, 0.780042),(0.428571*pct, 0.786796, 0.786796),
                        (0.444444*pct, 0.793492, 0.793492),(0.460317*pct, 0.800132, 0.800132),
                        (0.476190*pct, 0.806718, 0.806718),(0.492063*pct, 0.813250, 0.813250),
                        (0.507937*pct, 0.819730, 0.819730),(0.523810*pct, 0.826160, 0.826160),
                        (0.539683*pct, 0.832539, 0.832539),(0.555556*pct, 0.838870, 0.838870),
                        (0.571429*pct, 0.845154, 0.845154),(0.587302*pct, 0.851392, 0.851392),
                        (0.603175*pct, 0.857584, 0.857584),(0.619048*pct, 0.863731, 0.863731),
                        (0.634921*pct, 0.869835, 0.869835),(0.650794*pct, 0.875897, 0.875897),
                        (0.666667*pct, 0.881917, 0.881917),(0.682540*pct, 0.887896, 0.887896),
                        (0.698413*pct, 0.893835, 0.893835),(0.714286*pct, 0.899735, 0.899735),
                        (0.730159*pct, 0.905597, 0.905597),(0.746032*pct, 0.911421, 0.911421),
                        (0.761905*pct, 0.917208, 0.917208),(0.777778*pct, 0.922958, 0.922958),
                        (0.793651*pct, 0.928673, 0.928673),(0.809524*pct, 0.934353, 0.934353),
                        (0.825397*pct, 0.939999, 0.939999),(0.841270*pct, 0.945611, 0.945611),
                        (0.857143*pct, 0.951190, 0.951190),(0.873016*pct, 0.956736, 0.956736),
                        (0.888889*pct, 0.962250, 0.962250),(0.904762*pct, 0.967733, 0.967733),
                        (0.920635*pct, 0.973185, 0.973185),(0.936508*pct, 0.978607, 0.978607),
                        (0.952381*pct, 0.983999, 0.983999),(0.968254*pct, 0.989361, 0.989361),
                        (0.984127*pct, 0.994695, 0.994695),(1.0*pct, 1.0, 1.0),
            ((1-pct)*0.125+pct,
                0.87058824300765991, 0.87058824300765991), ((1-pct)*0.25+pct,
                0.7764706015586853, 0.7764706015586853), ((1-pct)*0.375+pct,
                0.61960786581039429, 0.61960786581039429), ((1-pct)*0.5+pct,
                0.41960784792900085, 0.41960784792900085), ((1-pct)*0.625+pct,
                0.25882354378700256, 0.25882354378700256), ((1-pct)*0.75+pct,
                0.12941177189350128, 0.12941177189350128), ((1-pct)*0.875+pct,
                0.031372550874948502, 0.031372550874948502), (1.0,
                0.031372550874948502, 0.031372550874948502)],




              'green': [(0., 0., 0.),(0.015873*pct, 0.102869, 0.102869),
                        (0.031746*pct, 0.145479, 0.145479),(0.047619*pct, 0.178174, 0.178174),
                        (0.063492*pct, 0.205738, 0.205738),(0.079365*pct, 0.230022, 0.230022),
                        (0.095238*pct, 0.251976, 0.251976),(0.111111*pct, 0.272166, 0.272166),
                        (0.126984*pct, 0.290957, 0.290957),(0.142857*pct, 0.308607, 0.308607),
                        (0.158730*pct, 0.325300, 0.325300),(0.174603*pct, 0.341178, 0.341178),
                        (0.190476*pct, 0.356348, 0.356348),(0.206349*pct, 0.370899, 0.370899),
                        (0.222222*pct, 0.384900, 0.384900),(0.238095*pct, 0.398410, 0.398410),
                        (0.253968*pct, 0.411476, 0.411476),(0.269841*pct, 0.424139, 0.424139),
                        (0.285714*pct, 0.436436, 0.436436),(0.301587*pct, 0.448395, 0.448395),
                        (0.317460*pct, 0.460044, 0.460044),(0.333333*pct, 0.471405, 0.471405),
                        (0.349206*pct, 0.482498, 0.482498),(0.365079*pct, 0.493342, 0.493342),
                        (0.380952*pct, 0.517549, 0.517549),(0.396825*pct, 0.540674, 0.540674),
                        (0.412698*pct, 0.562849, 0.562849),(0.428571*pct, 0.584183, 0.584183),
                        (0.444444*pct, 0.604765, 0.604765),(0.460317*pct, 0.624669, 0.624669),
                        (0.476190*pct, 0.643958, 0.643958),(0.492063*pct, 0.662687, 0.662687),
                        (0.507937*pct, 0.680900, 0.680900),(0.523810*pct, 0.698638, 0.698638),
                        (0.539683*pct, 0.715937, 0.715937),(0.555556*pct, 0.732828, 0.732828),
                        (0.571429*pct, 0.749338, 0.749338),(0.587302*pct, 0.765493, 0.765493),
                        (0.603175*pct, 0.781313, 0.781313),(0.619048*pct, 0.796819, 0.796819),
                        (0.634921*pct, 0.812029, 0.812029),(0.650794*pct, 0.826960, 0.826960),
                        (0.666667*pct, 0.841625, 0.841625),(0.682540*pct, 0.856040, 0.856040),
                        (0.698413*pct, 0.870216, 0.870216),(0.714286*pct, 0.884164, 0.884164),
                        (0.730159*pct, 0.897896, 0.897896),(0.746032*pct, 0.911421, 0.911421),
                        (0.761905*pct, 0.917208, 0.917208),(0.777778*pct, 0.922958, 0.922958),
                        (0.793651*pct, 0.928673, 0.928673),(0.809524*pct, 0.934353, 0.934353),
                        (0.825397*pct, 0.939999, 0.939999),(0.841270*pct, 0.945611, 0.945611),
                        (0.857143*pct, 0.951190, 0.951190),(0.873016*pct, 0.956736, 0.956736),
                        (0.888889*pct, 0.962250, 0.962250),(0.904762*pct, 0.967733, 0.967733),
                        (0.920635*pct, 0.973185, 0.973185),(0.936508*pct, 0.978607, 0.978607),
                        (0.952381*pct, 0.983999, 0.983999),(0.968254*pct, 0.989361, 0.989361),
                        (0.984127*pct, 0.994695, 0.994695),(1.0*pct, 1.0, 1.0), 
            ((1-pct)*0.125+pct,    0.92156863212585449, 0.92156863212585449), ((1-pct)*0.25+pct,
                0.85882353782653809, 0.85882353782653809), ((1-pct)*0.375+pct,
                0.7921568751335144, 0.7921568751335144), ((1-pct)*0.5+pct,
                0.68235296010971069, 0.68235296010971069), ((1-pct)*0.625+pct,
                0.57254904508590698, 0.57254904508590698), ((1-pct)*0.75+pct,
                0.44313725829124451, 0.44313725829124451), ((1-pct)*0.875+pct,
                0.31764706969261169, 0.31764706969261169), ((1-pct)*1.0+pct,
                0.18823529779911041, 0.18823529779911041)],


              'blue':  [(0.*pct, 0., 0.),(0.015873*pct, 0.102869, 0.102869),
                        (0.031746*pct, 0.145479, 0.145479),(0.047619*pct, 0.178174, 0.178174),
                        (0.063492*pct, 0.205738, 0.205738),(0.079365*pct, 0.230022, 0.230022),
                        (0.095238*pct, 0.251976, 0.251976),(0.111111*pct, 0.272166, 0.272166),
                        (0.126984*pct, 0.290957, 0.290957),(0.142857*pct, 0.308607, 0.308607),
                        (0.158730*pct, 0.325300, 0.325300),(0.174603*pct, 0.341178, 0.341178),
                        (0.190476*pct, 0.356348, 0.356348),(0.206349*pct, 0.370899, 0.370899),
                        (0.222222*pct, 0.384900, 0.384900),(0.238095*pct, 0.398410, 0.398410),
                        (0.253968*pct, 0.411476, 0.411476),(0.269841*pct, 0.424139, 0.424139),
                        (0.285714*pct, 0.436436, 0.436436),(0.301587*pct, 0.448395, 0.448395),
                        (0.317460*pct, 0.460044, 0.460044),(0.333333*pct, 0.471405, 0.471405),
                        (0.349206*pct, 0.482498, 0.482498),(0.365079*pct, 0.493342, 0.493342),
                        (0.380952*pct, 0.503953, 0.503953),(0.396825*pct, 0.514344, 0.514344),
                        (0.412698*pct, 0.524531, 0.524531),(0.428571*pct, 0.534522, 0.534522),
                        (0.444444*pct, 0.544331, 0.544331),(0.460317*pct, 0.553966, 0.553966),
                        (0.476190*pct, 0.563436, 0.563436),(0.492063*pct, 0.572750, 0.572750),
                        (0.507937*pct, 0.581914, 0.581914),(0.523810*pct, 0.590937, 0.590937),
                        (0.539683*pct, 0.599824, 0.599824),(0.555556*pct, 0.608581, 0.608581),
                        (0.571429*pct, 0.617213, 0.617213),(0.587302*pct, 0.625727, 0.625727),
                        (0.603175*pct, 0.634126, 0.634126),(0.619048*pct, 0.642416, 0.642416),
                        (0.634921*pct, 0.650600, 0.650600),(0.650794*pct, 0.658682, 0.658682),
                        (0.666667*pct, 0.666667, 0.666667),(0.682540*pct, 0.674556, 0.674556),
                        (0.698413*pct, 0.682355, 0.682355),(0.714286*pct, 0.690066, 0.690066),
                        (0.730159*pct, 0.697691, 0.697691),(0.746032*pct, 0.705234, 0.705234),
                        (0.761905*pct, 0.727166, 0.727166),(0.777778*pct, 0.748455, 0.748455),
                        (0.793651*pct, 0.769156, 0.769156),(0.809524*pct, 0.789314, 0.789314),
                        (0.825397*pct, 0.808969, 0.808969),(0.841270*pct, 0.828159, 0.828159),
                        (0.857143*pct, 0.846913, 0.846913),(0.873016*pct, 0.865261, 0.865261),
                        (0.888889*pct, 0.883229, 0.883229),(0.904762*pct, 0.900837, 0.900837),
                        (0.920635*pct, 0.918109, 0.918109),(0.936508*pct, 0.935061, 0.935061),
                        (0.952381*pct, 0.951711, 0.951711),(0.968254*pct, 0.968075, 0.968075),
                        (0.984127*pct, 0.984167, 0.984167),(1.0*pct, 1.0, 1.0),
            ((1-pct)*0.125+pct, 0.9686274528503418,
            0.9686274528503418), ((1-pct)*0.25+pct, 0.93725490570068359, 0.93725490570068359),
            ((1-pct)*0.375+pct, 0.88235294818878174, 0.88235294818878174), ((1-pct)*0.5+pct,
            0.83921569585800171, 0.83921569585800171), ((1-pct)*0.625+pct, 0.7764706015586853,
            0.7764706015586853), ((1-pct)*0.75+pct, 0.70980393886566162, 0.70980393886566162),
            ((1-pct)*0.875+pct, 0.61176472902297974, 0.61176472902297974), (1.0,
            0.41960784792900085, 0.41960784792900085)]}

   PkBl_coltbl = LinearSegmentedColormap('PKBL_COLTBL',PkBl_data)
   return PkBl_coltbl



def PuRdBlfloat(datarange):
   distance = max(datarange)-min(datarange)
   zeroloc = 0-min(datarange)
   pct = float(zeroloc)/float(distance)


   PuRdBl_data = {'blue': [
            (0.0*pct, 0.12156862765550613,0.12156862765550613),
            (0.125*pct,0.26274511218070984, 0.26274511218070984), 
            (0.25*pct,0.33725491166114807, 0.33725491166114807), 
            (0.375*pct, 0.54117649793624878, 0.54117649793624878), 
            (0.5*pct, 0.69019609689712524, 0.69019609689712524),
            (0.625*pct, 0.78039216995239258,0.78039216995239258), 
            (0.75*pct, 0.85490196943283081,0.85490196943283081), 
            (0.875*pct, 0.93725490570068359,0.93725490570068359), 
            (1.0*pct, 0.97647058963775635,0.97647058963775635), 
            ((1-pct)*0.125+pct, 0.9686274528503418,
            0.9686274528503418), ((1-pct)*0.25+pct, 0.93725490570068359, 0.93725490570068359),
            ((1-pct)*0.375+pct, 0.88235294818878174, 0.88235294818878174), ((1-pct)*0.5+pct,
            0.83921569585800171, 0.83921569585800171), ((1-pct)*0.625+pct, 0.7764706015586853,
            0.7764706015586853), ((1-pct)*0.75+pct, 0.70980393886566162, 0.70980393886566162),
            ((1-pct)*0.875+pct, 0.61176472902297974, 0.61176472902297974), (1.0,
            0.41960784792900085, 0.41960784792900085)],



            'green': [
            (0.0*pct, 0.0, 0.0),
            (0.125*pct, 0.0, 0.0),
            (0.25*pct, 0.070588238537311554, 0.070588238537311554), 
            (0.375*pct, 0.16078431904315948, 0.16078431904315948), 
            (0.5*pct, 0.3960784375667572, 0.3960784375667572), 
            (0.625*pct, 0.58039218187332153, 0.58039218187332153), 
            (0.75*pct, 0.72549021244049072, 0.72549021244049072), 
            (0.875*pct,0.88235294818878174, 0.88235294818878174), 
            (1.0*pct, 0.95686274766921997, 0.95686274766921997), 
            ((1-pct)*0.125+pct,    0.92156863212585449, 0.92156863212585449), ((1-pct)*0.25+pct,
                0.85882353782653809, 0.85882353782653809), ((1-pct)*0.375+pct,
                0.7921568751335144, 0.7921568751335144), ((1-pct)*0.5+pct,
                0.68235296010971069, 0.68235296010971069), ((1-pct)*0.625+pct,
                0.57254904508590698, 0.57254904508590698), ((1-pct)*0.75+pct,
                0.44313725829124451, 0.44313725829124451), ((1-pct)*0.875+pct,
                0.31764706969261169, 0.31764706969261169), ((1-pct)*1.0+pct,
                0.18823529779911041, 0.18823529779911041)],



            'red': [
            (0.0*pct, 0.40392157435417175, 0.40392157435417175),
            (0.125*pct, 0.59607845544815063, 0.59607845544815063), 
            (0.25*pct, 0.80784314870834351, 0.80784314870834351), 
            (0.375*pct, 0.90588235855102539, 0.90588235855102539), 
            (0.5*pct, 0.87450981140136719, 0.87450981140136719), 
            (0.625*pct, 0.78823530673980713, 0.78823530673980713), 
            (0.75*pct, 0.83137255907058716, 0.83137255907058716), 
            (0.875*pct, 0.90588235855102539, 0.90588235855102539), 
            (1.0*pct, 0.9686274528503418, 0.9686274528503418), 
            ((1-pct)*0.125+pct,
                0.87058824300765991, 0.87058824300765991), ((1-pct)*0.25+pct,
                0.7764706015586853, 0.7764706015586853), ((1-pct)*0.375+pct,
                0.61960786581039429, 0.61960786581039429), ((1-pct)*0.5+pct,
                0.41960784792900085, 0.41960784792900085), ((1-pct)*0.625+pct,
                0.25882354378700256, 0.25882354378700256), ((1-pct)*0.75+pct,
                0.12941177189350128, 0.12941177189350128), ((1-pct)*0.875+pct,
                0.031372550874948502, 0.031372550874948502), (1.0,
                0.031372550874948502, 0.031372550874948502)]}

   #print PuRdBl_data
   PuRdBl_coltbl = LinearSegmentedColormap('PURDBL_COLTBL',PuRdBl_data)
   return PuRdBl_coltbl


def RdBuWH():
   RdBuWH_data = {'blue': [(0.0, 0.3803921639919281, 0.3803921639919281),
                          (0.10000000000000002, 0.67450982332229614, 0.67450982332229614), 
                          (0.20000000000000004, 0.76470589637756348, 0.76470589637756348), 
                          (0.29999999999999996, 0.87058824300765991, 0.87058824300765991),
                          (0.39999999999999998, 0.94117647409439087, 0.94117647409439087),
                          (0.45, 1.0, 1.0),
                          (0.55, 1.0, 1.0),
                          (0.60000000000000002, 0.78039216995239258, 0.78039216995239258), 
                          (0.69999999999999999, 0.50980395078659058, 0.50980395078659058), 
                          (0.80000000000000001, 0.30196079611778259, 0.30196079611778259), 
                          (0.90000000000000001, 0.16862745583057404, 0.16862745583057404), 
                          (1.0, 0.12156862765550613, 0.12156862765550613)], 

                 'green': [(0.0, 0.18823529779911041, 0.18823529779911041),
                          (0.10000000000000002, 0.40000000596046448, 0.40000000596046448), 
                          (0.20000000000000004, 0.57647061347961426, 0.57647061347961426), 
                          (0.29999999999999996, 0.77254903316497803, 0.77254903316497803), 
                          (0.39999999999999998, 0.89803922176361084, 0.89803922176361084), 
                          (0.45, 1.0, 1.0),
                          (0.55, 1.0, 1.0),
                          (0.60000000000000002, 0.85882353782653809, 0.85882353782653809), 
                          (0.69999999999999999, 0.64705884456634521, 0.64705884456634521), 
                          (0.80000000000000001, 0.37647059559822083, 0.37647059559822083), 
                          (0.90000000000000001, 0.094117648899555206, 0.094117648899555206), 
                          (1.0, 0.0, 0.0)],

                   'red': [(0.0, 0.019607843831181526, 0.019607843831181526),
                          (0.10000000000000002, 0.12941177189350128, 0.12941177189350128),
                          (0.20000000000000004, 0.26274511218070984, 0.26274511218070984),
                          (0.29999999999999996, 0.57254904508590698, 0.57254904508590698),
                          (0.39999999999999998, 0.81960785388946533, 0.81960785388946533),
                          (0.45, 1.0, 1.0),
                          (0.55, 1.0, 1.0),
                          (0.60000000000000002, 0.99215686321258545, 0.99215686321258545),
                          (0.69999999999999999, 0.95686274766921997, 0.95686274766921997),
                          (0.80000000000000001, 0.83921569585800171, 0.83921569585800171),
                          (0.90000000000000001, 0.69803923368453979, 0.69803923368453979),
                          (1.0, 0.40392157435417175, 0.40392157435417175)]}
   RdBuWH_coltbl = LinearSegmentedColormap('RDBUWH_COLTBL',RdBuWH_data)
   return RdBuWH_coltbl

