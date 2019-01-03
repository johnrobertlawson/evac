""" Create a catalogue of objects, either in numpy array or dictionary format.

Catalogue may/may not allow for matching between different sets, so subcataloging?

May just be a function that returns a dataframe or labelled numpy array.
"""

class CombinedObjects:
    def __init__(self,):
        pass

    def compute_total_interest(self,propA=None,propB=None,cd=None,md=None,td=None,
                            cd_max=40.0,md_max=40.0,td_max=25.0):
        """
        Returns the total interest between object pairs.

        Args:
            propA, propB: dictionary? struct array? panda dataframe? slice of
                one object's properties

            cd (float): the centroid distance between two objects (km)
            md (float): the minimum distance between two objects (km)
            td (float): the time difference between the objects (min)

            cd_max (float): the limit for centroid distance (km)
            md_max (float): the limit for minimum distance (km)
        """
        # Sanity checks
        if propA is not None:
            assert propB is not None
            assert not all(cd,md,td)
            method = 1
        else:
            assert propB is None
            assert all(cd,md,td)
            method = 2
        assert all(cd_max,md_max,td_max)

        if method == 1:
            # Compute cd, md, td

        TI = ((td_max-td)/t_max)*0.5*(((cd_max-cd)/cd_max)+((md_max-md)/md_max))
        return TI
