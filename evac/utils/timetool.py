import os
import pdb

class TimeTool:
    """ A custom time object.

    Example:
        To convert from a number of seconds to the largest units
        possible, suitable for a WRF namelist for instance,
        try the following::
            
            secs = 10801
            TT = TimeTool(seconds=secs)
            d,h,m,s = TT.reduce()

        There!

    Note:
        * Not for dates - for that, use datetime.datetime.

    Args:
        years, months, days, hours, minutes, seconds (int): 
    """
    def __init__(self,days=None,hours=None,minutes=None,seconds=None,):
        self.days = days
        self.hours = hours
        self.minutes = minutes
        self.seconds = seconds

    def reduce(self):
        minutes,seconds = divmod(self.seconds,60)
        hours, minutes = divmod(minutes,60)
        days, hours = divmod(hours,24)
        return days, hours, minutes, seconds
