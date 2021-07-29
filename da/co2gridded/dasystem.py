"""CarbonTracker Data Assimilation Shell (CTDAS) Copyright (C) 2017 Wouter Peters. 
Users are recommended to contact the developers (wouter.peters@wur.nl) to receive
updates of the code. See also: http://www.carbontracker.eu. 

This program is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software Foundation, 
version 3. This program is distributed in the hope that it will be useful, but 
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. 

You should have received a copy of the GNU General Public License along with this 
program. If not, see <http://www.gnu.org/licenses/>."""
#!/usr/bin/env python
# control.py

"""
Author : peters 

Revision History:
File created on 26 Aug 2010.

"""

import logging


################### Begin Class CtDaSystem ###################

from da.baseclasses.dasystem import DaSystem

class CO2GriddedDaSystem(DaSystem):
    """ Information on the data assimilation system used. This is normally an rc-file with settings.
    """

    def __init__(self, rcfilename):
        """
        Initialization occurs from passed rc-file name, items in the rc-file will be added
        to the dictionary
        """

        self.ID = 'CarbonTracker Gridded CO2'    # the identifier gives the platform name
        self.load_rc(rcfilename)

        logging.debug('Data Assimilation System initialized: %s' % self.ID)

    def validate(self):
        """ 
        validate the contents of the rc-file given a dictionary of required keys
        """

        needed_rc_items = ['obs.input.dir',
                           'obs.input.fname',
                           'ocn.covariance',
                           'nparameters',
                           'deltaco2.prefix',
                           'regtype']


        for k, v in self.iteritems():
            if v == 'True' : self[k] = True
            if v == 'False': self[k] = False

        for key in needed_rc_items:
            if key not in self:
                msg = 'Missing a required value in rc-file : %s' % key
                logging.error(msg)
                raise IOError(msg)

        logging.debug('DA System Info settings have been validated succesfully')

################### End Class CtDaSystem ###################


if __name__ == "__main__":
    pass
