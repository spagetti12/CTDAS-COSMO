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
Author : aki

Revision History:
File created on 16 Nov 2013.

"""

import logging

################### Begin Class MethaneDaSystem ###################

from da.baseclasses.dasystem import DaSystem

class MethaneDaSystem(DaSystem):
    """ Information on the data assimilation system used. This is normally an rc-file with settings.
    """
    def validate(self):
        """ 
        validate the contents of the rc-file given a dictionary of required keys
        """

        needed_rc_items = ['obs.input.dir',
                           'obs.input.fname']


        for k, v in self.iteritems():
            if v == 'True' : 
                self[k] = True
            if v == 'False': 
                self[k] = False

        for key in needed_rc_items:
            if key not in self:
                logging.warning('Missing a required value in rc-file : %s' % key)
        logging.debug('DA System Info settings have been validated succesfully')

################### End Class MethaneDaSystem ###################


if __name__ == "__main__":
    pass
