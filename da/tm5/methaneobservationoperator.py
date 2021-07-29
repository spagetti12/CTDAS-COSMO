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
# tm5_tools.py

"""
Author : Aki 

File created on Feb 2010.
"""

import os
import sys
import logging
import shutil
import datetime
import subprocess
from string import join
import glob
sys.path.append(os.getcwd())
sys.path.append("../../")

import da.tools.rc as rc
from da.tools.general import create_dirs, to_datetime
from da.baseclasses.observationoperator import ObservationOperator
from da.tm5.observationoperator import TM5ObservationOperator

def validate_rc_methane(self):
    """
    Validate the contents of the tm_settings dictionary and add extra values. The required items for the TM5 rc-file
    are specified in the tm5_tools module, as dictionary variable "needed_rc_items".
    """

    if self.rcfiletype == 'pycasso':
        self.projectkey = 'my.project.dir'
        self.rundirkey = 'my.run.dir'
        self.outputdirkey = 'output.dir'
        self.savedirkey = 'restart.write.dir'
        self.timestartkey = 'timerange.start'
        self.timefinalkey = 'timerange.end'
        self.timelengthkey = 'jobstep.length'
        self.istartkey = 'istart'
        self.restartvalue = 33
        self.coldstartvalue = 5  #Aki: read from mmix file
        self.transitionvalue = 61
        self.mysourcedirs = 'my.source.dirs'
        self.addedsourcedirs = '  proj/levels/ml91/trunk proj/levels/ml91/tropo34/trunk'
    else:
        self.projectkey = 'runid'
        self.rundirkey = 'rundir'
        self.outputdirkey = 'outputdir'
        self.savedirkey = 'savedir'
        self.timestartkey = 'time.start'
        self.timefinalkey = 'time.final'
        self.timelengthkey = 'time.break.nday'
        self.istartkey = 'istart'
        self.restartvalue = 3
        self.coldstartvalue = 9
        self.transitionvalue = 61

    needed_rc_items = [
                        self.projectkey,
                        self.rundirkey,
                        self.outputdirkey,
                        self.savedirkey,
                        self.timestartkey,
                        self.timefinalkey,
                        self.timelengthkey,
                        self.istartkey
                      ]

    for k, v in self.tm_settings.items():
        if v == 'True' : self.tm_settings[k] = True
        if v == 'False': self.tm_settings[k] = False
        if 'date' in k : self.tm_settings[k] = to_datetime(v)
        if 'time.start' in k : 
            self.tm_settings[k] = to_datetime(v, fmt='TM5')
        if 'time.final' in k : 
            self.tm_settings[k] = to_datetime(v, fmt='TM5')
        if 'timerange.start' in k : 
            self.tm_settings[k] = to_datetime(v)
        if 'timerange.end' in k : 
            self.tm_settings[k] = to_datetime(v)

    for key in needed_rc_items:
        if key not in self.tm_settings:
            msg = 'Missing a required value in rc-file : %s' % key
            logging.error(msg)
            raise IOError(msg)
    logging.debug('rc-file has been validated succesfully')

TM5ObservationOperator.validate_rc = validate_rc_methane

