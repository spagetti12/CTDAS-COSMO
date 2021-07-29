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
# io.py

"""
Author : aki

Revision History:
File created on Apr 2016.

"""
import standardvariables
import datetime as dt
from numpy import array, arange
import os
import logging
import sys
sys.path.append('/stornext/store/carbon/CarbonTracker/netcdf4/python_packages_install_dir/lib/python2.7/site-packages/')
import netCDF4
sys.path.append('/stornext/store/carbon/CarbonTracker/pyhdf/pyhdf-0.8.3/lib/python2.7/site-packages/')
import pyhdf.SD as hdf
sys.path.append('../')
from da.tools.io4 import *


disclaimer = "This data belongs to the CarbonTracker project"
email = "aki.tsuruta@fmi.fi"
url = "http://en.ilmatieteenlaitos.fi/carbon-cycle-modelling"
institution = "Finnish Meteorological Institute, Climate research"
source 	= "CarbonTracker release 1.0" 
conventions = "CF-1.1"
historytext	= 'created on '+dt.datetime.now().strftime('%B %d, %Y')+' by %s'%os.environ['USER']


def add_tc_header2(self):

    #
    self.setncattr('Institution',institution)
    self.setncattr('Contact',email)
    self.setncattr('URL',url)
    self.setncattr('Source',source)
    self.setncattr('Convention',conventions)
    self.setncattr('Disclaimer',disclaimer)
    self.setncattr('History',historytext)

def standard_var2(self,varname):
    """ return properties of standard variables """
    import standardvariables
    
    if varname in standardvariables.standard_variables.keys():
        return standardvariables.standard_variables[varname]
    else:
        return standardvariables.standard_variables['unknown']



CT_CDF.add_tc_header = add_tc_header2
CT_CDF.standard_var = standard_var2
CT_HDF.standard_var = standard_var2
