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

#################################################################################################
# First order of business is always to make all other python modules accessible through the path
#################################################################################################

import sys
import os
import logging
sys.path.append(os.getcwd())

#################################################################################################
# Next, import the tools needed to initialize a data assimilation cycle
#################################################################################################

from da.tools.initexit import start_logger, validate_opts_args, parse_options, CycleControl
from da.tools.pipeline import ensemble_smoother_pipeline, header, footer
from da.platform.jet import JetPlatform 
from da.carbondioxide.dasystem import CO2DaSystem 
from da.carbondioxide.statevector import CO2StateVector 
from da.carbondioxide.obs import CO2Observations 
from da.tm5.observationoperator import TM5ObservationOperator 
from da.carbondioxide.optimizer import CO2Optimizer


#################################################################################################
# Parse and validate the command line options, start logging
#################################################################################################

start_logger()
opts, args = validate_opts_args(parse_options())

#################################################################################################
# Create the Cycle Control object for this job    
#################################################################################################

dacycle = CycleControl(opts, args)

platform = JetPlatform()
dasystem = CO2DaSystem(dacycle['da.system.rc'])
obsoperator = TM5ObservationOperator(dacycle['da.obsoperator.rc'])
samples = CO2Observations()
statevector = CO2StateVector()
optimizer = CO2Optimizer()

##########################################################################################
################### ENTER THE PIPELINE WITH THE OBJECTS PASSED BY THE USER ###############
##########################################################################################


logging.info(header + "Entering Pipeline " + footer) 

ensemble_smoother_pipeline(dacycle, platform, dasystem, samples, statevector, obsoperator, optimizer)


##########################################################################################
################### All done
##########################################################################################

