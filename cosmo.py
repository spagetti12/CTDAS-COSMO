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
from da.cosmo.pipeline import ensemble_smoother_pipeline, header, footer
from da.platform.maunaloa import MaunaloaPlatform 
from da.baseclasses.dasystem import DaSystem 
from da.cosmo.covariances import CO2StateVector 
from da.cosmo.observations_cosmo import ObsPackObservations
from da.cosmo.optimizer import CO2Optimizer
from da.cosmo.observationoperator_cosmo import ObservationOperator
#from da.cosmo.expand_fluxes import save_weekly_avg_1x1_data, save_weekly_avg_state_data, save_weekly_avg_tc_data, save_weekly_avg_ext_tc_data
#from da.analysis.expand_molefractions import write_mole_fractions


#################################################################################################
# Parse and validate the command line options, start logging
#################################################################################################

start_logger()
#start_logger(level=logging.DEBUG)
opts, args = parse_options()
opts, args = validate_opts_args(opts, args)

#################################################################################################
# Create the Cycle Control object for this job    
#################################################################################################

dacycle = CycleControl(opts, args)


platform = MaunaloaPlatform()
dasystem = DaSystem(dacycle['da.system.rc'])
obsoperator = ObservationOperator(dacycle['da.obsoperator.rc'])
#samples     = Obs()
samples     = ObsPackObservations()
statevector = CO2StateVector()
optimizer = CO2Optimizer()

##########################################################################################
################### ENTER THE PIPELINE WITH THE OBJECTS PASSED BY THE USER ###############
##########################################################################################


logging.info(header + "Entering Pipeline " + footer) 

ensemble_smoother_pipeline(dacycle, platform, dasystem, samples, statevector, obsoperator, optimizer)


##########################################################################################
################### All done, extra stuff can be added next, such as analysis
##########################################################################################

logging.info(header + "All done. God bless" + footer) 

sys.exit(0)
save_weekly_avg_1x1_data(dacycle, statevector)
save_weekly_avg_state_data(dacycle, statevector)
save_weekly_avg_tc_data(dacycle, statevector)
save_weekly_avg_ext_tc_data(dacycle)
write_mole_fractions(dacycle)

sys.exit(0)
