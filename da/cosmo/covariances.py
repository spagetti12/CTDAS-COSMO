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
# ct_statevector_tools.py

"""
Author : peters 

Revision History:
File created on 28 Jul 2010.

"""

import os
import sys
sys.path.append(os.getcwd())

import logging
import numpy as np
#from da.cosmo.statevector_uniform import StateVector, EnsembleMember
#from da.cosmo.statevector_read_from_output import StateVector, EnsembleMember
#from da.cosmo.statevector_mean import StateVector, EnsembleMember
#from da.cosmo.statevector import StateVector, EnsembleMember
from da.cosmo.statevector_mean import StateVector, EnsembleMember

import da.tools.io4 as io

identifier = 'CarbonTracker Statevector '
version = '0.0'

np.set_printoptions(threshold=sys.maxsize)

################### Begin Class CO2StateVector ###################

class CO2StateVector(StateVector):
    """ This is a StateVector object for CarbonTracker. It has a private method to make new ensemble members """

    def get_covariance(self, date, dacycle):
        """ Make a new ensemble from specified matrices, the attribute lag refers to the position in the state vector. 
            Note that lag=1 means an index of 0 in python, hence the notation lag-1 in the indexing below.
            The argument is thus referring to the lagged state vector as [1,2,3,4,5,..., nlag]

                                    0. Needleleaf Evergreen, Temperate
                                    1. Needleleaf Evergreen, Boreal
                                    2. Boardleaf Decidous, Temperate
                                    3. Boardleaf Decidous, Boreal
                                    4. Boardleaf Decidous Shrub, Temperate
                                    5. Boardleaf Decidous Shrub, Boreal
                                    6. C3 Arctic Grass
                                    7. C3 non-Arctic Grass
                                    8. C4 Grass
                                    9. Crop

"""

        fullcov = np.zeros(shape=(90,90))

#        partcov = np.array([ \
 #       (0.64, 0.36, 0.16, 0.16, 0.16, 0.16, 0.04, 0.04, 0.04, 0.01), \
  #      (0.36, 0.64, 0.16, 0.16, 0.16, 0.16, 0.04, 0.04, 0.04, 0.01), \
   #     (0.16, 0.16, 0.64, 0.36, 0.16, 0.16, 0.04, 0.04, 0.04, 0.01), \
    #    (0.16, 0.16, 0.36, 0.64, 0.16, 0.16, 0.04, 0.04, 0.04, 0.01), \
     #   (0.16, 0.16, 0.16, 0.16, 0.64, 0.36, 0.04, 0.04, 0.04, 0.01), \
      #  (0.16, 0.16, 0.16, 0.16, 0.36, 0.64, 0.04, 0.04, 0.04, 0.01), \
       # (0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.64, 0.16, 0.16, 0.16), \
        #(0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.16, 0.64, 0.16, 0.16), \
#        (0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.16, 0.16, 0.64, 0.16), \
 #       (0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.16, 0.16, 0.16, 0.64) ])

# 0.1089 subsituted with 0.64
        partcov = np.array([ \
        (0.64, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4), \
        (0.4, 0.64, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4), \
        (0.4, 0.4, 0.64, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4), \
        (0.4, 0.4, 0.4, 0.64, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4), \
        (0.4, 0.4, 0.4, 0.4, 0.64, 0.4, 0.4, 0.4, 0.4, 0.4), \
        (0.4, 0.4, 0.4, 0.4, 0.4, 0.64, 0.4, 0.4, 0.4, 0.4), \
        (0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.64, 0.4, 0.4, 0.4), \
        (0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.64, 0.4, 0.4), \
        (0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.64, 0.4), \
        (0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.64) ])

        partcov_dg = np.array([ \
        (0.64, 0., 0., 0., 0., 0., 0., 0., 0., 0.), \
        (0., 0.64, 0., 0., 0., 0., 0., 0., 0., 0.), \
        (0., 0., 0.64, 0., 0., 0., 0., 0., 0., 0.), \
        (0., 0., 0., 0.64, 0., 0., 0., 0., 0., 0.), \
        (0., 0., 0., 0., 0.64, 0., 0., 0., 0., 0.), \
        (0., 0., 0., 0., 0., 0.64, 0., 0., 0., 0.), \
        (0., 0., 0., 0., 0., 0., 0.64, 0., 0., 0.), \
        (0., 0., 0., 0., 0., 0., 0., 0.64, 0., 0.), \
        (0., 0., 0., 0., 0., 0., 0., 0., 0.64, 0.), \
        (0., 0., 0., 0., 0., 0., 0., 0., 0., 0.64) ])

# L = 300 km
#        L_matrix = np.array([\
 #       (1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000), \
  #      (0.000, 1.000, 0.715, 0.650, 0.733, 0.614, 0.909, 0.688, 0.544), \
   #     (0.000, 0.715, 1.000, 0.862, 0.751, 0.638, 0.695, 0.501, 0.575), \
    #    (0.000, 0.650, 0.862, 1.000, 0.776, 0.696, 0.651, 0.472, 0.639), \
     #   (0.000, 0.733, 0.751, 0.776, 1.000, 0.827, 0.774, 0.587, 0.732), \
      #  (0.000, 0.614, 0.638, 0.696, 0.827, 1.000, 0.660, 0.537, 0.885), \
       # (0.000, 0.909, 0.695, 0.651, 0.774, 0.660, 1.000, 0.721, 0.586), \
        #(0.000, 0.688, 0.501, 0.472, 0.587, 0.537, 0.721, 1.000, 0.489), \
#        (0.000, 0.544, 0.575, 0.639, 0.732, 0.885, 0.586, 0.489, 1.000) ])

# L = 100 km
        L_matrix = np.array([\
        (1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000), \
        (0.   , 1.   , 0.134, 0.075, 0.155, 0.054, 0.563, 0.106, 0.026), \
        (0.   , 0.134, 1.   , 0.411, 0.179, 0.068, 0.112, 0.016, 0.036), \
        (0.   , 0.075, 0.411, 1.   , 0.218, 0.114, 0.076, 0.011, 0.068), \
        (0.   , 0.155, 0.179, 0.218, 1.   , 0.319, 0.216, 0.041, 0.154), \
        (0.   , 0.054, 0.068, 0.114, 0.319, 1.   , 0.083, 0.024, 0.48 ), \
        (0.   , 0.563, 0.112, 0.076, 0.216, 0.083, 1.   , 0.14 , 0.041), \
        (0.   , 0.106, 0.016, 0.011, 0.041, 0.024, 0.14 , 1.   , 0.014), \
        (0.   , 0.026, 0.036, 0.068, 0.154, 0.48 , 0.041, 0.014, 1.   ), ])

        zero_reg = ([3, 6])

        for i in range(9):
            for j in range(9):
                if i in zero_reg or j in zero_reg:
                    fullcov[i*10:(i+1)*10,j*10:(j+1)*10] = partcov_dg * L_matrix[i,j]
                else:
                    fullcov[i*10:(i+1)*10,j*10:(j+1)*10] = partcov * L_matrix[i,j]
#                fullcov[i*10:(i+1)*10,j*10:(j+1)*10] = partcov_cdgh * L_matrix[i,j]

        return fullcov

    def read_from_legacy_file(self, filename, qual='opt'):
        """ 
        :param filename: the full filename for the input NetCDF file
        :param qual: a string indicating whether to read the 'prior' or 'opt'(imized) StateVector from file
        :rtype: None

        Read the StateVector information from a NetCDF file and put in a StateVector object
        In principle the input file will have only one four datasets inside 
        called:
            * `meanstate_prior`, dimensions [nlag, nparamaters]
            * `ensemblestate_prior`, dimensions [nlag,nmembers, nparameters]
            * `meanstate_opt`, dimensions [nlag, nparamaters]
            * `ensemblestate_opt`, dimensions [nlag,nmembers, nparameters]

        This NetCDF information can be written to file using 
        :meth:`~da.baseclasses.statevector.StateVector.write_to_file`

        """
        

        f = io.ct_read(filename, 'read')

        for n in range(self.nlag):
            if qual == 'opt':
                meanstate = f.get_variable('statevectormean_opt')
                EnsembleMembers = f.get_variable('statevectorensemble_opt')

            elif qual == 'prior':
                meanstate = f.get_variable('statevectormean_prior')
                EnsembleMembers = f.get_variable('statevectorensemble_prior')

            if not self.ensemble_members[n] == []:
                self.ensemble_members[n] = []
                logging.warning('Existing ensemble for lag=%d was removed to make place for newly read data' % (n + 1))

            for m in range(self.nmembers):
                newmember = EnsembleMember(m)
                newmember.param_values = EnsembleMembers[m, :].flatten() + meanstate  # add the mean to the deviations to hold the full parameter values
        #        newmember.param_values[newmember.param_values<0.] = 0
                self.ensemble_members[n].append(newmember)

        f.close()

        logging.info('Successfully read the State Vector from file (%s) ' % filename)
    
    def read_from_file_exceptsam(self, filename, qual='opt'):
        """ 
        :param filename: the full filename for the input NetCDF file
        :param qual: a string indicating whether to read the 'prior' or 'opt'(imized) StateVector from file
        :rtype: None

        Read the StateVector information from a NetCDF file and put in a StateVector object
        In principle the input file will have only one four datasets inside 
        called:
            * `meanstate_prior`, dimensions [nlag, nparamaters]
            * `ensemblestate_prior`, dimensions [nlag,nmembers, nparameters]
            * `meanstate_opt`, dimensions [nlag, nparamaters]
            * `ensemblestate_opt`, dimensions [nlag,nmembers, nparameters]

        This NetCDF information can be written to file using 
        :meth:`~da.baseclasses.statevector.StateVector.write_to_file`

        """

        f = io.ct_read(filename, 'read')
        
        meanstate = f.get_variable('statevectormean_' + qual)
#        meanstate[:,39:77] = 1
        ensmembers = f.get_variable('statevectorensemble_' + qual)
        f.close()

        for n in range(self.nlag):
            if not self.ensemble_members[n] == []:
                self.ensemble_members[n] = []
                logging.warning('Existing ensemble for lag=%d was removed to make place for newly read data' % (n + 1))

            for m in range(self.nmembers):
                newmember = EnsembleMember(m)
                newmember.param_values = ensmembers[n, m, :].flatten() + meanstate[n]  # add the mean to the deviations to hold the full parameter values
                self.ensemble_members[n].append(newmember)

        logging.info('Successfully read the State Vector from file (%s) ' % filename)
#        logging.info('State Vector set to 1 for South American regions')

################### End Class CO2StateVector ###################


if __name__ == "__main__":
    pass
