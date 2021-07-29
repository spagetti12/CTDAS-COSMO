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
from da.baseclasses.statevector import StateVector, EnsembleMember
import da.tools.io4 as io
from datetime import timedelta

identifier = 'SF6 Statevector '
version = '0.0'

################### Begin Class SF6StateVector ###################

class SF6StateVector(StateVector):
    """ This is a StateVector object for CarbonTracker. It has a private method to make new ensemble members """

    def get_covariance(self, date, dacycle):
        """ Make a new ensemble from specified matrices, the attribute lag refers to the position in the state vector. 
            Note that lag=1 means an index of 0 in python, hence the notation lag-1 in the indexing below.
            The argument is thus referring to the lagged state vector as [1,2,3,4,5,..., nlag]
        """    
        try:
            import matplotlib.pyplot as plt
        except:
            pass

        # Get the needed matrices from the specified covariance files

        fullcov = np.zeros((self.nparams, self.nparams), float)

        for n in range(self.nparams):
            fullcov[n,n] = 0.5**2

        return fullcov

    def setup(self, dacycle):
        """
        Initialize the object by specifying the dimensions. 
        There are two major requirements for each statvector that you want to build:
        
            (1) is that the statevector can map itself onto a regular grid
            (2) is that the statevector can map itself (mean+covariance) onto TransCom regions

        An example is given below.
        """

        self.nlag = int(dacycle['time.nlag'])
        self.nmembers = int(dacycle['da.optimizer.nmembers'])
        self.nparams = int(dacycle.dasystem['nparameters'])
        self.nobs = 0
        
        self.obs_to_assimilate = ()  # empty containter to hold observations to assimilate later on

        # These list objects hold the data for each time step of lag in the system. Note that the ensembles for each time step consist 
        # of lists of EnsembleMember objects, we define member 0 as the mean of the distribution and n=1,...,nmembers as the spread.

        self.ensemble_members = range(self.nlag)

        for n in range(self.nlag):
            self.ensemble_members[n] = []


        # This specifies the file to read with the gridded mask at 1x1 degrees. Each gridbox holds a number that specifies the parametermember
        #  that maps onto it. From this map, a dictionary is created that allows a reverse look-up so that we can map parameters to a grid.

        mapfile = os.path.join(dacycle.dasystem['regionsfile'])
        ncf = io.ct_read(mapfile, 'read')
        self.tcmap = ncf.get_variable('transcom_regions')
        ncf.close()

        self.gridmap = np.ones((180,360),'float')

        logging.debug("A TransCom  map on 1x1 degree was read from file %s" % dacycle.dasystem['regionsfile'])
        logging.debug("A parameter map on 1x1 degree was created")

        # Create a dictionary for state <-> gridded map conversions

        nparams = self.gridmap.max()
        self.griddict = {}
        for r in range(1, int(nparams) + 1):
            sel = (self.gridmap.flat == r).nonzero()
            if len(sel[0]) > 0: 
                self.griddict[r] = sel

        logging.debug("A dictionary to map grids to states and vice versa was created")

        # Create a matrix for state <-> TransCom conversions

        self.tcmatrix = np.zeros((self.nparams, 23), 'float') 

        logging.debug("A matrix to map states to TransCom regions and vice versa was created")

        # Create a mask for species/unknowns

        self.make_species_mask()

    def make_species_mask(self):

        """

        This method creates a dictionary with as key the name of a tracer, and as values an array of 0.0/1.0 values 
        specifying which StateVector elements are constrained by this tracer. This mask can be used in 
        the optimization to ensure that certain types of osbervations only update certain unknowns.

        An example would be that the tracer '14CO2' can be allowed to only map onto fossil fuel emissions in the state

        The form of the mask is:

        {'co2': np.ones(self.nparams), 'co2c14', np.zeros(self.nparams)  }

        so that 'co2' maps onto all parameters, and 'co2c14' on none at all. These arrays are used in the Class 
        optimizer when state updates are actually performed

        """
        self.speciesdict = {'sf6': np.ones(self.nparams)}
        logging.debug("A species mask was created, only the following species are recognized in this system:")
        for k in self.speciesdict.keys(): 
            logging.debug("   ->    %s" % k)

    def propagate(self,dummy):
        """ 
        :rtype: None

        Propagate the parameter values in the StateVector to the next cycle. This means a shift by one cycle 
        step for all states that will
        be optimized once more, and the creation of a new ensemble for the time step that just 
        comes in for the first time (step=nlag). 
        In the future, this routine can incorporate a formal propagation of the statevector.

        """

        self.ensemble_members.append([])
        cov = self.get_covariance(None,None)
        newmean = np.ones(self.nparams)
        self.make_new_ensemble(self.nlag+1, newmean = newmean, covariancematrix=cov)
        self.ensemble_members.pop(0)

        logging.info('The state vector remains the same in the SF6 run')
        logging.info('The state vector has been propagated by one cycle')

    def make_new_ensemble(self, lag, newmean=None, covariancematrix=None):
        """ 
        :param lag: an integer indicating the time step in the lag order
        :param covariancematrix: a matrix to draw random values from
        :rtype: None
    
        Make a new ensemble, the attribute lag refers to the position in the state vector. 
        Note that lag=1 means an index of 0 in python, hence the notation lag-1 in the indexing below.
        The argument is thus referring to the lagged state vector as [1,2,3,4,5,..., nlag]

        The optional covariance object to be passed holds a matrix of dimensions [nparams, nparams] which is
        used to draw ensemblemembers from. If this argument is not passed it will ne substituted with an 
        identity matrix of the same dimensions.

        """    

        if newmean == None: 
            newmean = np.ones(self.nparams)

        if covariancematrix == None: 
            covariancematrix = np.identity(self.nparams)

        # Make a cholesky decomposition of the covariance matrix


        _, s, _ = np.linalg.svd(covariancematrix)
        dof = np.sum(s) ** 2 / sum(s ** 2)
        C = np.linalg.cholesky(covariancematrix)

        logging.debug('Cholesky decomposition has succeeded ')
        logging.info('Appr. degrees of freedom in covariance matrix is %s' % (int(dof)))


        # Create the first ensemble member with a deviation of 0.0 and add to list

        newmember = EnsembleMember(0)
        newmember.param_values = newmean.flatten()  # no deviations
        self.ensemble_members[lag-1].append(newmember)

        # Create members 1:nmembers and add to EnsembleMembers list

        for member in range(1, self.nmembers):
            rands = np.random.randn(self.nparams)

            newmember = EnsembleMember(member)
            # routine to avoids that members < 0.0
            dev   = np.dot(C, rands)                        #VDV
            dummy = np.zeros(len(dev))                      #VDV   
            for i in range(len(dev)):                       #VDV
                if dev[i] < 0.0:                            #VDV 
                    dummy[i] = newmean[i]*np.exp(dev[i])    #VDV
                else:                                       #VDV
                    dummy[i] = newmean[i]*(1+dev[i])        #VDV
            newmember.param_values = dummy               #VDV
            #newmember.ParameterValues = np.dot(C, rands) + newmean 
            self.ensemble_members[lag-1].append(newmember)

        logging.debug('%d new ensemble members were added to the state vector # %d' % (self.nmembers, lag))

    def write_members_to_file(self, lag, outdir):
        """ 
           :param: lag: Which lag step of the filter to write, must lie in range [1,...,nlag]
           :rtype: None

           Write ensemble member information to a NetCDF file for later use. The standard output filename is 
           *parameters.DDD.nc* where *DDD* is the number of the ensemble member. Standard output file location 
           is the `dir.input` of the dacycle object. In principle the output file will have only two datasets inside 
           called `parametervalues` which is of dimensions `nparameters` and `parametermap` which is of dimensions (180,360). 
           This dataset can be read and used by a :class:`~da.baseclasses.observationoperator.ObservationOperator` object. 

           .. note:: if more, or other information is needed to complete the sampling of the ObservationOperator you
                     can simply inherit from the StateVector baseclass and overwrite this write_members_to_file function.

        """

        # These import statements caused a crash in netCDF4 on MacOSX. No problems on Jet though. Solution was
        # to do the import already at the start of the module, not just in this method.
           
        #import da.tools.io as io
        #import da.tools.io4 as io

        members = self.ensemble_members[lag]

        for mem in members:
            filename = os.path.join(outdir, 'parameters.%03d.nc' % mem.membernumber)
            ncf = io.CT_CDF(filename, method='create')
            dimparams = ncf.add_params_dim(self.nparams)
            dimgrid = ncf.add_latlon_dim()
            logging.warning('Parameters for TM5 are NOT taken with exponential')

            # Explicitly set maximum allowable value to 0.0
            data = np.where(mem.param_values < 0, 0.0, mem.param_values)

            savedict = io.std_savedict.copy()
            savedict['name'] = "parametervalues"
            savedict['long_name'] = "parameter_values_for_member_%d" % mem.membernumber
            savedict['units'] = "unitless"
            savedict['dims'] = dimparams 
            savedict['values'] = data
            savedict['comment'] = 'These are parameter values to use for member %d, note: they result from an exponential function' % mem.membernumber
            ncf.add_data(savedict)

            griddata = self.vector2grid(vectordata=data)

            savedict = io.std_savedict.copy()
            savedict['name'] = "parametermap"
            savedict['long_name'] = "parametermap_for_member_%d" % mem.membernumber
            savedict['units'] = "unitless"
            savedict['dims'] = dimgrid 
            savedict['values'] = griddata.tolist()
            savedict['comment'] = 'These are gridded parameter values to use for member %d, note: they result from an exponential function' % mem.membernumber
            ncf.add_data(savedict)

            ncf.close()

            logging.debug('Successfully wrote data from ensemble member %d to file (%s) ' % (mem.membernumber, filename))




################### End Class SF6StateVector ###################

