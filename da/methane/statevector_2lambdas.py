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
Author :  aki

Revision History:
File created on 18 Nov 2013.

"""

import os
import sys
sys.path.append(os.getcwd())

import logging
import numpy as np
from da.baseclasses.statevector import StateVector, EnsembleMember

from datetime import timedelta
import da.methane.io4 as io

identifier = 'CarbonTracker Statevector '
version = '0.0'

################### Begin Class CH4StateVector ###################

class MethaneStateVector(StateVector):
    """ This is a StateVector object for CarbonTracker. It has a private method to make new ensemble members """

    def setup(self, dacycle):
        """
        setup the object by specifying the dimensions. 
        There are two major requirements for each statvector that you want to build:
        
            (1) is that the statevector can map itself onto a regular grid
            (2) is that the statevector can map itself (mean+covariance) onto TransCom regions

        An example is given below.
        """

        self.nlag = int(dacycle['time.nlag'])
        self.nmembers = int(dacycle['da.optimizer.nmembers'])
        #self.distribution = dacycle['da.distribution']
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
        self.gridmap = ncf.get_variable('regions')
        self.tcmap = ncf.get_variable('transcom_regions')
        self.nparams = int(self.gridmap.max())
        self.nparams_bio = int(ncf.get_variable('nparams_bio'))
        ncf.close()

        logging.info ("Regional information read from file %s" % dacycle.dasystem['regionsfile'])
        logging.debug("A TransCom  map on 1x1 degree was read from file %s" % dacycle.dasystem['regionsfile'])
        logging.debug("A parameter map on 1x1 degree was read from file %s" % dacycle.dasystem['regionsfile'])

        # Create a dictionary for state <-> gridded map conversions

        nparams = self.gridmap.max()
        self.griddict = {}
        for r in range(1, int(nparams) + 1):
            sel = (self.gridmap.flat == r).nonzero()
            if len(sel[0]) > 0: 
                self.griddict[r] = sel

        logging.debug("A dictionary to map grids to states and vice versa was created")

        # Create a matrix for state <-> TransCom conversions

        self.tcmatrix = np.zeros((self.nparams, int(self.tcmap.max())), 'float') 

        for r in range(1, self.nparams + 1):
            sel = (self.gridmap.flat == r).nonzero()
            if len(sel[0]) < 1: 
                continue
            else:
                n_tc = set(self.tcmap.flatten().take(sel[0]))
                if len(n_tc) > 1: 
                    logging.error("Parameter %d seems to map to multiple TransCom regions (%s), I do not know how to handle this" % (r, n_tc))
                    raise ValueError
                self.tcmatrix[r - 1, n_tc.pop() - 1] = 1.0

        logging.debug("A matrix to map states to TransCom regions and vice versa was created")

        # Create a mask for species/unknowns

        self.make_species_mask()


    def make_species_mask(self):

        self.speciesdict = {'ch4': np.ones(self.nparams)}
        logging.debug("A species mask was created, only the following species are recognized in this system:")
        for k in self.speciesdict.keys(): 
            logging.debug("   ->    %s" % k)

    def get_covariance(self, date, dacycle):
        """ Make a new ensemble from specified matrices, the attribute lag refers to the position in the state vector. 
            Note that lag=1 means an index of 0 in python, hence the notation lag-1 in the indexing below.
            The argument is thus referring to the lagged state vector as [1,2,3,4,5,..., nlag]
        """    
        try:
            import matplotlib.pyplot as plt
        except:
            pass

        #----Arbitray covariance matrix ----#
        #fullcov = np.zeros((self.nparams, self.nparams), float)

        #for i in xrange(self.nparams):
        #    #fullcov[i,i] = 1.  # Identity matrix
        #    fullcov[i,i] = 0.08
        #fullcov[self.nparams-1,self.nparams-1] = 1.e-10

        #---- Covariance read from file -----#
        cov_file = dacycle.dasystem['covfile']
        nfcr = io.CT_CDF(cov_file, 'read')
        logging.info('covariance read from file: %s' %cov_file)
        fullcov = nfcr.variables['covariance_matrix'][:]

#        try:
#            plt.imshow(fullcov)
#            plt.colorbar()
#            plt.savefig('fullcovariancematrix.png')
#            plt.close('all')
#            logging.debug("Covariance matrix visualized for inspection")
#        except:
#            pass

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
                meanstate = f.get_variable('xac_%02d' % (n + 1))
                EnsembleMembers = f.get_variable('adX_%02d' % (n + 1))

            elif qual == 'prior':
                meanstate = f.get_variable('xpc_%02d' % (n + 1))
                EnsembleMembers = f.get_variable('pdX_%02d' % (n + 1))

            if not self.ensemble_members[n] == []:
                self.ensemble_members[n] = []
                logging.warning('Existing ensemble for lag=%d was removed to make place for newly read data' % (n + 1))

            for m in range(self.nmembers):
                newmember = EnsembleMember(m)
                newmember.param_values = EnsembleMembers[m, :].flatten() + meanstate  # add the mean to the deviations to hold the full parameter values
                self.ensemble_members[n].append(newmember)


        f.close()

        logging.info('Successfully read the State Vector from file (%s) ' % filename)

    def write_members_to_file(self, lag, outdir, dacycle):
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

        #--- regionfile ---#
        region_land_file = dacycle.dasystem['regionsfile']
        nfcr = io.CT_CDF(region_land_file, 'read')
        logging.debug('region land file read: %s' %region_land_file)
        bio_land = nfcr.variables['bio_land'][:]
        fossil_land = nfcr.variables['fossil_land'][:]


        for mem in members:
            filename = os.path.join(outdir, 'parameters.%03d.nc' % mem.membernumber)
            ncf = io.CT_CDF(filename, method='create')
            dimparams = ncf.add_params_dim(self.nparams)
            dimgrid = ncf.add_latlon_dim()

            data = mem.param_values

            savedict = io.std_savedict.copy()
            savedict['name'] = "parametervalues"
            savedict['long_name'] = "parameter_values_for_member_%d" % mem.membernumber
            savedict['units'] = "unitless"
            savedict['dims'] = dimparams 
            savedict['values'] = data
            savedict['comment'] = 'These are parameter values to use for member %d' % mem.membernumber
            ncf.add_data(savedict)

#            #--- bio ---#
#            dimparams_bio = ncf.add_params_dim(self.nparams_bio)
#            data_bio = data[0:self.nparams_bio]
#
#            savedict = io.std_savedict.copy()
#            savedict['name'] = "parametervalues_bio"
#            savedict['long_name'] = "parameter_bio_values_for_member_%d" % mem.membernumber
#            savedict['units'] = "unitless"
#            savedict['dims'] = dimparams_bio 
#            savedict['values'] = data_bio
#            savedict['comment'] = 'These are parameter bio values to use for member %d' % mem.membernumber
#            ncf.add_data(savedict)
#
#            #--- fossil ---#
#            nparams_fossil = self.nparams - self.nparams_bio
#            dimparams_fossil = ncf.add_params_dim(nparams_fossil)
#            data_fossil = data[self.nparams_bio:nparams_fossil]
#
#            savedict = io.std_savedict.copy()
#            savedict['name'] = "parametervalues_ff"
#            savedict['long_name'] = "parameter_ff_values_for_member_%d" % mem.membernumber
#            savedict['units'] = "unitless"
#            savedict['dims'] = dimparams_fossil 
#            savedict['values'] = data_fossil
#            savedict['comment'] = 'These are parameter fossil values to use for member %d' % mem.membernumber
#            ncf.add_data(savedict)

            #--- All parameters, gridded ---#
            griddata = self.vector2grid(vectordata=data)

            savedict = io.std_savedict.copy()
            savedict['name'] = "parametermap"
            savedict['long_name'] = "parametermap_for_member_%d" % mem.membernumber
            savedict['units'] = "unitless"
            savedict['dims'] = dimgrid 
            savedict['values'] = griddata.tolist()
            savedict['comment'] = 'These are gridded parameter values to use for member %d' % mem.membernumber
            ncf.add_data(savedict)

            #--- bio parameters, gridded ---#
            griddata_bio = griddata[:]
            w = np.where(bio_land==0)
            griddata_bio[w] = 1.0

            savedict = io.std_savedict.copy()
            savedict['name'] = "parametermap_bio"
            savedict['long_name'] = "parametermap_bio_for_member_%d" % mem.membernumber
            savedict['units'] = "unitless"
            savedict['dims'] = dimgrid 
            savedict['values'] = griddata_bio.tolist()
            savedict['comment'] = 'These are gridded parameter bio values to use for member %d' % mem.membernumber
            ncf.add_data(savedict)

            #--- All parameters, gridded ---#
            griddata = self.vector2grid(vectordata=data)
            griddata_fossil = griddata[:]
            w = np.where(fossil_land==0)
            griddata_fossil[w] = 1.0

            savedict = io.std_savedict.copy()
            savedict['name'] = "parametermap_ff"
            savedict['long_name'] = "parametermap_ff_for_member_%d" % mem.membernumber
            savedict['units'] = "unitless"
            savedict['dims'] = dimgrid 
            savedict['values'] = griddata_fossil.tolist()
            savedict['comment'] = 'These are gridded parameter fossil values to use for member %d' % mem.membernumber
            ncf.add_data(savedict)

            ncf.close()

            logging.debug('Successfully wrote data from ensemble member %d to file (%s) ' % (mem.membernumber, filename))

#    def write_to_file(self, filename, qual):
#        """
#        :param filename: the full filename for the output NetCDF file
#        :rtype: None
#
#        Write the StateVector information to a NetCDF file for later use. 
#        In principle the output file will have only one two datasets inside 
#        called:
#            * `meanstate`, dimensions [nlag, nparamaters]
#            * `ensemblestate`, dimensions [nlag,nmembers, nparameters]
#
#        This NetCDF information can be read back into a StateVector object using 
#        :meth:`~da.baseclasses.statevector.StateVector.read_from_file`
#
#        """
#        #import da.tools.io4 as io
#        #import da.tools.io as io
#
#        if qual == 'prior':
#            f = io.CT_CDF(filename, method='create')
#            logging.debug('Creating new StateVector output file (%s)' % filename)
#            #qual = 'prior'
#        else:
#            f = io.CT_CDF(filename, method='write')
#            logging.debug('Opening existing StateVector output file (%s)' % filename)
#            #qual = 'opt'
#
#        dimparams = f.add_params_dim(self.nparams)
#        dimmembers = f.add_members_dim(self.nmembers)
#        dimlag = f.add_lag_dim(self.nlag, unlimited=True)
#
#        for n in range(self.nlag):
#            members = self.ensemble_members[n]
#            mean_state = members[0].param_values
#
#            savedict = f.standard_var(varname='meanstate_%s' % qual)
#            savedict['dims'] = dimlag + dimparams 
#            savedict['values'] = mean_state
#            savedict['count'] = n
#            savedict['comment'] = 'this represents the mean of the ensemble'
#            f.add_data(savedict)
#
#            savedict = f.standard_var(varname='meanstate_bio_%s' % qual)
#            savedict['dims'] = dimlag + dimparams_bio
#            savedict['values'] = mean_state_bio
#            savedict['count'] = n
#            savedict['comment'] = 'this represents the mean of the ensemble'
#            f.add_data(savedict)
#
#            savedict = f.standard_var(varname='meanstate_fossil_%s' % qual)
#            savedict['dims'] = dimlag + dimparams_fossil
#            savedict['values'] = mean_state_fossil
#            savedict['count'] = n
#            savedict['comment'] = 'this represents the mean of the ensemble'
#            f.add_data(savedict)
#
#            members = self.ensemble_members[n]
#            devs = np.asarray([m.param_values.flatten() for m in members])
#            data = devs - np.asarray(mean_state)
#
#            savedict = f.standard_var(varname='ensemblestate_%s' % qual)
#            savedict['dims'] = dimlag + dimmembers + dimparams 
#            savedict['values'] = data
#            savedict['count'] = n
#            savedict['comment'] = 'this represents deviations from the mean of the ensemble'
#            f.add_data(savedict)
#        f.close()
#
#        logging.info('Successfully wrote the State Vector to file (%s) ' % filename)


    def make_new_ensemble(self, lag, dacycle, covariancematrix=None):
        """ 
        :param lag: an integer indicating the time step in the lag order
        :param covariancematrix: a matrix to draw random values from
        :rtype: None
    
        Make a new ensemble, the attribute lag refers to the position in the state vector. 
        Note that lag=1 means an index of 0 in python, hence the notation lag-1 in the indexing below.
        The argument is thus referring to the lagged state vector as [1,2,3,4,5,..., nlag]

        The optional covariance object to be passed holds a matrix of dimensions [nparams, nparams] which is
        used to draw ensemblemembers from. If this argument is not passed it will be substituted with an 
        identity matrix of the same dimensions.

        """    

        if covariancematrix == None: 
            covariancematrix = np.identity(self.nparams)

        # Make a cholesky decomposition of the covariance matrix


        _, s, _ = np.linalg.svd(covariancematrix)
        dof = np.sum(s) ** 2 / sum(s ** 2)
        C = np.linalg.cholesky(covariancematrix)

        logging.debug('Cholesky decomposition has succeeded ')
        logging.info('Appr. degrees of freedom in covariance matrix is %s' % (int(dof)))


        # Create mean values 

        newmean = np.ones(self.nparams, float) # standard value for a new time step is 1.0

        # If this is not the start of the filter, average previous two optimized steps into the mix

        if lag == self.nlag - 1 and self.nlag >= 3:
            newmean += self.ensemble_members[lag - 1][0].param_values + \
                                           self.ensemble_members[lag - 2][0].param_values 
            newmean = newmean / 3.0

        # Create the first ensemble member with a deviation of 0.0 and add to list

        newmember = EnsembleMember(0)
        newmember.param_values = newmean.flatten()  # no deviations
        self.ensemble_members[lag].append(newmember)

        # Create members 1:nmembers and add to ensemble_members list
        #dist = self.distribution

        for member in range(1, self.nmembers):
           # if dist == 'normal':
           #    rands = np.random.randn(self.nparams)
           # elif dist == 'gamma11':
           #    rands = np.random.exponential(1,size=self.nparams) - 1 #substruct mean=1.
           # elif dist == 'lognormal':
           #    rands = np.random.lognormal(0,0.5,size=self.nparams) - 1 #substruct median=1.
           # else:
           #    logging.error('Distribution (%s)to generate from not known' %dist)
           #    sys.exit(2)

            rands = np.random.randn(self.nparams)
            newmember = EnsembleMember(member)
            newmember.param_values = np.dot(C, rands) + newmean
            self.ensemble_members[lag].append(newmember)

        logging.info('%d new ensemble members generated from %s were added to the state vector # %d'\
                      #%(self.nmembers,dist, (lag + 1)))
                      %(self.nmembers,'normal', (lag + 1)))

    def propagate(self, dacycle):
        """ 
        :rtype: None

        Propagate the parameter values in the StateVector to the next cycle. This means a shift by one cycle 
        step for all states that will
        be optimized once more, and the creation of a new ensemble for the time step that just 
        comes in for the first time (step=nlag). 
        In the future, this routine can incorporate a formal propagation of the statevector.

        """
        
        # Remove State Vector n=1 by simply "popping" it from the list and appending a new empty list at the front. This empty list will
        # hold the new ensemble for the new cycle 

        self.ensemble_members.pop(0)
        self.ensemble_members.append([])

        # And now create a new time step of mean + members for n=nlag
        date = dacycle['time.start'] + timedelta(days=(self.nlag - 0.5) * int(dacycle['time.cycle']))
        cov = self.get_covariance(date, dacycle)
        self.make_new_ensemble(self.nlag - 1, dacycle, cov)

        logging.info('The state vector has been propagated by one cycle')

################### End Class MethaneStateVector ###################


if __name__ == "__main__":
    pass
