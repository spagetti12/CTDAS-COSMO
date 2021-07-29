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
# obs.py

"""
Author : peters 

Revision History:
File created on 28 Jul 2010.

"""
import os
import sys
import logging
        
import datetime as dtm
#from string import strip
from numpy import array, logical_and, sqrt
sys.path.append(os.getcwd())
sys.path.append('../../')

identifier = 'CarbonTracker CO2 mole fractions'
version = '0.0'

from da.baseclasses.obs import Observations
import da.tools.io4 as io
import da.tools.rc as rc
################### Begin Class ObsPackObservations ###################

class ObsPackObservations(Observations):
    """ an object that holds data + methods and attributes needed to manipulate mole fraction values """

    def setup(self, dacycle):

        self.startdate = dacycle['time.sample.start']
        self.enddate = dacycle['time.sample.end']

        op_id = dacycle.dasystem['obspack.input.id']
        op_dir = dacycle.dasystem['obspack.input.dir']

        if not os.path.exists(op_dir):
            msg = 'Could not find  the required ObsPack distribution (%s) ' % op_dir
            logging.error(msg)
            raise IOError(msg)
        else:
            self.obspack_dir = op_dir
            self.obspack_id = op_id

        self.datalist = []

    def add_observations(self):
        """ Returns a MoleFractionList holding individual MoleFractionSample objects for all obs in a file
      
            The ObsPack mole fraction files are provided as time series per site with all dates in sequence. 
            We will loop over all site files in the ObsPackage, and subset each to our needs
            
        """

        # Step 1: Read list of available site files in package

        infile = os.path.join(self.obspack_dir, 'list.txt')
        f = open(infile, 'r')
        lines = f.readlines()
        f.close()

        ncfilelist = []
        for line in lines:
            if not line.startswith('# dataset:'): continue

            items = line.split(':')
            #ncfile, lab , start_date, stop_date, data_comparison = items[0:5]
            #ncfile, lab , start_date, stop_date, data_comparison= line[:105].split()
            ncfile = items[1].strip()

            ncfilelist += [ncfile]

        logging.debug("ObsPack dataset info read, proceeding with %d netcdf files" % len(ncfilelist))

        for ncfile in ncfilelist:

            infile = os.path.join(self.obspack_dir, ncfile + '.nc')
            ncf = io.ct_read(infile, 'read')
            idates = ncf.get_variable('time_components')
            dates = array([dtm.datetime(*d) for d in idates])

            subselect = logical_and(dates >= self.startdate , dates <= self.enddate).nonzero()[0]

            dates = dates.take(subselect, axis=0)
 
            if 'merge_num' in ncf.variables:
                obspacknum = ncf.get_variable('merge_num').take(subselect)
            else:
                obspacknum = ncf.get_variable('obspack_num').take(subselect)
            if 'ccggAllData' in ncfile:
                obspackid = ncf.get_variable('id').take(subselect, axis=0)
            else:
                obspackid = ncf.get_variable('obspack_id').take(subselect, axis=0)
            obspackid = [s.tostring().lower() for s in obspackid]
            obspackid = list(map(str.strip,str(obspackid)))
            datasetname = ncfile  # use full name of dataset to propagate for clarity
            lats = ncf.get_variable('lat').take(subselect, axis=0)
            lons = ncf.get_variable('lon').take(subselect, axis=0)
            alts = ncf.get_variable('altitude').take(subselect, axis=0)
            obs = ncf.get_variable('value').take(subselect, axis=0)
            species = 'co2' #ncf.get_attribute('dataset_parameter')
            flags = ncf.get_variable('obs_flag').take(subselect, axis=0)
            ncf.close()

            for n in range(len(dates)): 
                self.datalist.append(MoleFractionSample(obspacknum[n], dates[n], datasetname, obs[n], 0.0, 0.0, 0.0, 0.0, flags[n], alts[n], lats[n], lons[n], obspackid[n], species, 1, 0.0, infile))

            logging.debug("Added %d observations from file (%s) to the Data list" % (len(dates), ncfile)) 

        logging.info("Observations list now holds %d values" % len(self.datalist))

    def add_simulations(self, filename, silent=False):
        """ Adds model simulated values to the mole fraction objects """


        if not os.path.exists(filename):
            msg = "Sample output filename for observations could not be found : %s" % filename 
            logging.error(msg)
            logging.error("Did the sampling step succeed?")
            logging.error("...exiting")
            raise IOError(msg)

        ncf = io.ct_read(filename, method='read')
        ids = ncf.get_variable('obs_num')
        simulated = ncf.get_variable('flask')
        ncf.close()
        logging.info("Successfully read data from model sample file (%s)" % filename)

        obs_ids = self.getvalues('id').tolist()
        ids = list(map(int, ids))

        missing_samples = []

        for idx, val in zip(ids, simulated): 
            if idx in obs_ids:
                index = obs_ids.index(idx)

                self.datalist[index].simulated = val  # in mol/mol
            else:     
                missing_samples.append(idx)

        if not silent and missing_samples != []:
            logging.warning('Model samples were found that did not match any ID in the observation list. Skipping them...')
            #msg = '%s'%missing_samples ; logging.warning(msg)

        logging.debug("Added %d simulated values to the Data list" % (len(ids) - len(missing_samples)))

    def write_sample_coords(self, obsinputfile):
        """ 
            Write the information needed by the observation operator to a file. Return the filename that was written for later use

        """

        if len(self.datalist) == 0:
            #f.close()
            #return obsinputfile
            logging.debug("No observations found for this time period, nothing written to obs file")
        else:
            f = io.CT_CDF(obsinputfile, method='create')
            logging.debug('Creating new observations file for ObservationOperator (%s)' % obsinputfile)

            dimid = f.add_dim('obs', len(self.datalist))
            dim200char = f.add_dim('string_of200chars', 200)
            dim10char = f.add_dim('string_of10chars', 10)
            dimcalcomp = f.add_dim('calendar_components', 6)

            for key, value in self.site_move.items():
                msg = "Site is moved by %3.2f degrees latitude and %3.2f degrees longitude" % value 
                f.add_attribute(key, msg)

            data = self.getvalues('id')

            savedict = io.std_savedict.copy() 
            savedict['name'] = "obs_num"
            savedict['dtype'] = "int"
            savedict['long_name'] = "Unique_Dataset_observation_index_number"
            savedict['units'] = ""
            savedict['dims'] = dimid
            savedict['values'] = data.tolist()
            savedict['comment'] = "Unique index number within this dataset ranging from 0 to UNLIMITED."
            f.add_data(savedict)

            data = [[d.year, d.month, d.day, d.hour, d.minute, d.second] for d in self.getvalues('xdate') ]

            savedict = io.std_savedict.copy() 
            savedict['dtype'] = "int"
            savedict['name'] = "date_components"
            savedict['units'] = "integer components of UTC date/time"
            savedict['dims'] = dimid + dimcalcomp
            savedict['values'] = data
            savedict['missing_value'] = -9
            savedict['comment'] = "Calendar date components as integers. Times and dates are UTC." 
            savedict['order'] = "year, month, day, hour, minute, second"
            f.add_data(savedict)

            data = self.getvalues('lat')

            savedict = io.std_savedict.copy() 
            savedict['name'] = "latitude"
            savedict['units'] = "degrees_north"
            savedict['dims'] = dimid
            savedict['values'] = data.tolist()
            savedict['missing_value'] = -999.9
            f.add_data(savedict)

            data = self.getvalues('lon')

            savedict = io.std_savedict.copy() 
            savedict['name'] = "longitude"
            savedict['units'] = "degrees_east"
            savedict['dims'] = dimid
            savedict['values'] = data.tolist()
            savedict['missing_value'] = -999.9
            f.add_data(savedict)

            data = self.getvalues('height')

            savedict = io.std_savedict.copy() 
            savedict['name'] = "altitude"
            savedict['units'] = "meters_above_sea_level"
            savedict['dims'] = dimid
            savedict['values'] = data.tolist()
            savedict['missing_value'] = -999.9
            f.add_data(savedict)

            data = self.getvalues('samplingstrategy')

            savedict = io.std_savedict.copy() 
            savedict['dtype'] = "int"
            savedict['name'] = "sampling_strategy"
            savedict['units'] = "NA"
            savedict['dims'] = dimid
            savedict['values'] = data.tolist()
            savedict['missing_value'] = -9
            f.add_data(savedict)

            data = self.getvalues('evn')

            savedict = io.std_savedict.copy() 
            savedict['dtype'] = "char"
            savedict['name'] = "obs_id"
            savedict['units'] = "ObsPack datapoint identifier"
            savedict['dims'] = dimid + dim200char
            savedict['values'] = data
            savedict['missing_value'] = '!'
            f.add_data(savedict)

            data = self.getvalues('obs')

            savedict = io.std_savedict.copy()
            savedict['name'] = "observed"
            savedict['long_name'] = "observedvalues"
            savedict['units'] = "mol mol-1"
            savedict['dims'] = dimid
            savedict['values'] = data.tolist()
            savedict['comment'] = 'Observations used in optimization'
            f.add_data(savedict)
    
            data = self.getvalues('mdm')
    
            savedict = io.std_savedict.copy()
            savedict['name'] = "modeldatamismatch"
            savedict['long_name'] = "modeldatamismatch"
            savedict['units'] = "[mol mol-1]"
            savedict['dims'] = dimid
            savedict['values'] = data.tolist()
            savedict['comment'] = 'Standard deviation of mole fractions resulting from model-data mismatch'
            f.add_data(savedict)
            f.close()

            logging.debug("Successfully wrote data to obs file")
            logging.info("Sample input file for obs operator now in place [%s]" % obsinputfile)

        

    def add_model_data_mismatch(self, filename):
        """ 
            Get the model-data mismatch values for this cycle.

                (1) Open a sites_weights file
                (2) Parse the data
                (3) Compare site list against data
                (4) Take care of double sites, etc

        """    

        if not os.path.exists(filename):
            msg = 'Could not find  the required sites.rc input file (%s) ' % filename
            logging.error(msg)
            raise IOError(msg)
        else:
            self.sites_file = filename

        sites_weights = rc.read(self.sites_file)

        self.rejection_threshold = int(sites_weights['obs.rejection.threshold'])
        self.global_R_scaling = float(sites_weights['global.R.scaling'])
        self.n_site_categories = int(sites_weights['n.site.categories'])

        logging.debug('Model-data mismatch rejection threshold: %d ' % self.rejection_threshold)
        logging.warning('Model-data mismatch scaling factor     : %f ' % self.global_R_scaling)
        logging.debug('Model-data mismatch site categories    : %d ' % self.n_site_categories)
   
        cats = [k for k in sites_weights.keys() if 'site.category' in k] 

        site_categories = {}
        for key in cats:
            name, error, may_localize, may_reject = sites_weights[key].split(';')
            name = name.strip().lower()
            error = float(error)
            may_reject = ("TRUE" in may_reject.upper())
            may_localize = ("TRUE" in may_localize.upper())
            site_categories[name] = {'category': name, 'error': error, 'may_localize': may_localize, 'may_reject': may_reject}
            print('VAZNO', site_categories[name])

        site_info = {}
        site_move = {}
        site_incalt = {} # option to increase sampling altitude for sites specified in sites and weights file 
        for key, value in sites_weights.items():
#            if 'co2_' in key or 'sf6' in key:  # to be fixed later, do not yet know how to parse valid keys from rc-files yet.... WP
            print(key, value)
            sitename, sitecategory = key, value
            sitename = sitename.strip()
            print(sitename)
            sitecategory = sitecategory.split()[0].strip().lower()
            print(sitecategory)
            print(site_categories[sitecategory])
            site_info[sitename] = site_categories[sitecategory]
            if 'site.move' in key:
                identifier, latmove, lonmove = value.split(';')
                site_move[identifier.strip()] = (float(latmove), float(lonmove))
            if 'site.incalt' in key:
                identifier, incalt = value.split(';')
                site_incalt[identifier.strip()] = (int(incalt))

        for obs in self.datalist:  # first loop over all available data points to set flags correctly

            obs.mdm = 1000.0  # default is very high model-data-mismatch, until explicitly set by script
            if obs.flag == 1: # flag is taken from the gv+ datasets: 1=background/representative, 0=local. 
                obs.flag = 0
            elif obs.flag == 0:
                obs.flag = 99 # 99 means: do-not-use
            else: obs.flag = 99    
                    
        for obs in self.datalist:  # second loop over all available data points to set mdm

            identifier = obs.code
            site = identifier.split('_')
            #species, site, method, lab, datasetnr = identifier.split('_')

            print('IDEN', identifier)
            print('SITE_I', site_info)

            if identifier in site_info:
                if site_info[identifier]['category'] == 'do-not-use' or obs.flag == 99:
                    logging.warning("Observation found (%s, %d), but not used in assimilation." % (identifier, obs.id))
                    obs.mdm = site_info[identifier]['error'] * self.global_R_scaling
                    obs.may_localize = site_info[identifier]['may_localize']
                    obs.may_reject = site_info[identifier]['may_reject']
                    obs.flag = 99
                else:
                    if site_info[identifier]['category'] == 'aircraft':
                        nr_obs_per_day = 1
                    else:
                        nr_obs_per_day = len([c.code for c in self.datalist if c.code == obs.code and c.xdate.day == obs.xdate.day and c.flag == 0])
                    logging.debug("Observation found (%s, %d), mdm category is: %0.2f, scaled with number of observations per day (%i), final mdm applied is: %0.2f." % (identifier, obs.id, site_info[identifier]['error'],nr_obs_per_day,site_info[identifier]['error']*sqrt(nr_obs_per_day)))
                    obs.mdm = site_info[identifier]['error'] * sqrt(nr_obs_per_day) * self.global_R_scaling
                    obs.may_localize = site_info[identifier]['may_localize']
                    obs.may_reject = site_info[identifier]['may_reject']
                    obs.flag = 0
            else:
                logging.warning("Observation NOT found (%s, %d), please check sites.rc file (%s)  !!!" % (identifier, obs.id, self.sites_file))

            if identifier in site_move:

                movelat, movelon = site_move[identifier]
                obs.lat = obs.lat + movelat
                obs.lon = obs.lon + movelon

                logging.warning("Observation location for (%s, %d), is moved by %3.2f degrees latitude and %3.2f degrees longitude" % (identifier, obs.id, movelat, movelon))

            if identifier in site_incalt:

                incalt = site_incalt[identifier]
                obs.height = obs.height + incalt

                logging.warning("Observation location for (%s, %d), is moved by %3.2f meters in altitude" % (identifier, obs.id, incalt))


        # Add site_info dictionary to the Observations object for future use

        self.site_info = site_info
        self.site_move = site_move
        self.site_incalt = site_incalt

        logging.debug("Added Model Data Mismatch to all samples ")

    def write_sample_auxiliary(self, auxoutputfile):
        """ 
            Write selected information contained in the Observations object to a file. 

        """

        f = io.CT_CDF(auxoutputfile, method='create')
        logging.debug('Creating new auxiliary sample output file for postprocessing (%s)' % auxoutputfile)

        dimid = f.add_dim('obs', len(self.datalist))
        dim200char = f.add_dim('string_of200chars', 200)
        dim10char = f.add_dim('string_of10chars', 10)
        dimcalcomp = f.add_dim('calendar_components', 6)

        if len(self.datalist) == 0:
            f.close()
            #return outfile

        for key, value in self.site_move.items():
            msg = "Site is moved by %3.2f degrees latitude and %3.2f degrees longitude" % value 
            f.add_attribute(key, msg)

        data = self.getvalues('id')

        savedict = io.std_savedict.copy() 
        savedict['name'] = "obs_num"
        savedict['dtype'] = "int"
        savedict['long_name'] = "Unique_Dataset_observation_index_number"
        savedict['units'] = ""
        savedict['dims'] = dimid
        savedict['values'] = data.tolist()
        savedict['comment'] = "Unique index number within this dataset ranging from 0 to UNLIMITED."
        f.add_data(savedict)

        data = [[d.year, d.month, d.day, d.hour, d.minute, d.second] for d in self.getvalues('xdate')]

        savedict = io.std_savedict.copy() 
        savedict['dtype'] = "int"
        savedict['name'] = "date_components"
        savedict['units'] = "integer components of UTC date/time"
        savedict['dims'] = dimid + dimcalcomp
        savedict['values'] = data
        savedict['missing_value'] = -9
        savedict['comment'] = "Calendar date components as integers. Times and dates are UTC." 
        savedict['order'] = "year, month, day, hour, minute, second"
        f.add_data(savedict)

        data = self.getvalues('obs')

        savedict = io.std_savedict.copy()
        savedict['name'] = "observed"
        savedict['long_name'] = "observedvalues"
        savedict['units'] = "mol mol-1"
        savedict['dims'] = dimid
        savedict['values'] = data.tolist()
        savedict['comment'] = 'Observations used in optimization'
        f.add_data(savedict)

        data = self.getvalues('mdm')

        savedict = io.std_savedict.copy()
        savedict['name'] = "modeldatamismatch"
        savedict['long_name'] = "modeldatamismatch"
        savedict['units'] = "[mol mol-1]"
        savedict['dims'] = dimid
        savedict['values'] = data.tolist()
        savedict['comment'] = 'Standard deviation of mole fractions resulting from model-data mismatch'
        f.add_data(savedict)

        data = self.getvalues('simulated') 

        dimmembers = f.add_dim('members', data.shape[1])

        savedict = io.std_savedict.copy()
        savedict['name'] = "modelsamples"
        savedict['long_name'] = "modelsamples for all ensemble members"
        savedict['units'] = "mol mol-1"
        savedict['dims'] = dimid + dimmembers
        savedict['values'] = data.tolist()
        savedict['comment'] = 'simulated mole fractions based on optimized state vector'
        f.add_data(savedict)

        data = self.getvalues('fromfile') 

        savedict = io.std_savedict.copy()
        savedict['name'] = "inputfilename"
        savedict['long_name'] = "name of file where original obs data was taken from"
        savedict['dtype'] = "char"
        savedict['dims'] = dimid + dim200char
        savedict['values'] = data
        savedict['missing_value'] = '!'
        f.add_data(savedict)

        f.close()

        logging.debug("Successfully wrote data to auxiliary sample output file (%s)" % auxoutputfile)

        #return outfile



################### End Class CtObservations ###################



################### Begin Class MoleFractionSample ###################

class MoleFractionSample(object):
    """ 
        Holds the data that defines a mole fraction Sample in the data assimilation framework. Sor far, this includes all
        attributes listed below in the __init__ method. One can additionally make more types of data, or make new
        objects for specific projects.

    """

    def __init__(self, idx, xdate, code='XXX', obs=0.0, simulated=0.0, resid=0.0, hphr=0.0, mdm=0.0, flag=0, height=0.0, lat= -999., lon= -999., evn='0000', species='co2', samplingstrategy=1, sdev=0.0, fromfile='none.nc'):
        self.code = code.strip()      # dataset identifier, i.e., co2_lef_tower_insitu_1_99
        self.xdate = xdate             # Date of obs
        self.obs = obs               # Value observed
        self.simulated = simulated         # Value simulated by model
        self.resid = resid             # Mole fraction residuals
        self.hphr = hphr              # Mole fraction prior uncertainty from fluxes and (HPH) and model data mismatch (R)
        self.mdm = mdm               # Model data mismatch
        self.may_localize = True           # Whether sample may be localized in optimizer
        self.may_reject = True              # Whether sample may be rejected if outside threshold
        self.flag = flag              # Flag
        self.height = height            # Sample height in masl
        self.lat = lat               # Sample lat
        self.lon = lon               # Sample lon
        self.id = idx               # Obspack ID within distrution (integer), e.g., 82536
        self.evn = evn               # Obspack Number within distrution (string), e.g., obspack_co2_1_PROTOTYPE_v0.9.2_2012-07-26_99_82536
        self.sdev = sdev              # standard deviation of ensemble
        self.masl = True              # Sample is in Meters Above Sea Level
        self.mag = not self.masl     # Sample is in Meters Above Ground
        self.species = species.strip()
        self.samplingstrategy = samplingstrategy
        self.fromfile = fromfile   # netcdf filename inside ObsPack distribution, to write back later

################### End Class MoleFractionSample ###################


if __name__ == "__main__":
    pass



