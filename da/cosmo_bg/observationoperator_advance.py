#!/usr/bin/env python
# model.py

import logging
import os
import sys
import subprocess
import da.cosmo.io4 as io
import numpy as np
from netCDF4 import Dataset
from datetime import datetime, timedelta
from dateutil import rrule
from cdo import *
from . import site_height
from da.cosmo.icbc4ctdas import ct
from itertools import repeat
from multiprocessing import Pool
from da.tools.general import to_datetime

identifier = 'ObservationOperator'
version = '10'

cdo = Cdo()

################### Begin Class ObservationOperator ###################
class ObservationOperator(object):
    def __init__(self, dacycle=None):
        self.ID = identifier
        self.version = version
        self.restart_filelist = []
        self.output_filelist = []
        self.outputdir = None # Needed for opening the samples.nc files created 

        logging.info('Observation Operator object initialized: %s' % self.ID)

        if dacycle != None:
            self.dacycle = dacycle
        else:
            self.dacycle = {}

    
    def get_initial_data(self):
        """ This method places all initial data needed by an ObservationOperator in the proper folder for the model """

    def setup(self,dacycle):
        """ Perform all steps necessary to start the observation operator through a simple Run() call """

        self.dacycle = dacycle
        self.outputdir = dacycle['dir.output']

    def prepare_run(self):
        """ Prepare the running of the actual forecast model, for example compile code """

        # Define the name of the file that will contain the modeled output of each observation

        self.simulated_file = os.path.join(self.outputdir, 'samples_simulated.%s.nc' % self.dacycle['time.sample.stamp'])
        self.forecast_nmembers = int(self.dacycle['da.optimizer.nmembers'])

    def run(self,lag,dacycle,statevector,advance=False):
        members = statevector.ensemble_members[lag]
        self.forecast_nmembers = int(self.dacycle['da.optimizer.nmembers'])
        self.nparams = int(self.dacycle['nparameters'])
        absolute_start_time = str((to_datetime(dacycle['abs.time.start'])).strftime('%Y%m%d%H'))
        absolute_start_time_ch = str((to_datetime(dacycle['abs.time.start'])).strftime('%Y-%m-%d'))
        starth = abs((to_datetime(dacycle['abs.time.start'])-dacycle['time.start']).days)*24
        endh = abs((to_datetime(dacycle['abs.time.start'])-dacycle['time.finish']).days)*24

        f = io.CT_CDF(self.simulated_file, method='create')
        logging.debug('Creating new simulated observation file in ObservationOperator (%s)' % self.simulated_file)

        dimid = f.createDimension('obs_num', size=None)
        dimid = ('obs_num',)
        savedict = io.std_savedict.copy() 
        savedict['name'] = "obs_num"
        savedict['dtype'] = "int"
        savedict['long_name'] = "Unique_Dataset_observation_index_number"
        savedict['units'] = ""
        savedict['dims'] = dimid
        savedict['comment'] = "Unique index number within this dataset ranging from 0 to UNLIMITED."
        f.add_data(savedict,nsets=0)

        dimmember = f.createDimension('nmembers', size=self.forecast_nmembers)
        dimmember = ('nmembers',)
        savedict = io.std_savedict.copy() 
        savedict['name'] = "flask"
        savedict['dtype'] = "float"
        savedict['long_name'] = "mole_fraction_of_trace_gas_in_dry_air"
        savedict['units'] = "ppm"
        savedict['dims'] = dimid + dimmember
        savedict['comment'] = "Simulated model value created by COSMO"

        f.add_data(savedict,nsets=0)

        # Open file with x,y,z,t of model samples that need to be sampled

        f_in = io.ct_read(self.dacycle['ObsOperator.inputfile'],method='read')

        # Get simulated values and ID

        ids = f_in.get_variable('obs_num')
        obs = f_in.get_variable('observed')
        mdm = f_in.get_variable('modeldatamismatch')

        f_in.close()

        model_data = np.empty(shape=(self.forecast_nmembers,mdm.size))   # 3x7

        self.lambda_file = os.path.join(self.outputdir, 'lambda.%s.nc' % self.dacycle['time.sample.stamp'])
        ofile = Dataset(self.lambda_file, mode='w')
        opar = ofile.createDimension('nparam', self.nparams)
        omem = ofile.createDimension('nensembles', self.forecast_nmembers)#len(members.nmembers))

        l = ofile.createVariable('lambda', np.float32, ('nensembles','nparam'),fill_value=-999.99)
        co2 = np.empty(shape=(self.forecast_nmembers,self.nparams))

        for m in range(0,self.forecast_nmembers):
            co2[m,:] = members[m].param_values
#        co2[co2<0] = 0.
        l[:] = co2
        ofile.close()
        os.system('cp '+self.lambda_file+' '+dacycle['da.vprm']+'/lambdas.nc')
        if advance:
            os.chdir(dacycle['da.obsoperator.home'])

            if os.path.exists(dacycle['dir.da_run']+'/'+absolute_start_time+"_"+str(starth+lag*168)+"_"+str(endh+lag*168)+"/cosmo/output/"):
                if os.path.exists(dacycle['dir.da_run']+"/non_opt_"+absolute_start_time+"_"+str(starth+lag*168)+"_"+str(endh+lag*168)+"/cosmo/output/"):
                    os.rename(dacycle['dir.da_run']+"/"+absolute_start_time+"_"+str(starth+lag*168)+"_"+str(endh+lag*168), dacycle['dir.da_run']+"/old_non_opt_"+dacycle['time.start'].strftime('%Y%m%d%H')+"_"+str(starth+lag*168)+"_"+str(endh+lag*168))
                else:
                    os.rename(dacycle['dir.da_run']+"/"+absolute_start_time+"_"+str(starth+lag*168)+"_"+str(endh+lag*168), dacycle['dir.da_run']+"/non_opt_"+dacycle['time.start'].strftime('%Y%m%d%H')+"_"+str(starth+lag*168)+"_"+str(endh+lag*168))

            os.system('python run_chain.py '+self.dacycle['run.name']+' '+absolute_start_time_ch+' '+str(starth+lag*168)+' '+str(endh+lag*168)+' -j meteo icbc int2lm post_int2lm oae octe online_vprm cosmo')
            logging.info('COSMO done!')
            os.chdir(dacycle['dir.da_run'])

            args = [
                (dacycle, starth+168*lag, endh+168*lag-1, n)
                for n in range(1,self.forecast_nmembers+1)
            ]

            with Pool(self.forecast_nmembers) as pool:
                pool.starmap(self.extract_model_data, args)

        for i in range(0,self.forecast_nmembers):
            idx = str(i+1).zfill(3)
            cosmo_file = os.path.join('/store/empa/em05/parsenov/cosmo_data/model_'+idx+'_%s.nc' % dacycle['time.sample.stamp'])
            ifile = Dataset(cosmo_file, mode='r')
            model_data[i,:] = (np.squeeze(ifile.variables['CO2'][:])*29./44.01)*1E6   # in ppm
            ifile.close()

        for j,data in enumerate(zip(ids,obs,mdm)):
            f.variables['obs_num'][j] = data[0]
            f.variables['flask'][j,:] = model_data[:,j]
        f.close()

        logging.info('ObservationOperator finished successfully, output file written (%s)' % self.simulated_file)

    def run_forecast_model(self, lag, dacycle, statevector, advance):
        self.prepare_run()
        self.run(lag, dacycle, statevector, advance)

    def extract_model_data(self, dacycle, hstart, hstop, ensnum):

        self.dacycle = dacycle
        time_stamp = dacycle['time.sample.stamp']
        abs_start_time = str((to_datetime(dacycle['abs.time.start'])).strftime('%Y%m%d%H'))

        cosmo_out = dacycle['dir.da_run']+"/"+abs_start_time+"_"+str(hstart)+"_"+str(hstop+1)+"/cosmo/output/"
        hhl_cosmo_out = dacycle['dir.da_run']+"/"+abs_start_time+"_0_168/cosmo/output/"
        cosmo_save = "/store/empa/em05/parsenov/cosmo_data/"
        hhl_fn = hhl_cosmo_out+'lffd'+abs_start_time+'c.nc'

        ens = str(ensnum).zfill(3)
        files2cat_jfj=[]
        files2cat_lhw=[]
        files2cat_brm=[]
        files2cat_ssl=[]

        if ens == "001":
            cdo.selname("HHL", input = hhl_fn, output = cosmo_out+"hhl.nc")
            cdo.remapnn("lon=7.99_lat=46.54,", input = cosmo_out+"hhl.nc", output = cosmo_out+"hhl_jfj.nc")
            cdo.remapnn("lon=8.40_lat=47.48,", input = cosmo_out+"hhl.nc", output = cosmo_out+"hhl_lhw.nc")
            cdo.remapnn("lon=8.18_lat=47.19,", input = cosmo_out+"hhl.nc", output = cosmo_out+"hhl_brm.nc")
            cdo.remapnn("lon=7.92_lat=47.92,", input = cosmo_out+"hhl.nc", output = cosmo_out+"hhl_ssl.nc")

        for dt in rrule.rrule(rrule.HOURLY, dtstart=to_datetime(dacycle['abs.time.start'])+timedelta(hours=hstart), until=to_datetime(dacycle['abs.time.start'])+timedelta(hours=hstop)):
            dt=dt.strftime('%Y%m%d%H')
            if ens == "001":
                logging.info('Extracting output for time %s' % (str(dt)))
            co2_in_fn = cosmo_out+'lffd'+dt+'.nc'
            co2_out_jfj = cosmo_out+'CO2_jfj_'+ens+'_'+dt+'.nc'
            co2_out_lhw = cosmo_out+'CO2_lhw_'+ens+'_'+dt+'.nc'
            co2_out_brm = cosmo_out+'CO2_brm_'+ens+'_'+dt+'.nc'
            co2_out_ssl = cosmo_out+'CO2_ssl_'+ens+'_'+dt+'.nc'
            cdo.expr("'CO2=(CO2_BG"+ens+"-CO2_GPP"+ens+"+CO2_RA"+ens+"+CO2_A)/(1.-QV)'", input = "-remapnn,lon=7.99_lat=46.54 -selname,QV,CO2_BG"+ens+",CO2_GPP"+ens+",CO2_RA"+ens+",CO2_A "+co2_in_fn, output = co2_out_jfj)
            cdo.expr("'CO2=(CO2_BG"+ens+"-CO2_GPP"+ens+"+CO2_RA"+ens+"+CO2_A)/(1.-QV)'", input = "-remapnn,lon=8.40_lat=47.48 -selname,QV,CO2_BG"+ens+",CO2_GPP"+ens+",CO2_RA"+ens+",CO2_A "+co2_in_fn, output = co2_out_lhw)
            cdo.expr("'CO2=(CO2_BG"+ens+"-CO2_GPP"+ens+"+CO2_RA"+ens+"+CO2_A)/(1.-QV)'", input = "-remapnn,lon=8.18_lat=47.19 -selname,QV,CO2_BG"+ens+",CO2_GPP"+ens+",CO2_RA"+ens+",CO2_A "+co2_in_fn, output = co2_out_brm)
            cdo.expr("'CO2=(CO2_BG"+ens+"-CO2_GPP"+ens+"+CO2_RA"+ens+"+CO2_A)/(1.-QV)'", input = "-remapnn,lon=7.92_lat=47.92 -selname,QV,CO2_BG"+ens+",CO2_GPP"+ens+",CO2_RA"+ens+",CO2_A "+co2_in_fn, output = co2_out_ssl)
            files2cat_jfj.append(co2_out_jfj)
            files2cat_lhw.append(co2_out_lhw)
            files2cat_brm.append(co2_out_brm)
            files2cat_ssl.append(co2_out_ssl)

        cdo.cat(input = files2cat_jfj, output = cosmo_out+"CO2_jfj_"+ens+"_"+time_stamp+".nc")
        cdo.cat(input = files2cat_lhw, output = cosmo_out+"CO2_lhw_"+ens+"_"+time_stamp+".nc")
        cdo.cat(input = files2cat_brm, output = cosmo_out+"CO2_brm_"+ens+"_"+time_stamp+".nc")
        cdo.cat(input = files2cat_ssl, output = cosmo_out+"CO2_ssl_"+ens+"_"+time_stamp+".nc")

        sites = ("lhw","brm","jfj","ssl")
        for s,ss in enumerate(sites):
            site_height.main(cosmo_out, str(ens), ss, time_stamp)

        cdo.intlevel("860", input = cosmo_out+"CO2_60lev_"+ens+"_lhw_"+time_stamp+".nc", output = cosmo_out+"modelled_"+ens+"_lhw_"+time_stamp+".nc")
        cdo.intlevel("1009", input = cosmo_out+"CO2_60lev_"+ens+"_brm_"+time_stamp+".nc", output = cosmo_out+"modelled_"+ens+"_brm_"+time_stamp+".nc")
        cdo.intlevel("3580", input = cosmo_out+"CO2_60lev_"+ens+"_jfj_"+time_stamp+".nc", output = cosmo_out+"modelled_"+ens+"_jfj_"+time_stamp+".nc")
        cdo.intlevel("1205", input = cosmo_out+"CO2_60lev_"+ens+"_ssl_"+time_stamp+".nc", output = cosmo_out+"modelled_"+ens+"_ssl_"+time_stamp+".nc")

        cdo.cat(input = cosmo_out+"modelled_"+ens+"_brm_"+time_stamp+".nc "+cosmo_out+"modelled_"+ens+"_jfj_"+time_stamp+".nc "+cosmo_out+"modelled_"+ens+"_lhw_"+time_stamp+".nc "+cosmo_out+"modelled_"+ens+"_ssl_"+time_stamp+".nc ", output = cosmo_save+"model_"+ens+"_"+time_stamp+".nc")
        logging.info('Extracting done for ens %s' % (ens))



################### End Class ObservationOperator ###################

class RandomizerObservationOperator(ObservationOperator):
    """ This class holds methods and variables that are needed to use a random number generated as substitute
        for a true observation operator. It takes observations and returns values for each obs, with a specified 
        amount of white noise added 
    """



if __name__ == "__main__":
    pass
