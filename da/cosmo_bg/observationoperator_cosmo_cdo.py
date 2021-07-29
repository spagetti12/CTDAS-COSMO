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
from itertools import repeat
from multiprocessing import Pool
from da.tools.general import to_datetime

identifier = 'ObservationOperator'
version = '10'

#cdo = Cdo()
cdo = Cdo(logging=True, logFile='cdo_commands.log')

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

        model_data = np.empty(shape=(self.forecast_nmembers,obs.size))   # 3x7

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
            model_data[i,:] = (np.squeeze(ifile.variables['CO2'][:])*29./44.01)*1E6   # in ppm DRY
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

        files2cat_albs = []
        files2cat_bntg = []
        files2cat_brm = []
        files2cat_chri = []
        files2cat_due1 = []
        files2cat_esmo = []
        files2cat_frob = []
        files2cat_gimm = []
        files2cat_hae = []
        files2cat_laeg = []
        files2cat_magn = []
        files2cat_payn = []
        files2cat_reck = []
        files2cat_rig = []
        files2cat_save = []
        files2cat_semp = []
        files2cat_sott = []
        files2cat_ssal = []
        files2cat_taen = []
        files2cat_zhbr = []
        files2cat_zsch = []
        files2cat_zue = []

        if ens == "001":
            cdo.selname("HHL", input = hhl_fn, output = cosmo_out+"hhl.nc")
            cdo.remapnn("lon=8.51_lat=47.31,", input = cosmo_out+"hhl.nc", output = cosmo_out+"hhl_albs.nc")
            cdo.remapnn("lon=7.53_lat=46.98,", input = cosmo_out+"hhl.nc", output = cosmo_out+"hhl_bntg.nc")
            cdo.remapnn("lon=8.18_lat=47.19,", input = cosmo_out+"hhl.nc", output = cosmo_out+"hhl_brm.nc")
            cdo.remapnn("lon=7.69_lat=47.57,", input = cosmo_out+"hhl.nc", output = cosmo_out+"hhl_chri.nc")
            cdo.remapnn("lon=8.61_lat=47.40,", input = cosmo_out+"hhl.nc", output = cosmo_out+"hhl_due1.nc")
            cdo.remapnn("lon=8.57_lat=47.52,", input = cosmo_out+"hhl.nc", output = cosmo_out+"hhl_esmo.nc")
            cdo.remapnn("lon=7.90_lat=47.38,", input = cosmo_out+"hhl.nc", output = cosmo_out+"hhl_frob.nc")
            cdo.remapnn("lon=7.25_lat=47.05,", input = cosmo_out+"hhl.nc", output = cosmo_out+"hhl_gimm.nc")
            cdo.remapnn("lon=7.82_lat=47.31,", input = cosmo_out+"hhl.nc", output = cosmo_out+"hhl_hae.nc")
            cdo.remapnn("lon=8.40_lat=47.48,", input = cosmo_out+"hhl.nc", output = cosmo_out+"hhl_laeg.nc")
            cdo.remapnn("lon=8.93_lat=46.16,", input = cosmo_out+"hhl.nc", output = cosmo_out+"hhl_magn.nc")
            cdo.remapnn("lon=6.94_lat=46.81,", input = cosmo_out+"hhl.nc", output = cosmo_out+"hhl_payn.nc")
            cdo.remapnn("lon=8.52_lat=47.43,", input = cosmo_out+"hhl.nc", output = cosmo_out+"hhl_reck.nc")
            cdo.remapnn("lon=8.46_lat=47.07,", input = cosmo_out+"hhl.nc", output = cosmo_out+"hhl_rig.nc")
            cdo.remapnn("lon=7.36_lat=46.24,", input = cosmo_out+"hhl.nc", output = cosmo_out+"hhl_save.nc")
            cdo.remapnn("lon=8.21_lat=47.12,", input = cosmo_out+"hhl.nc", output = cosmo_out+"hhl_semp.nc")
            cdo.remapnn("lon=6.74_lat=46.66,", input = cosmo_out+"hhl.nc", output = cosmo_out+"hhl_sott.nc")
            cdo.remapnn("lon=8.95_lat=45.98,", input = cosmo_out+"hhl.nc", output = cosmo_out+"hhl_ssal.nc")
            cdo.remapnn("lon=8.90_lat=47.48,", input = cosmo_out+"hhl.nc", output = cosmo_out+"hhl_taen.nc")
            cdo.remapnn("lon=8.57_lat=47.38,", input = cosmo_out+"hhl.nc", output = cosmo_out+"hhl_zhbr.nc")
            cdo.remapnn("lon=8.52_lat=47.37,", input = cosmo_out+"hhl.nc", output = cosmo_out+"hhl_zsch.nc")
            cdo.remapnn("lon=8.53_lat=47.38,", input = cosmo_out+"hhl.nc", output = cosmo_out+"hhl_zue.nc")

        for dt in rrule.rrule(rrule.HOURLY, dtstart=to_datetime(dacycle['abs.time.start'])+timedelta(hours=hstart), until=to_datetime(dacycle['abs.time.start'])+timedelta(hours=hstop)):
            dt=dt.strftime('%Y%m%d%H')
            if ens == "001":
                logging.info('Extracting output for time %s' % (str(dt)))
            co2_in_fn = cosmo_out+'lffd'+dt+'.nc'

            co2_out_albs = cosmo_out+'CO2_albs_'+ens+'_'+dt+'.nc'
            co2_out_bntg = cosmo_out+'CO2_bntg_'+ens+'_'+dt+'.nc'
            co2_out_brm = cosmo_out+'CO2_brm_'+ens+'_'+dt+'.nc'
            co2_out_chri = cosmo_out+'CO2_chri_'+ens+'_'+dt+'.nc'
            co2_out_due1 = cosmo_out+'CO2_due1_'+ens+'_'+dt+'.nc'
            co2_out_esmo = cosmo_out+'CO2_esmo_'+ens+'_'+dt+'.nc'
            co2_out_frob = cosmo_out+'CO2_frob_'+ens+'_'+dt+'.nc'
            co2_out_gimm = cosmo_out+'CO2_gimm_'+ens+'_'+dt+'.nc'
            co2_out_hae = cosmo_out+'CO2_hae_'+ens+'_'+dt+'.nc'
            co2_out_laeg = cosmo_out+'CO2_laeg_'+ens+'_'+dt+'.nc'
            co2_out_magn = cosmo_out+'CO2_magn_'+ens+'_'+dt+'.nc'
            co2_out_payn = cosmo_out+'CO2_payn_'+ens+'_'+dt+'.nc'
            co2_out_reck = cosmo_out+'CO2_reck_'+ens+'_'+dt+'.nc'
            co2_out_rig = cosmo_out+'CO2_rig_'+ens+'_'+dt+'.nc'
            co2_out_save = cosmo_out+'CO2_save_'+ens+'_'+dt+'.nc'
            co2_out_semp = cosmo_out+'CO2_semp_'+ens+'_'+dt+'.nc'
            co2_out_sott = cosmo_out+'CO2_sott_'+ens+'_'+dt+'.nc'
            co2_out_ssal = cosmo_out+'CO2_ssal_'+ens+'_'+dt+'.nc'
            co2_out_taen = cosmo_out+'CO2_taen_'+ens+'_'+dt+'.nc'
            co2_out_zhbr = cosmo_out+'CO2_zhbr_'+ens+'_'+dt+'.nc'
            co2_out_zsch = cosmo_out+'CO2_zsch_'+ens+'_'+dt+'.nc'
            co2_out_zue = cosmo_out+'CO2_zue_'+ens+'_'+dt+'.nc'

            cdo.expr("'CO2=(CO2_BG"+ens+"-CO2_GPP"+ens+"+CO2_RA"+ens+"+CO2_A)/(1.-QV)'", input = "-remapnn,lon=8.51_lat=47.31 -selname,QV,CO2_BG"+ens+",CO2_GPP"+ens+",CO2_RA"+ens+",CO2_A "+co2_in_fn, output = co2_out_albs)
            cdo.expr("'CO2=(CO2_BG"+ens+"-CO2_GPP"+ens+"+CO2_RA"+ens+"+CO2_A)/(1.-QV)'", input = "-remapnn,lon=7.53_lat=46.98 -selname,QV,CO2_BG"+ens+",CO2_GPP"+ens+",CO2_RA"+ens+",CO2_A "+co2_in_fn, output = co2_out_bntg)
            cdo.expr("'CO2=(CO2_BG"+ens+"-CO2_GPP"+ens+"+CO2_RA"+ens+"+CO2_A)/(1.-QV)'", input = "-remapnn,lon=8.18_lat=47.19 -selname,QV,CO2_BG"+ens+",CO2_GPP"+ens+",CO2_RA"+ens+",CO2_A "+co2_in_fn, output = co2_out_brm)
            cdo.expr("'CO2=(CO2_BG"+ens+"-CO2_GPP"+ens+"+CO2_RA"+ens+"+CO2_A)/(1.-QV)'", input = "-remapnn,lon=7.69_lat=47.57 -selname,QV,CO2_BG"+ens+",CO2_GPP"+ens+",CO2_RA"+ens+",CO2_A "+co2_in_fn, output = co2_out_chri)
            cdo.expr("'CO2=(CO2_BG"+ens+"-CO2_GPP"+ens+"+CO2_RA"+ens+"+CO2_A)/(1.-QV)'", input = "-remapnn,lon=8.61_lat=47.40 -selname,QV,CO2_BG"+ens+",CO2_GPP"+ens+",CO2_RA"+ens+",CO2_A "+co2_in_fn, output = co2_out_due1)
            cdo.expr("'CO2=(CO2_BG"+ens+"-CO2_GPP"+ens+"+CO2_RA"+ens+"+CO2_A)/(1.-QV)'", input = "-remapnn,lon=8.57_lat=47.52 -selname,QV,CO2_BG"+ens+",CO2_GPP"+ens+",CO2_RA"+ens+",CO2_A "+co2_in_fn, output = co2_out_esmo)
            cdo.expr("'CO2=(CO2_BG"+ens+"-CO2_GPP"+ens+"+CO2_RA"+ens+"+CO2_A)/(1.-QV)'", input = "-remapnn,lon=7.90_lat=47.38 -selname,QV,CO2_BG"+ens+",CO2_GPP"+ens+",CO2_RA"+ens+",CO2_A "+co2_in_fn, output = co2_out_frob)
            cdo.expr("'CO2=(CO2_BG"+ens+"-CO2_GPP"+ens+"+CO2_RA"+ens+"+CO2_A)/(1.-QV)'", input = "-remapnn,lon=7.25_lat=47.05 -selname,QV,CO2_BG"+ens+",CO2_GPP"+ens+",CO2_RA"+ens+",CO2_A "+co2_in_fn, output = co2_out_gimm)
            cdo.expr("'CO2=(CO2_BG"+ens+"-CO2_GPP"+ens+"+CO2_RA"+ens+"+CO2_A)/(1.-QV)'", input = "-remapnn,lon=7.82_lat=47.31 -selname,QV,CO2_BG"+ens+",CO2_GPP"+ens+",CO2_RA"+ens+",CO2_A "+co2_in_fn, output = co2_out_hae)
            cdo.expr("'CO2=(CO2_BG"+ens+"-CO2_GPP"+ens+"+CO2_RA"+ens+"+CO2_A)/(1.-QV)'", input = "-remapnn,lon=8.40_lat=47.48 -selname,QV,CO2_BG"+ens+",CO2_GPP"+ens+",CO2_RA"+ens+",CO2_A "+co2_in_fn, output = co2_out_laeg)
            cdo.expr("'CO2=(CO2_BG"+ens+"-CO2_GPP"+ens+"+CO2_RA"+ens+"+CO2_A)/(1.-QV)'", input = "-remapnn,lon=8.93_lat=46.16 -selname,QV,CO2_BG"+ens+",CO2_GPP"+ens+",CO2_RA"+ens+",CO2_A "+co2_in_fn, output = co2_out_magn)
            cdo.expr("'CO2=(CO2_BG"+ens+"-CO2_GPP"+ens+"+CO2_RA"+ens+"+CO2_A)/(1.-QV)'", input = "-remapnn,lon=6.94_lat=46.81 -selname,QV,CO2_BG"+ens+",CO2_GPP"+ens+",CO2_RA"+ens+",CO2_A "+co2_in_fn, output = co2_out_payn)
            cdo.expr("'CO2=(CO2_BG"+ens+"-CO2_GPP"+ens+"+CO2_RA"+ens+"+CO2_A)/(1.-QV)'", input = "-remapnn,lon=8.52_lat=47.43 -selname,QV,CO2_BG"+ens+",CO2_GPP"+ens+",CO2_RA"+ens+",CO2_A "+co2_in_fn, output = co2_out_reck)
            cdo.expr("'CO2=(CO2_BG"+ens+"-CO2_GPP"+ens+"+CO2_RA"+ens+"+CO2_A)/(1.-QV)'", input = "-remapnn,lon=8.46_lat=47.07 -selname,QV,CO2_BG"+ens+",CO2_GPP"+ens+",CO2_RA"+ens+",CO2_A "+co2_in_fn, output = co2_out_rig)
            cdo.expr("'CO2=(CO2_BG"+ens+"-CO2_GPP"+ens+"+CO2_RA"+ens+"+CO2_A)/(1.-QV)'", input = "-remapnn,lon=7.36_lat=46.24 -selname,QV,CO2_BG"+ens+",CO2_GPP"+ens+",CO2_RA"+ens+",CO2_A "+co2_in_fn, output = co2_out_save)
            cdo.expr("'CO2=(CO2_BG"+ens+"-CO2_GPP"+ens+"+CO2_RA"+ens+"+CO2_A)/(1.-QV)'", input = "-remapnn,lon=8.21_lat=47.12 -selname,QV,CO2_BG"+ens+",CO2_GPP"+ens+",CO2_RA"+ens+",CO2_A "+co2_in_fn, output = co2_out_semp)
            cdo.expr("'CO2=(CO2_BG"+ens+"-CO2_GPP"+ens+"+CO2_RA"+ens+"+CO2_A)/(1.-QV)'", input = "-remapnn,lon=6.74_lat=46.66 -selname,QV,CO2_BG"+ens+",CO2_GPP"+ens+",CO2_RA"+ens+",CO2_A "+co2_in_fn, output = co2_out_sott)
            cdo.expr("'CO2=(CO2_BG"+ens+"-CO2_GPP"+ens+"+CO2_RA"+ens+"+CO2_A)/(1.-QV)'", input = "-remapnn,lon=8.95_lat=45.98 -selname,QV,CO2_BG"+ens+",CO2_GPP"+ens+",CO2_RA"+ens+",CO2_A "+co2_in_fn, output = co2_out_ssal)
            cdo.expr("'CO2=(CO2_BG"+ens+"-CO2_GPP"+ens+"+CO2_RA"+ens+"+CO2_A)/(1.-QV)'", input = "-remapnn,lon=8.90_lat=47.48 -selname,QV,CO2_BG"+ens+",CO2_GPP"+ens+",CO2_RA"+ens+",CO2_A "+co2_in_fn, output = co2_out_taen)
            cdo.expr("'CO2=(CO2_BG"+ens+"-CO2_GPP"+ens+"+CO2_RA"+ens+"+CO2_A)/(1.-QV)'", input = "-remapnn,lon=8.57_lat=47.38 -selname,QV,CO2_BG"+ens+",CO2_GPP"+ens+",CO2_RA"+ens+",CO2_A "+co2_in_fn, output = co2_out_zhbr)
            cdo.expr("'CO2=(CO2_BG"+ens+"-CO2_GPP"+ens+"+CO2_RA"+ens+"+CO2_A)/(1.-QV)'", input = "-remapnn,lon=8.52_lat=47.37 -selname,QV,CO2_BG"+ens+",CO2_GPP"+ens+",CO2_RA"+ens+",CO2_A "+co2_in_fn, output = co2_out_zsch)
            cdo.expr("'CO2=(CO2_BG"+ens+"-CO2_GPP"+ens+"+CO2_RA"+ens+"+CO2_A)/(1.-QV)'", input = "-remapnn,lon=8.53_lat=47.38 -selname,QV,CO2_BG"+ens+",CO2_GPP"+ens+",CO2_RA"+ens+",CO2_A "+co2_in_fn, output = co2_out_zue)

            files2cat_albs.append(co2_out_albs)
            files2cat_bntg.append(co2_out_bntg)
            files2cat_brm.append(co2_out_brm)
            files2cat_chri.append(co2_out_chri)
            files2cat_due1.append(co2_out_due1)
            files2cat_esmo.append(co2_out_esmo)
            files2cat_frob.append(co2_out_frob)
            files2cat_gimm.append(co2_out_gimm)
            files2cat_hae.append(co2_out_hae)
            files2cat_laeg.append(co2_out_laeg)
            files2cat_magn.append(co2_out_magn)
            files2cat_payn.append(co2_out_payn)
            files2cat_reck.append(co2_out_reck)
            files2cat_rig.append(co2_out_rig)
            files2cat_save.append(co2_out_save)
            files2cat_semp.append(co2_out_semp)
            files2cat_sott.append(co2_out_sott)
            files2cat_ssal.append(co2_out_ssal)
            files2cat_taen.append(co2_out_taen)
            files2cat_zhbr.append(co2_out_zhbr)
            files2cat_zsch.append(co2_out_zsch)
            files2cat_zue.append(co2_out_zue)

        cdo.cat(input = files2cat_albs, output = cosmo_out+"CO2_albs_"+ens+"_"+time_stamp+".nc")
        cdo.cat(input = files2cat_bntg, output = cosmo_out+"CO2_bntg_"+ens+"_"+time_stamp+".nc")
        cdo.cat(input = files2cat_brm, output = cosmo_out+"CO2_brm_"+ens+"_"+time_stamp+".nc")
        cdo.cat(input = files2cat_chri, output = cosmo_out+"CO2_chri_"+ens+"_"+time_stamp+".nc")
        cdo.cat(input = files2cat_due1, output = cosmo_out+"CO2_due1_"+ens+"_"+time_stamp+".nc")
        cdo.cat(input = files2cat_esmo, output = cosmo_out+"CO2_esmo_"+ens+"_"+time_stamp+".nc")
        cdo.cat(input = files2cat_frob, output = cosmo_out+"CO2_frob_"+ens+"_"+time_stamp+".nc")
        cdo.cat(input = files2cat_gimm, output = cosmo_out+"CO2_gimm_"+ens+"_"+time_stamp+".nc")
        cdo.cat(input = files2cat_hae, output = cosmo_out+"CO2_hae_"+ens+"_"+time_stamp+".nc")
        cdo.cat(input = files2cat_laeg, output = cosmo_out+"CO2_laeg_"+ens+"_"+time_stamp+".nc")
        cdo.cat(input = files2cat_magn, output = cosmo_out+"CO2_magn_"+ens+"_"+time_stamp+".nc")
        cdo.cat(input = files2cat_payn, output = cosmo_out+"CO2_payn_"+ens+"_"+time_stamp+".nc")
        cdo.cat(input = files2cat_reck, output = cosmo_out+"CO2_reck_"+ens+"_"+time_stamp+".nc")
        cdo.cat(input = files2cat_rig, output = cosmo_out+"CO2_rig_"+ens+"_"+time_stamp+".nc")
        cdo.cat(input = files2cat_save, output = cosmo_out+"CO2_save_"+ens+"_"+time_stamp+".nc")
        cdo.cat(input = files2cat_semp, output = cosmo_out+"CO2_semp_"+ens+"_"+time_stamp+".nc")
        cdo.cat(input = files2cat_sott, output = cosmo_out+"CO2_sott_"+ens+"_"+time_stamp+".nc")
        cdo.cat(input = files2cat_ssal, output = cosmo_out+"CO2_ssal_"+ens+"_"+time_stamp+".nc")
        cdo.cat(input = files2cat_taen, output = cosmo_out+"CO2_taen_"+ens+"_"+time_stamp+".nc")
        cdo.cat(input = files2cat_zhbr, output = cosmo_out+"CO2_zhbr_"+ens+"_"+time_stamp+".nc")
        cdo.cat(input = files2cat_zsch, output = cosmo_out+"CO2_zsch_"+ens+"_"+time_stamp+".nc")
        cdo.cat(input = files2cat_zue, output = cosmo_out+"CO2_zue_"+ens+"_"+time_stamp+".nc")


        sites = ("albs", "bntg", "brm", "chri", "due1", "esmo", "frob", "gimm", "hae", "laeg", "magn", "payn", "reck", "rig", "save", "semp", "sott", "ssal", "taen", "zhbr", "zsch", "zue")
        for s,ss in enumerate(sites):
            site_height.main(cosmo_out, str(ens), ss, time_stamp)

        cdo.intlevel("857.25", input = cosmo_out+"CO2_60lev_"+ens+"_albs_"+time_stamp+".nc", output = cosmo_out+"modelled_"+ens+"_albs_"+time_stamp+".nc")
        cdo.intlevel("995.35", input = cosmo_out+"CO2_60lev_"+ens+"_bntg_"+time_stamp+".nc", output = cosmo_out+"modelled_"+ens+"_bntg_"+time_stamp+".nc")
        cdo.intlevel("1009.7", input = cosmo_out+"CO2_60lev_"+ens+"_brm_"+time_stamp+".nc", output = cosmo_out+"modelled_"+ens+"_brm_"+time_stamp+".nc")
        cdo.intlevel("592.8", input = cosmo_out+"CO2_60lev_"+ens+"_chri_"+time_stamp+".nc", output = cosmo_out+"modelled_"+ens+"_chri_"+time_stamp+".nc")
        #cdo.intlevel("432.4", input = cosmo_out+"CO2_60lev_"+ens+"_due1_"+time_stamp+".nc", output = cosmo_out+"modelled_"+ens+"_due1_"+time_stamp+".nc")
        cdo.sellevidx("60", input = cosmo_out+"CO2_60lev_"+ens+"_due1_"+time_stamp+".nc", output = cosmo_out+"modelled_"+ens+"_due1_"+time_stamp+".nc")
        cdo.intlevel("560.9", input = cosmo_out+"CO2_60lev_"+ens+"_esmo_"+time_stamp+".nc", output = cosmo_out+"modelled_"+ens+"_esmo_"+time_stamp+".nc")
        cdo.intlevel("897.3", input = cosmo_out+"CO2_60lev_"+ens+"_frob_"+time_stamp+".nc", output = cosmo_out+"modelled_"+ens+"_frob_"+time_stamp+".nc")
        cdo.intlevel("478.5", input = cosmo_out+"CO2_60lev_"+ens+"_gimm_"+time_stamp+".nc", output = cosmo_out+"modelled_"+ens+"_gimm_"+time_stamp+".nc")
        #cdo.intlevel("430.2", input = cosmo_out+"CO2_60lev_"+ens+"_hae_"+time_stamp+".nc", output = cosmo_out+"modelled_"+ens+"_hae_"+time_stamp+".nc")
        cdo.sellevidx("60", input = cosmo_out+"CO2_60lev_"+ens+"_hae_"+time_stamp+".nc", output = cosmo_out+"modelled_"+ens+"_hae_"+time_stamp+".nc")
        cdo.intlevel("855", input = cosmo_out+"CO2_60lev_"+ens+"_laeg_"+time_stamp+".nc", output = cosmo_out+"modelled_"+ens+"_laeg_"+time_stamp+".nc")
        #cdo.intlevel("208", input = cosmo_out+"CO2_60lev_"+ens+"_magn_"+time_stamp+".nc", output = cosmo_out+"modelled_"+ens+"_magn_"+time_stamp+".nc")
        cdo.sellevidx("60", input = cosmo_out+"CO2_60lev_"+ens+"_magn_"+time_stamp+".nc", output = cosmo_out+"modelled_"+ens+"_magn_"+time_stamp+".nc")
        cdo.intlevel("488.7", input = cosmo_out+"CO2_60lev_"+ens+"_payn_"+time_stamp+".nc", output = cosmo_out+"modelled_"+ens+"_payn_"+time_stamp+".nc")
        #cdo.intlevel("443.1", input = cosmo_out+"CO2_60lev_"+ens+"_reck_"+time_stamp+".nc", output = cosmo_out+"modelled_"+ens+"_reck_"+time_stamp+".nc")
        cdo.sellevidx("60", input = cosmo_out+"CO2_60lev_"+ens+"_reck_"+time_stamp+".nc", output = cosmo_out+"modelled_"+ens+"_reck_"+time_stamp+".nc")
        cdo.intlevel("1030.5", input = cosmo_out+"CO2_60lev_"+ens+"_rig_"+time_stamp+".nc", output = cosmo_out+"modelled_"+ens+"_rig_"+time_stamp+".nc")
        cdo.intlevel("789.4", input = cosmo_out+"CO2_60lev_"+ens+"_save_"+time_stamp+".nc", output = cosmo_out+"modelled_"+ens+"_save_"+time_stamp+".nc")
        cdo.intlevel("583.3", input = cosmo_out+"CO2_60lev_"+ens+"_semp_"+time_stamp+".nc", output = cosmo_out+"modelled_"+ens+"_semp_"+time_stamp+".nc")
        cdo.intlevel("774.9", input = cosmo_out+"CO2_60lev_"+ens+"_sott_"+time_stamp+".nc", output = cosmo_out+"modelled_"+ens+"_sott_"+time_stamp+".nc")
        cdo.intlevel("876.1", input = cosmo_out+"CO2_60lev_"+ens+"_ssal_"+time_stamp+".nc", output = cosmo_out+"modelled_"+ens+"_ssal_"+time_stamp+".nc")
        #cdo.intlevel("543.1", input = cosmo_out+"CO2_60lev_"+ens+"_taen_"+time_stamp+".nc", output = cosmo_out+"modelled_"+ens+"_taen_"+time_stamp+".nc")
        cdo.sellevidx("60", input = cosmo_out+"CO2_60lev_"+ens+"_taen_"+time_stamp+".nc", output = cosmo_out+"modelled_"+ens+"_taen_"+time_stamp+".nc")
        cdo.intlevel("615.8", input = cosmo_out+"CO2_60lev_"+ens+"_zhbr_"+time_stamp+".nc", output = cosmo_out+"modelled_"+ens+"_zhbr_"+time_stamp+".nc")
        #cdo.intlevel("417.1", input = cosmo_out+"CO2_60lev_"+ens+"_zsch_"+time_stamp+".nc", output = cosmo_out+"modelled_"+ens+"_zsch_"+time_stamp+".nc")
        cdo.sellevidx("60", input = cosmo_out+"CO2_60lev_"+ens+"_zsch_"+time_stamp+".nc", output = cosmo_out+"modelled_"+ens+"_zsch_"+time_stamp+".nc")
        #cdo.intlevel("408.8", input = cosmo_out+"CO2_60lev_"+ens+"_zue_"+time_stamp+".nc", output = cosmo_out+"modelled_"+ens+"_zue_"+time_stamp+".nc")
        cdo.sellevidx("60", input = cosmo_out+"CO2_60lev_"+ens+"_zue_"+time_stamp+".nc", output = cosmo_out+"modelled_"+ens+"_zue_"+time_stamp+".nc")

        cdo.cat(input = cosmo_out+"modelled_"+ens+"_*_"+time_stamp+".nc ", output = cosmo_save+"model_"+ens+"_"+time_stamp+".nc")
        logging.info('Extracting done for ens %s' % (ens))



################### End Class ObservationOperator ###################

class RandomizerObservationOperator(ObservationOperator):
    """ This class holds methods and variables that are needed to use a random number generated as substitute
        for a true observation operator. It takes observations and returns values for each obs, with a specified 
        amount of white noise added 
    """



if __name__ == "__main__":
    pass
