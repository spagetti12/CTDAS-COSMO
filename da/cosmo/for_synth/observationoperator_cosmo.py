#!/usr/bin/env python
# model.py

import logging, os, sys, subprocess, csv, cfgrib, os.path
import da.cosmo.io4 as io
import numpy as np
from netCDF4 import Dataset
from scipy import interpolate
from datetime import datetime, timedelta
from dateutil import rrule
from itertools import repeat
from multiprocessing import Pool
from da.tools.general import to_datetime
import amrs.misc.transform as transform

identifier = 'ObservationOperator'
version = '10'

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

    def run(self, lag, dacycle, statevector):
        members = statevector.ensemble_members[lag]
        self.forecast_nmembers = int(self.dacycle['da.optimizer.nmembers'])
        self.nparams = int(self.dacycle['nparameters'])
        self.days = int(dacycle['time.cycle'])
        abs_start_time = str((to_datetime(dacycle['abs.time.start'])).strftime('%Y%m%d%H'))
        abs_start_time_ch = str((to_datetime(dacycle['abs.time.start'])).strftime('%Y-%m-%d'))
        starth = abs((to_datetime(dacycle['abs.time.start'])-dacycle['time.start']).days)*24
        endh = abs((to_datetime(dacycle['abs.time.start'])-dacycle['time.finish']).days)*24

        start = dacycle['time.start']
        end = dacycle['time.finish']

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

        savedict = io.std_savedict.copy() 
        savedict['name'] = "synthetic"
        savedict['dtype'] = "float"
        savedict['long_name'] = "sythetic observations"
        savedict['units'] = "ppm"
        savedict['dims'] = dimid
        savedict['comment'] = "Simulated observations by COSMO"
        f.add_data(savedict,nsets=0)

        # Open file with x,y,z,t of model samples that need to be sampled

        f_in = io.ct_read(self.dacycle['ObsOperator.inputfile'],method='read')

        # Get simulated values and ID

        ids = f_in.get_variable('obs_num')
        obs = f_in.get_variable('observed')
        mdm = f_in.get_variable('modeldatamismatch')

        f_in.close()

        memplus1 = int(self.forecast_nmembers + 1)

        model_data = np.empty(shape=(self.forecast_nmembers, obs.size))   # +1 for OBS!!!!
        l_bg = np.empty(shape=memplus1)

        for m in range(0,int(self.forecast_nmembers)):
            l_bg[m] = members[m].param_values[-1]
        l_bg[-1] = 1.

       # LAMBDA FILE 1
        self.lambda_file1 = os.path.join(self.outputdir, 'lambda_1.%s.nc' % self.dacycle['time.sample.stamp'])
        ofile = Dataset(self.lambda_file1, mode='w')
        opar = ofile.createDimension('nparam', self.nparams)
        omem = ofile.createDimension('nensembles', memplus1/2)

        l = ofile.createVariable('lambda', np.float32, ('nensembles','nparam'),fill_value=-999.99)
        co2 = np.empty(shape=(int(memplus1/2),self.nparams))

        for m in range(0,int(memplus1/2)):
            co2[m,:] = members[m].param_values

        l[:] = co2
        ofile.close()

       # LAMBDA FILE 2
        self.lambda_file2 = os.path.join(self.outputdir, 'lambda_2.%s.nc' % self.dacycle['time.sample.stamp'])
        ofile = Dataset(self.lambda_file2, mode='w')
        opar = ofile.createDimension('nparam', self.nparams)
        omem = ofile.createDimension('nensembles', int(memplus1/2))

        l = ofile.createVariable('lambda', np.float32, ('nensembles','nparam'),fill_value=-999.99)
        co2 = np.empty(shape=(int(memplus1/2),self.nparams))

        for m in range(int(memplus1/2),int(memplus1-1)):
            co2[int(m-self.forecast_nmembers/2),:] = members[m].param_values

        # True lambdas!
        co2[-1,:] = np.hstack((0.7*np.ones(180),0.987))

        l[:] = co2
        ofile.close()

        os.system('cp '+self.lambda_file1+' '+dacycle['da.vprm']+'/lambdas_1.nc')
        os.system('cp '+self.lambda_file2+' '+dacycle['da.vprm']+'/lambdas_2.nc')

        os.chdir(dacycle['da.obsoperator.home'])

        logging.info('Starting COSMO')

        runnames = [self.dacycle['run.name1'], self.dacycle['run.name2']]
        args = [
            (name, abs_start_time_ch, str(starth+lag*168), str(endh+lag*168))
            for name in runnames
        ]

        with Pool(self.days) as pool:
            pool.starmap(self.start_cosmo, args)

        logging.info('COSMO done!')

        # Here the extraction of COSMO output starts
        dicts = self.read_csv(dacycle)
        rlat, rlon, dicts, path_in1, path_in2 = self.get_hhl_data(dacycle, lag, 'lffd'+abs_start_time+'c', dicts, starth, endh)
        logging.info('Starting parallel extraction \m/')

        args = [
            (dacycle, dacycle['time.sample.start']+timedelta(hours = 24*n), dicts, rlat, rlon, path_in1, path_in2, l_bg)
            for n in range(self.days)
        ]

        with Pool(self.days) as pool:
            pool.starmap(self.get_cosmo_data, args)

        logging.info('Finished parallel extraction \m/')

        self.cat_cosmo_data(dacycle)

        for i in range(memplus1):
            idx = str(i+1).zfill(3)
            cosmo_file = os.path.join(self.dacycle['dir.ct_save'], 'Hx_'+idx+'_%s.nc' % dacycle['time.sample.stamp'])
            ifile = Dataset(cosmo_file, mode='r')
            if i == memplus1-1:
                obs = np.squeeze(ifile.variables['CO2'][:])
            else:
                model_data[i,:] = np.squeeze(ifile.variables['CO2'][:])
            ifile.close()

        for j,data in enumerate(zip(ids,obs,mdm)):
            f.variables['obs_num'][j] = data[0]
            f.variables['flask'][j,:] = model_data[:,j]
            f.variables['synthetic'][j] = obs[j] + np.random.normal(loc=0, scale=np.sqrt(data[2]))
        f.close()

        logging.info('ObservationOperator finished successfully, output file written (%s)' % self.simulated_file)

    def run_forecast_model(self, lag, dacycle, statevector):
        self.prepare_run()
        self.run(lag, dacycle, statevector)

    def read_csv(self, dacycle):
        """Reads csv file where information about stations is written"""
        ensnum = int(dacycle['da.optimizer.nmembers'])+1
        csvfile = dacycle['locations']
        dicts = []
        with open(csvfile) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            for row in csv_reader:
                for e in range(ensnum):
                    e = str(e+1).zfill(3)
                    dicts = np.append(dicts, {'ensnum':e, 'name':row[1], 'lon':row[2], 'lat':row[3], \
                                              'rlon':None, 'rlat':None, 'h1':None, 'h2':None, \
                                              'hidx1':None, 'hidx2':None, \
                                              'alt':float(row[4])+float(row[5]), 'time':None, \
                                              'run':None, 'co2':[], \
                                              'co2_bg':[], 'co2_bio':[], 'co2_a':[]})
        return dicts

    def get_hhl_data(self, dacycle, lag, cfile, dicts, starth, endh):
        abs_start_time = str((to_datetime(dacycle['abs.time.start'])).strftime('%Y%m%d%H'))
        path_in1 = os.path.join(dacycle['dir.da_run1'], abs_start_time+'_'+str(starth+lag*168)+'_'+str(endh+lag*168), "cosmo/output/")
        path_in2 = os.path.join(dacycle['dir.da_run2'], abs_start_time+'_'+str(starth+lag*168)+'_'+str(endh+lag*168), "cosmo/output/")
        hhl = np.empty(shape=(60))
        hhl60 = np.empty(shape=(60, 150, 225))
        if not os.path.exists(path_in1+cfile):
            logging.info('Simulation file doesnt exist: %s' % path_in1+cfile)
            sys.exit(1)
        Filter = {'filter_by_keys':{'typeOfLevel':'hybrid'}}

        ds = cfgrib.open_dataset(path_in1+cfile, backend_kwargs=Filter)

        pollat = -1 * ds.HH.attrs['GRIB_latitudeOfSouthernPoleInDegrees']
        pollon = ds.HH.attrs['GRIB_longitudeOfSouthernPoleInDegrees'] - 180

        nx = ds.HH.attrs['GRIB_Nx']
        ny = ds.HH.attrs['GRIB_Ny']

        startlat = ds.HH.attrs['GRIB_latitudeOfFirstGridPointInDegrees']
        endlat = ds.HH.attrs['GRIB_latitudeOfLastGridPointInDegrees']

        startlon = ds.HH.attrs['GRIB_longitudeOfFirstGridPointInDegrees']
        endlon = ds.HH.attrs['GRIB_longitudeOfLastGridPointInDegrees']

        rlat = np.linspace(startlat,endlat,ny)
        rlon = np.linspace(startlon,endlon,nx)

        hhl_3d = ds.HH.data

        for h in range(0,60):
            hhl60[h, :, :]=(hhl_3d[h, :, :]+hhl_3d[h+1, :, :])/2.

        for station in dicts:
            myrlon, myrlat = transform.rotpole2wgs(float(station['lon']), \
                                                   float(station['lat']), \
                                                   pollon, pollat, inverse=True)
            for h in range(0,60):
                hhl[h] = interpolate.interpn((rlat,rlon), hhl60[h,:,:], [myrlat,myrlon], method='linear')

            if float(station['alt']) < hhl[59]:
                station['h1'] = hhl[59]
                station['h2'] = hhl[59]
                station['hidx1'] = 59
                station['hidx2'] = 59
                station['alt'] = hhl[59]

            else:
                for l, ll in enumerate(hhl):
                    if float(station['alt']) < ll:
                        station['h1'] = hhl[l]
                        station['h2'] = hhl[l+1]
                        station['hidx1'] = l
                        station['hidx2'] = l+1
                        # The following line: we interpolate on the middle of alt - lower level
            #            print(station['name'], station['alt'], float(station['alt']) - (float(station['alt']) - hhl[l+1])/2.)
                #        station['alt'] = float(station['alt']) - (float(station['alt']) - hhl[l+1])/2.
                station['alt'] = float(station['alt']) - (float(station['alt']) - float(station['h2']))/2.

            station['rlon'] = myrlon
            station['rlat'] = myrlat

        return rlat, rlon, dicts, path_in1, path_in2

    def get_cosmo_data(self, dacycle, date_begin, dicts, rlat, rlon, path_in1, path_in2, l_bg):
        hours = ['12', '13', '14', '15']
        qv_int = np.empty(shape=(60))
        co2_bg_int = np.empty(shape=(60))
        co2_bio_int = np.empty(shape=(60))
        co2_a_int = np.empty(shape=(60))
        qv_int = np.empty(shape=(60))

        co2_bg_hourly = []
        co2_bio_hourly = []
        co2_a_hourly = []
        co2_hourly = []

        for hrs in hours:
            fn1 = (path_in1+'lffd'+date_begin.strftime("%Y%m%d")+hrs)
            fn2 = (path_in2+'lffd'+date_begin.strftime("%Y%m%d")+hrs)

            ds = [cfgrib.open_datasets(fn1+'_1')[0], cfgrib.open_datasets(fn1+'_2')[0], cfgrib.open_datasets(fn2+'_1')[0], cfgrib.open_datasets(fn2+'_2')[0]]
#            logging.info('Reading of the files done for %s' % date_begin.strftime("%Y%m%d")+hrs)

            qv = ds[0]['q'].values
            co2_a = ds[0]['CO2_A'].values
            co2_bg001 = ds[0]['CO2_BG001'].values

            for station in dicts:
                myrlat = station['rlat']
                myrlon = station['rlon']
                e  = station['ensnum']
                h1 = station['h1']
                h2 = station['h2']
                i1 = station['hidx1']
                i2 = station['hidx2']
                if int(e) <= 6:
                    open_fn = 0
                    relative_e = e
                elif int(e) >= 7 and int(e) <= 21:
                    open_fn = 1
                    relative_e = e
                elif int(e) >= 22 and int(e) <= 27:
                    open_fn = 2
                    relative_e = str(int(e) - 21).zfill(3)
                elif int(e) >= 28:
                    open_fn = 3
                    relative_e = str(int(e) - 21).zfill(3)

                co2_bg  = co2_bg001 * l_bg[int(e)-1]
                co2_bio = ds[open_fn]['CO2_BIO'+relative_e].values

                for h in range(60):
                    qv_int[h] = interpolate.interpn((rlat,rlon), qv[h,:,:], [myrlat,myrlon], method='linear')
                    co2_bg_int[h] = interpolate.interpn((rlat,rlon), co2_bg[h,:,:], [myrlat,myrlon], method='linear')
                    co2_bio_int[h] = interpolate.interpn((rlat,rlon), co2_bio[h,:,:], [myrlat,myrlon], method='linear')
                    co2_a_int[h] = interpolate.interpn((rlat,rlon), co2_a[h,:,:], [myrlat,myrlon], method='linear')

                co2_bg1 = co2_bg_int[i1]
                co2_bg2 = co2_bg_int[i2]
                co2_bio1 = co2_bio_int[i1]
                co2_bio2 = co2_bio_int[i2]
                co2_a1 = co2_a_int[i1]
                co2_a2 = co2_a_int[i2]
                qv1 = qv_int[i1]
                qv2 = qv_int[i2]
                if h1 == h2:
                    co2_bg_final = co2_bg1
                    co2_bio_final = co2_bio1
                    co2_a_final = co2_a1
                    qv_final = qv1
                else:
                    co2_bg_final = co2_bg1 + (float(station['alt'])-h1)*(co2_bg2-co2_bg1)/(h2-h1)
                    co2_bio_final = co2_bio1 + (float(station['alt'])-h1)*(co2_bio2-co2_bio1)/(h2-h1)
                    co2_a_final = co2_a1 + (float(station['alt'])-h1)*(co2_a2-co2_a1)/(h2-h1)
                    qv_final = qv1 + (float(station['alt'])-h1)*(qv2-qv1)/(h2-h1)

                kgkg2ppm = 658941.149738696 / (1-qv_final)

                co2_bg_final = kgkg2ppm*co2_bg_final
                co2_bio_final = kgkg2ppm*co2_bio_final
                co2_a_final = kgkg2ppm*co2_a_final
                co2_final = co2_bio_final + co2_a_final

                station['co2_bg'].append(co2_bg_final)
                station['co2_bio'].append(co2_bio_final)
                station['co2_a'].append(co2_a_final)
                station['co2'].append(co2_final)

        for station in dicts:
            station['co2_bg'] = np.mean(np.asarray(station['co2_bg']))
            station['co2_bio'] = np.mean(np.asarray(station['co2_bio']))
            station['co2_a'] = np.mean(np.asarray(station['co2_a']))
            station['co2'] = np.mean(np.asarray(station['co2']))

            station['time'] = (date_begin - datetime(1970,1,1)).total_seconds()

        self.write_cosmo_data(dacycle, dicts, date_begin)

    def write_cosmo_data(self, dacycle, dicts, date):
        co2_all = []
        time_all = []
        date = date.strftime("%Y%m%d")
        for station in dicts:
            e = station['ensnum']
            if not os.path.exists(os.path.join(dacycle['dir.ct_save'])):
                os.mkdir(os.path.join(dacycle['dir.ct_save']))
            if not os.path.exists(os.path.join(dacycle['dir.ct_save'], "hourly")):
                os.mkdir(os.path.join(dacycle['dir.ct_save'], "hourly"))

            filename = os.path.join(dacycle['dir.ct_save'], 'hourly', station['name']+'_'+e+'_'+date+'.nc')
            with Dataset(filename, mode='w') as ofile:
                olev = ofile.createDimension('alt', 1)
                olat = ofile.createDimension('lat', 1)
                olon = ofile.createDimension('lon', 1)
                otime = ofile.createDimension('time', 1)

                olat = ofile.createVariable('lat', np.float64, ('lat',))
                olon = ofile.createVariable('lon', np.float64, ('lon',))
                olev = ofile.createVariable('alt', np.float64, ('alt',))
                otime = ofile.createVariable('time', int, ('time',))

                oco2_bg = ofile.createVariable('CO2_BG', np.float32, ('time','alt','lat','lon'),fill_value=-999.99)
                oco2_bio = ofile.createVariable('CO2_BIO', np.float32, ('time','alt','lat','lon'),fill_value=-999.99)
                oco2_a = ofile.createVariable('CO2_A', np.float32, ('time','alt','lat','lon'),fill_value=-999.99)
                oco2 = ofile.createVariable('CO2', np.float32, ('time','alt','lat','lon'),fill_value=-999.99)

                otime.units = 'seconds since 20190301 00:00:00'
                otime.calendar = 'proleptic_gregorian'

                olat[:] = float(station['lat'])
                olon[:] = float(station['lon'])
                olev[:] = float(station['alt'])
                otime[:] = (np.asarray(station['time'])).astype(int)

                oco2_bg[:] = np.asarray(station['co2_bg']).astype(float)
                oco2_bio[:] = np.asarray(station['co2_bio']).astype(float)
                oco2_a[:] = np.asarray(station['co2_a']).astype(float)
                oco2[:] = np.asarray(station['co2']).astype(float)
        logging.info('Written data for %s' % station['time'])


    def cat_cosmo_data(self, dacycle):

        date_begin = dacycle['time.sample.start']
        date_end = dacycle['time.sample.end']
        path_in = os.path.join(dacycle['dir.ct_save'], 'hourly')
        path_out = dacycle['dir.ct_save']
        csvfile = dacycle['locations']
        ensnum = int(dacycle['da.optimizer.nmembers'])+1
        for e in range(ensnum):
            co2_bg_all = []
            co2_bio_all = []
            co2_a_all = []
            co2_all = []
            e = str(e+1).zfill(3)
            with open(csvfile) as csv_file:
                csv_reader = csv.reader(csv_file, delimiter=',')
                for row in csv_reader:
                    for dt in rrule.rrule(rrule.DAILY, dtstart=date_begin, until=date_end):
                        date = dt.strftime("%Y%m%d")
                        with Dataset(os.path.join(path_in, row[1]+'_'+e+'_'+date+'.nc')) as nc:
                            co2_bg = np.squeeze(nc.variables['CO2_BG'][:])
                            co2_bio = np.squeeze(nc.variables['CO2_BIO'][:])
                            co2_a = np.squeeze(nc.variables['CO2_A'][:])
                            co2 = np.squeeze(nc.variables['CO2'][:])
                            co2_bg_all.append(co2_bg)
                            co2_bio_all.append(co2_bio)
                            co2_a_all.append(co2_a)
                            co2_all.append(co2)

            co2_bg_all = np.asarray(co2_bg_all)
            co2_bio_all = np.asarray(co2_bio_all)
            co2_a_all = np.asarray(co2_a_all)
            co2_all = np.asarray(co2_all)

            co2_bg_all = co2_bg_all.flatten()
            co2_bio_all = co2_bio_all.flatten()
            co2_a_all = co2_a_all.flatten()
            co2_all = co2_all.flatten()

            filename = os.path.join(path_out, 'Hx_'+e+'_'+dacycle['time.sample.stamp']+'.nc')
            filename_old = os.path.join(path_out, 'Hx_'+e+'_'+dacycle['time.sample.stamp']+'_once_opt.nc')
            if os.path.isfile(filename):
                os.system('mv '+filename+' '+filename_old)

            with Dataset(filename, mode='w') as ofile:
                otime = ofile.createDimension('obs', len(co2_all))
                otime = ofile.createVariable('obs', int, ('obs',))
                oco2_bg = ofile.createVariable('CO2_BG', np.float32, ('obs'),fill_value=-999.99)
                oco2_bio = ofile.createVariable('CO2_BIO', np.float32, ('obs'),fill_value=-999.99)
                oco2_a = ofile.createVariable('CO2_A', np.float32, ('obs'),fill_value=-999.99)
                oco2 = ofile.createVariable('CO2', np.float32, ('obs'),fill_value=-999.99)

                otime[:] = np.linspace(1, len(co2_all), len(co2_all))
                oco2_bg[:] = co2_bg_all
                oco2_bio[:] = co2_bio_all
                oco2_a[:] = co2_a_all
                oco2[:] = co2_all

    def start_cosmo(self, runname, abs_start_time_ch, starthour, endhour):
        os.system('python run_chain.py '+runname+' '+abs_start_time_ch+' '+starthour+' '+endhour+' -j meteo icbc int2lm post_int2lm oae octe online_vprm cosmo -f')



################### End Class ObservationOperator ###################

class RandomizerObservationOperator(ObservationOperator):
    """ This class holds methods and variables that are needed to use a random number generated as substitute
        for a true observation operator. It takes observations and returns values for each obs, with a specified 
        amount of white noise added 
    """



if __name__ == "__main__":
    pass
