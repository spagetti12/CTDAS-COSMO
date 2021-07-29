import sys
import random
import numpy as np
from netCDF4 import Dataset

def main(dir, ens, station, times):

    ifile_co2 = Dataset(dir+'CO2_'+station+'_'+ens+'_'+times+'.nc', mode='r')
    ifile_hhl = Dataset(dir+'hhl_'+station+'.nc', mode='r')
    rlat = ifile_co2.variables['lat'][:]
    rlon = ifile_co2.variables['lon'][:]
    time = ifile_co2.variables['time'][:]
    co2 = ifile_co2.variables['CO2'][:]   #*29./44. # not needed, it's in the observationoperator class
    hhl = ifile_hhl.variables['HHL'][:]
    ifile_co2.close()
    ifile_hhl.close()

    hhl60 = np.empty(shape=60)
    for h in range(0,60):
        hhl60[h]=(hhl[0,h,0,0]+hhl[0,h+1,0,0])/2.

    co2_e = np.empty(shape=168)
    ofile = Dataset(dir+'CO2_60lev_'+ens+'_'+station+'_'+times+'.nc', mode='w')

    orlat = ofile.createDimension('rlat', 1)
    orlon = ofile.createDimension('rlon', 1)
    oz = ofile.createDimension('h', len(hhl60))
    otime = ofile.createDimension('time', len(time))

    orlat = ofile.createVariable('rlat', np.float64, ('rlat',))
    orlon = ofile.createVariable('rlon', np.float64, ('rlon',))
    oz = ofile.createVariable('h', np.float64, ('h',))
    otime = ofile.createVariable('time', np.float64, ('time',))

    oco2 = ofile.createVariable('CO2', np.float32, ('time','rlat','rlon','h'),fill_value=-999.99)

    oz.units = 'm'
    oz.axis = 'Z'

    orlat[:] = rlat
    orlon[:] = rlon
    oz[:] = hhl60
    otime[:] = time
    co2_e = np.empty(shape=(len(time),60,1,1))
    for t in range(len(time)):
        co2_e[t,...] = co2[t,...]
    oco2[:] = co2_e[:]
    ofile.close()
