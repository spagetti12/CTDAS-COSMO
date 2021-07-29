import sys
import numpy as np
from netCDF4 import Dataset

idir = '/store/empa/em05/parsenov/cosmo_input/icbc/'
odir = '/store/empa/em05/parsenov/cosmo_input/icbc/processed/ensemble/'

half_a=(0,7.367743,210.39389,855.361755,2063.779785,3850.91333,6144.314941,8802.356445,11632.75879,14411.12402,16899.46875,18864.75,20097.40234,20429.86328,19755.10938,18045.18359,15379.80566,12077.44629,8765.053711,6018.019531,3960.291504,1680.640259,713.218079,298.495789,95.636963,0)
half_b=(1,0.99401945,0.97966272,0.95182151,0.90788388,0.84737492,0.77159661,0.68326861,0.58616841,0.48477158,0.38389215,0.28832296,0.2024759,0.13002251,0.07353383,0.03412116,0.01114291,0.00181516,0.00007582,0,0,0,0,0,0,0)

def ct(date, co2_bg):
    date = str(date)
    co2_bg = np.array(co2_bg)
    ifile = Dataset(idir+'3d_molefractions_1x1_'+date+'.nc', mode='r')
    levs = ifile.variables['levels'][:]
    lats = ifile.variables['latitude'][:]
    lons = ifile.variables['longitude'][:]

    co2 = np.squeeze(ifile.variables['co2'][:])*(0.04401/0.02896)
    p_bound = np.squeeze(ifile.variables['pressure'][:])
    ifile.close()

    lats_tf=(lats<=66.5)&(lats>=33.5)
    lons_tf=(lons<=39.5)&(lons>=-13.5)
    lats=lats[np.logical_and(lats>=33.5,lats<=66.5)]
    lons=lons[np.logical_and(lons>=-13.5,lons<=39.5)]

    #lats_tf=(lats<=63.)&(lats>=35.)
    #lons_tf=(lons<=34.5)&(lons>=-10.5)
    #lats=lats[np.logical_and(lats>=35.,lats<=63.)]
    #lons=lons[np.logical_and(lons>=-10.5,lons<=34.5)]

    p_bound=p_bound[:,:,lats_tf,:][:,:,:,lons_tf]
    p=np.empty(shape=(8,25,34,54))
#    p=np.empty(shape=(8,25,15,16))    # was like this before
    for i in range(0,24):
        p[:,i,:,:]=(p_bound[:,i,:,:]+p_bound[:,i+1,:,:])/2.
    sp=np.squeeze(p[:,0,:,:])

    co2=co2[:,:,lats_tf,:][:,:,:,lons_tf]
    co2_ens = np.empty(shape=(co2_bg.size,8,25,34,54))
    for i in range(0, co2_bg.size):
        co2_ens[i,:,:,:,:] = co2[:,:,:,:]*co2_bg[i]

    hyam=np.empty(25)
    hybm=np.empty(25)
    for i in range(0,25):
        hyam[i]=(half_a[i]+half_a[i+1])/2.
        hybm[i]=(half_b[i]+half_b[i+1])/2.

    ttt=("00","03","06","09","12","15","18","21")

    for ti,tt in enumerate(ttt):
        ofile = Dataset(odir+'ct_'+date+tt+'.nc', mode='w')

        olev = ofile.createDimension('level', len(levs))
        olat = ofile.createDimension('lat', len(lats))
        olon = ofile.createDimension('lon', len(lons))
        odate = ofile.createDimension('date', 1)

        olat = ofile.createVariable('lat', np.float64, ('lat',))
        olon = ofile.createVariable('lon', np.float64, ('lon',))
        olev = ofile.createVariable('level', np.float64, ('level',))
        otime = ofile.createVariable('time', np.float64, ('date',))
        odate = ofile.createVariable('date', np.float64, ('date',))

        ohyam = ofile.createVariable('hyam', np.float32, ('level',))
        ohybm = ofile.createVariable('hybm', np.float32, ('level',))

        op = ofile.createVariable('pressure', np.float32, ('date','level','lat','lon'),fill_value=-999.99)
        osp = ofile.createVariable('PSURF', np.float32, ('date','lat','lon'),fill_value=-999.99)
        op0 = ofile.createVariable('P0', np.float32, ('date'),fill_value=-999.99)

        odate.comment = 'time-interval average, centered on times in the date axis'
        odate.long_name = 'UTC dates and times'
        odate.units = 'days since '+date[:4]+'-'+date[4:6]+'-'+date[6:8]+' 00:00:00'
#        otime.dtype = 'double'

        otime.units = 'seconds since '+date[:4]+'-'+date[4:6]+'-'+date[6:8]+' 00:00:00'
        otime.calendar = 'proleptic_gregorian'

        olat.standard_name = 'latitude'
        olat.long_name = 'latitude'
        olat.units = 'degree_north'
        olat.axis = 'Y'

        olon.standard_name = 'longitude'
        olon.long_name = 'longitude'
        olon.units = 'degree_east'
        olon.axis = 'X'

        ohyam.long_name = 'hybrid A coefficient at layer midpoints'
        ohyam.units = 'Pa'

        ohybm.long_name = 'hybrid B coefficient at layer midpoints'
        ohybm.units = '1'

        olev.positive = 'up'
        olev.units = 'levels'

        op.long_name = 'pressure_at_center_levels'
        op.units = 'Pa'
        op.standard_name = 'air pressure'

        op0.units = 'Pa'

        osp.cell_methods = 'level:mean'
        osp.units = 'Pa'
        osp.long_name = 'surface pressure'
        osp.table = '128'
        osp.lev = '1'

        olat[:] = lats
        olon[:] = lons
        olev[:] = levs
        odate[:] = 3.*ti/24.
        otime[:] = ti*3*3600

        osp[:] = sp[ti,:]
        op[:] = p[ti,:]
        op0[:] = 1.
        ohyam[:] = hyam
        ohybm[:] = hybm

        for e in range(0, co2_bg.size):
            ens = str(e).zfill(3)
            oco2 = ofile.createVariable('BG_'+ens, np.float64, ('date','level','lat','lon'),fill_value=-999.99)
            oco2.standard_name = 'mass_fraction_of_carbon_dioxide_in_air'
            oco2.long_name = 'mass mixing ratio of CO2 from outside Europe'
            oco2.units = 'kg kg-1'
            oco2[:] = co2_ens[e,ti,:,:,:]

        ofile.close()
