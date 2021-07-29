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
# test_optimizer.py

"""
Author : peters 

Revision History:
File created on 04 Aug 2010.

"""

def serial_py_against_serial_fortran():
    """ Test the solution of the serial algorithm against the CT cy2 fortran generated one """

    # get data from the savestate.hdf file from the first cycle of CarbonTracker 2009 release

    print "WARNING: The optimization algorithm has changed from the CT2009 release because of a bug"
    print "WARNING: in the fortran code. Hence, the two solutions calculated are no longer the same."
    print "WARNING: To change the python algorithm so that it corresponds to the fortran, change the"
    print "WARNING: loop from m=n+1,nfcast to m=1,nfcast"

    savefile = '/data/CO2/peters/carbontracker/raw/ct09rc0i/20000108/savestate.hdf'
    print savefile

    f = Nio.open_file(savefile, 'r')
    obs = f.variables['co2_obs_fcast'].get_value()
    sel_obs = obs.shape[0]

    dims = (int(dacycle.da_settings['time.nlag']),
                  int(dacycle.da_settings['forecast.nmembers']),
                  int(dacycle.dasystem.da_settings['nparameters']),
                  sel_obs,)

    nlag, nmembers, nparams, nobs = dims

    optserial = CtOptimizer(dims)
    opt = optserial

    opt.set_localization('CT2007')

    obs = f.variables['co2_obs_fcast'].get_value()[0:nobs]
    opt.obs = obs
    sim = f.variables['co2_sim_fcast'].get_value()[0:nobs]
    opt.Hx = sim
    error = f.variables['error_sim_fcast'].get_value()[0:nobs]
    flags = f.variables['flag_sim_fcast'].get_value()[0:nobs]
    opt.flags = flags
    simana = f.variables['co2_sim_ana'].get_value()[0:nobs]

    for n in range(nobs): opt.R[n, n] = np.double(error[n] ** 2)

    xac = []
    adX = []
    for lag in range(nlag):
        xpc = f.variables['xpc_%02d' % (lag + 1)].get_value()
        opt.x[lag * nparams:(lag + 1) * nparams] = xpc
        X = f.variables['pdX_%02d' % (lag + 1)].get_value()
        opt.X_prime[lag * nparams:(lag + 1) * nparams, :] = np.transpose(X)
        HX = f.variables['dF'][:, 0:sel_obs]
        opt.HX_prime[:, :] = np.transpose(HX)

        # Also create arrays of the analysis of the fortran code for later comparison

        xac.extend (f.variables['xac_%02d' % (lag + 1)].get_value())
        adX.append (f.variables['adX_%02d' % (lag + 1)].get_value())

    xac = np.array(xac)
    X_prime = np.array(adX).swapaxes(1, 2).reshape((opt.nparams * opt.nlag, opt.nmembers))

    opt.serial_minimum_least_squares()

    print "Maximum differences and correlation of 2 state vectors:"
    print np.abs(xac - opt.x).max(), np.corrcoef(xac, opt.x)[0, 1]
       
    plt.figure(1)
    plt.plot(opt.x, label='SerialPy')
    plt.plot(xac, label='SerialFortran')
    plt.grid(True)
    plt.legend(loc=0)
    plt.title('Analysis of state vector')

    print "Maximum differences of 2 state vector deviations:"
    print np.abs(X_prime - opt.X_prime).max()

    plt.figure(2)
    plt.plot(opt.X_prime.flatten(), label='SerialPy')
    plt.plot(X_prime.flatten(), label='SerialFortran')
    plt.grid(True)
    plt.legend(loc=0)
    plt.title('Analysis of state vector deviations')

    print "Maximum differences and correlation of 2 simulated obs vectors:"
    print np.abs(simana - opt.Hx).max(), np.corrcoef(simana, opt.Hx)[0, 1]

    plt.figure(3)
    plt.plot(opt.Hx, label='SerialPy')
    plt.plot(simana, label='SerialFortran')
    plt.grid(True)
    plt.legend(loc=0)
    plt.title('Analysis of CO2 mole fractions')
    plt.show()

    f.close()

def serial_vs_bulk():
    """ A test of the two algorithms currently implemented: serial vs bulk solution """    

    # get data from the savestate.hdf file from the first cycle of CarbonTracker 2009 release

    savefile = '/data/CO2/peters/carbontracker/raw/ct09rc0i/20000108/savestate.hdf'
    print savefile

    f = Nio.open_file(savefile, 'r')
    obs = f.variables['co2_obs_fcast'].get_value()

    nobs = 77

    dims = (int(dacycle.da_settings['time.nlag']),
                  int(dacycle.da_settings['forecast.nmembers']),
                  int(dacycle.dasystem.da_settings['nparameters']),
                  nobs,)

    nlag, nmembers, nparams, nobs = dims

    optbulk = CtOptimizer(dims)
    optserial = CtOptimizer(dims)

    for o, opt in enumerate([optbulk, optserial]):

        opt.set_localization('CT2007')

        obs = f.variables['co2_obs_fcast'].get_value()[0:nobs]
        opt.obs = obs
        sim = f.variables['co2_sim_fcast'].get_value()[0:nobs]
        opt.Hx = sim
        error = f.variables['error_sim_fcast'].get_value()[0:nobs]
        flags = f.variables['flag_sim_fcast'].get_value()[0:nobs]
        opt.flags = flags

        for n in range(nobs): 
            opt.R[n, n] = np.double(error[n] ** 2)

        xac = []
        for lag in range(nlag):
            xpc = f.variables['xpc_%02d' % (lag + 1)].get_value()
            opt.x[lag * nparams:(lag + 1) * nparams] = xpc
            X = f.variables['pdX_%02d' % (lag + 1)].get_value()
            opt.X_prime[lag * nparams:(lag + 1) * nparams, :] = np.transpose(X)
            HX = f.variables['dF'][:, 0:nobs]
            opt.HX_prime[:, :] = np.transpose(HX)

        if o == 0:
            opt.bulk_minimum_least_squares()
            x1 = opt.x
            xp1 = opt.X_prime
            hx1 = opt.Hx
            hxp1 = opt.HX_prime
            hphr1 = opt.HPHR
            k1 = opt.KG
        if o == 1:
            opt.serial_minimum_least_squares()
            x2 = opt.x
            xp2 = opt.X_prime
            hx2 = opt.Hx
            hxp2 = opt.HX_prime
            hphr2 = opt.HPHR
            k2 = opt.KG
           
    plt.figure()

    print "Maximum differences and correlation of 2 state vectors:"
    print np.abs(x2 - x1).max(), np.corrcoef(x2, x1)[0, 1]
       
    plt.figure(1)
    plt.plot(x1, label='Serial')
    plt.plot(x2, label='Bulk')
    plt.grid(True)
    plt.legend(loc=0)
    plt.title('Analysis of state vector')

    print "Maximum differences of 2 state vector deviations:"
    print np.abs(xp2 - xp1).max()

    plt.figure(2)
    plt.plot(xp1.flatten(), label='Serial')
    plt.plot(xp2.flatten(), label='Bulk')
    plt.grid(True)
    plt.legend(loc=0)
    plt.title('Analysis of state vector deviations')

    print "Maximum differences and correlation of 2 simulated obs vectors:"
    print np.abs(hx2 - hx1).max(), np.corrcoef(hx2, hx1)[0, 1]

    plt.figure(3)
    plt.plot(hx1, label='Serial')
    plt.plot(hx2, label='Bulk')
    plt.title('Analysis of CO2 mole fractions')
    plt.grid(True)
    plt.legend(loc=0)

    plt.show()

    f.close()



if __name__ == "__main__":
    pass
