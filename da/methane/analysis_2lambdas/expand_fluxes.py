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
# expand_fluxes.py
import sys
sys.path.append('../../')
import os
from datetime import datetime, timedelta

import logging
import numpy as np
from da.tools.general import date2num, num2date
import da.methane.io4 as io
from da.methane.analysis.tools_regions import globarea,state_to_grid
from da.tools.general import create_dirs
    


"""
Author: aki

Revision History:
File created on April 2016 by Aki T.

"""

def save_weekly_avg_1x1_data(dacycle, statevector):
    """
        Function creates a NetCDF file with output on 1x1 degree grid. It uses the flux data written by the 
        :class:`~da.baseclasses.obsoperator.ObsOperator.py`, and multiplies these with the mapped parameters and
        variance (not covariance!) from the :class:`~da.baseclasses.statevector.StateVector`.
        
           :param dacycle: a :class:`~da.tools.initexit.CycleControl` object
           :param statevector: a :class:`~da.baseclasses.statevector.StateVector`
           :rtype: None
    """
#
    dirname = create_dirs(os.path.join(dacycle['dir.analysis'], 'data_flux1x1_weekly'))
#
# Some help variables
#
    dectime0 = date2num(datetime(2000, 1, 1))
    dt = dacycle['cyclelength']
    startdate = dacycle['time.start'] 
    enddate = dacycle['time.end'] 
    nlag = statevector.nlag

    logging.debug("DA Cycle start date is %s" % startdate.strftime('%Y-%m-%d %H:%M'))
    logging.debug("DA Cycle end   date is %s" % enddate.strftime('%Y-%m-%d %H:%M'))

#
# Region bio_land and fossil_land
# 

    region_land_file = dacycle.dasystem['regionsfile']
    nfcr = io.CT_CDF(region_land_file, 'read')
    logging.info('region land file read: %s' %region_land_file)
    bio_land = nfcr.variables['bio_land'][:]
    fossil_land = nfcr.variables['fossil_land'][:]
    #logging.debug("bio_land %s" %bio_land[25,205])
    #logging.debug("ff_land %s" %fossil_land[25,205])
   
#
# Create or open NetCDF output file
#
    logging.debug("Create NetCDF output file: flux_1x1.%s.nc" % startdate.strftime('%Y-%m-%d'))
    saveas = os.path.join(dirname, 'flux_1x1.%s.nc' % startdate.strftime('%Y-%m-%d'))
    ncf = io.CT_CDF(saveas, 'write')

#
# Create dimensions and lat/lon grid
#
    dimgrid = ncf.add_latlon_dim()
    dimensemble = ncf.add_dim('members', statevector.nmembers)
    dimdate = ncf.add_date_dim()
#
# set title and tell GMT that we are using "pixel registration"
#
    setattr(ncf, 'Title', 'CarbonTracker fluxes')
    setattr(ncf, 'node_offset', 1)
#
# skip dataset if already in file
#
    ncfdate = date2num(startdate) - dectime0 + dt.days / 2.0
    skip = ncf.has_date(ncfdate)
    if skip:
        logging.warning('Skipping writing of data for date %s : already present in file %s' % (startdate.strftime('%Y-%m-%d'), saveas))
    else:
        
#
# if not, process this cycle. Start by getting flux input data from CTDAS
#
        filename = os.path.join(dacycle['dir.output'], 'flux1x1_%s_%s.nc' % (startdate.strftime('%Y%m%d%H'), enddate.strftime('%Y%m%d%H')))

        file = io.ct_read(filename, 'read')
        bio = np.array(file.get_variable(dacycle.dasystem['background.ch4.bio.flux']))
        ocean = np.array(file.get_variable(dacycle.dasystem['background.ch4.ocean.flux']))
        fire = np.array(file.get_variable(dacycle.dasystem['background.ch4.fires.flux']))
        fossil = np.array(file.get_variable(dacycle.dasystem['background.ch4.fossil.flux']))
        term = np.array(file.get_variable(dacycle.dasystem['background.ch4.term.flux']))
        #mapped_parameters   = np.array(file.get_variable(dacycle.dasystem['final.param.mean.1x1']))
        file.close()

        next = ncf.inq_unlimlen()[0]


# Start adding datasets from here on, both prior and posterior datasets for bio and ocn

        for prior in [True, False]:
#
# Now fill the statevector with the prior values for this time step. Note that the prior value for this time step
# occurred nlag time steps ago, so we make a shift in the output directory, but only if we are more than nlag cycle away from the start date..
#

            if prior:
                logging.debug("prior")
                qual_short = 'prior'
                for n in range(nlag, 0, -1):
                    priordate = enddate - timedelta(dt.days * n)
                    savedir = dacycle['dir.output'].replace(startdate.strftime('%Y%m%d'), priordate.strftime('%Y%m%d'))
                    filename = os.path.join(savedir, 'savestate_%s.nc' % priordate.strftime('%Y%m%d'))
                    if os.path.exists(filename):
                        statevector.read_from_file(filename, qual=qual_short)
                        gridmean, gridensemble = statevector.state_to_grid(lag=n)

# Replace the mean statevector by all ones (assumed priors)

                        gridmean = statevector.vector2grid(vectordata=np.ones(statevector.nparams,))

                        logging.debug('Read prior dataset from file %s, sds %d: ' % (filename, n))
                        break
            else:
                logging.debug("opt")
                qual_short = 'opt'
                savedir = dacycle['dir.output']
                filename = os.path.join(savedir, 'savestate_%s.nc' % startdate.strftime('%Y%m%d'))
                statevector.read_from_file(filename, qual=qual_short)
                gridmean, gridensemble = statevector.state_to_grid(lag=1)

                logging.debug('Read posterior dataset from file %s, sds %d: ' % (filename, 1))
#
# if prior, do not multiply fluxes with parameters, otherwise do
#
            logging.debug("gridmean %s" %gridmean[25,205])
            logging.debug("bio land %s" %bio_land[25,205])
            w = np.where(bio_land == 0)
            s = gridmean.shape
            gridmean_bio = np.zeros(shape=s)
            gridmean_bio[:,:] = gridmean[:,:]
            gridmean_bio[w] = 1.0
            biomapped = bio * gridmean_bio
            s = gridensemble.shape
            gridensemble_bio = np.zeros(shape=s)
            gridensemble_bio[:,:,:] = gridensemble[:,:,:]
            gridensemble_bio[:, w[0],w[1]] = 0.0
            biovarmapped = bio * gridensemble_bio

            w = np.where(fossil_land == 0)
            gridmean_fossil = gridmean
            gridmean_fossil[w] = 1.0
            fossilmapped = fossil * gridmean_fossil
            gridensemble_fossil = gridensemble
            gridensemble_fossil[:, w[0],w[1]] = 0.0
            fossilvarmapped = fossil * gridensemble_fossil

            logging.debug("gridmean %s" %gridmean[25,205])
            logging.debug("gridmean bio %s" %gridmean_bio[25,205])
            logging.debug("gridmean ff %s" %gridmean_fossil[25,205])
#
#
#  For each dataset, get the standard definitions from the module mysettings, add values, dimensions, and unlimited count, then write
#
            savedict = ncf.standard_var(varname='bio_flux_' + qual_short)
            savedict['values'] = biomapped.tolist()
            savedict['dims'] = dimdate + dimgrid
            savedict['count'] = next
            ncf.add_data(savedict)
#
            savedict = ncf.standard_var(varname='fossil_flux_' + qual_short)
            savedict['values'] = fossilmapped.tolist()
            savedict['dims'] = dimdate + dimgrid
            savedict['count'] = next
            ncf.add_data(savedict)

            #print biovarmapped.shape
            savedict = ncf.standard_var(varname='bio_flux_%s_ensemble' % qual_short)
            savedict['values'] = biovarmapped.tolist()
            savedict['dims'] = dimdate + dimensemble + dimgrid
            savedict['count'] = next
            ncf.add_data(savedict)
#
            savedict = ncf.standard_var(varname='fossil_flux_%s_ensemble' % qual_short)
            savedict['values'] = fossilvarmapped.tolist()
            savedict['dims'] = dimdate + dimensemble + dimgrid
            savedict['count'] = next
            ncf.add_data(savedict)

        # End prior/posterior block

        savedict = ncf.standard_var(varname='fire_flux_imp')
        savedict['values'] = fire.tolist()
        savedict['dims'] = dimdate + dimgrid
        savedict['count'] = next
        ncf.add_data(savedict)
#
        savedict = ncf.standard_var(varname='ocn_flux_imp')
        savedict['values'] = ocean.tolist()
        savedict['dims'] = dimdate + dimgrid
        savedict['count'] = next
        ncf.add_data(savedict)

        savedict = ncf.standard_var(varname='term_flux_imp')
        savedict['values'] = term.tolist()
        savedict['dims'] = dimdate + dimgrid
        savedict['count'] = next
        ncf.add_data(savedict)

        #area = globarea()
        #savedict = ncf.standard_var(varname='cell_area')
        #savedict['values'] = area.tolist()
        #savedict['dims'] = dimgrid
        #ncf.add_data(savedict)
#
        savedict = ncf.standard_var(varname='date')
        savedict['values'] = date2num(startdate) - dectime0 + dt.days / 2.0
        savedict['dims'] = dimdate
        savedict['count'] = next
        ncf.add_data(savedict)

        sys.stdout.write('.')
        sys.stdout.flush()
#
#   Done, close the new NetCDF file
#
    ncf.close()
#
#   Return the full name of the NetCDF file so it can be processed by the next routine
#
    logging.info("Gridded weekly average fluxes now written")

    return saveas

def save_weekly_avg_state_data(dacycle, statevector):
    """
        Function creates a NetCDF file with output for all parameters. It uses the flux data written by the 
        :class:`~da.baseclasses.obsoperator.ObsOperator.py`, and multiplies these with the mapped parameters and
        variance (not covariance!) from the :class:`~da.baseclasses.statevector.StateVector`.
        
           :param dacycle: a :class:`~da.tools.initexit.CycleControl` object
           :param statevector: a :class:`~da.baseclasses.statevector.StateVector`
           :rtype: None
    """
    logging.debug('start: save weekly avg state data')
    dirname = create_dirs(os.path.join(dacycle['dir.analysis'], 'data_state_weekly'))
#
# Some help variables
#
    dectime0 = date2num(datetime(2000, 1, 1))
    dt = dacycle['cyclelength']
    startdate = dacycle['time.start'] 
    enddate = dacycle['time.end'] 
    nlag = statevector.nlag

    area = globarea()
    vectorarea = statevector.grid2vector(griddata=area, method='sum')

    logging.debug("DA Cycle start date is %s" % startdate.strftime('%Y-%m-%d %H:%M'))
    logging.debug("DA Cycle end   date is %s" % enddate.strftime('%Y-%m-%d %H:%M'))

#
# Create or open NetCDF output file
#
    saveas = os.path.join(dirname, 'statefluxes.nc')
    ncf = io.CT_CDF(saveas, 'write')

    logging.debug('save weekly avg state data to file %s' %saveas)
#
# Create dimensions and lat/lon grid
#
    dimregs = ncf.add_dim('nparameters', statevector.nparams)
    dimmembers = ncf.add_dim('nmembers', statevector.nmembers)
    dimdate = ncf.add_date_dim()
#
# set title and tell GMT that we are using "pixel registration"
#
    setattr(ncf, 'Title', 'CarbonTracker fluxes')
    setattr(ncf, 'node_offset', 1)
#
# skip dataset if already in file
#
    ncfdate = date2num(startdate) - dectime0 + dt.days / 2.0
    skip = ncf.has_date(ncfdate)
    if skip:
        logging.warning('Skipping writing of data for date %s : already present in file %s' % (startdate.strftime('%Y-%m-%d'), saveas))
    else:
        next = ncf.inq_unlimlen()[0]
        logging.debug('Writing of data for date %s: file %s' % (startdate.strftime('%Y-%m-%d'), saveas))

#
# if not, process this cycle. Start by getting flux input data from CTDAS
#
        filename = os.path.join(dacycle['dir.output'], 'flux1x1_%s_%s.nc' % (startdate.strftime('%Y%m%d%H'), enddate.strftime('%Y%m%d%H')))

        file = io.ct_read(filename, 'read')
        bio = np.array(file.get_variable(dacycle.dasystem['background.ch4.bio.flux']))
        ocean = np.array(file.get_variable(dacycle.dasystem['background.ch4.ocean.flux']))
        fire = np.array(file.get_variable(dacycle.dasystem['background.ch4.fires.flux']))
        fossil = np.array(file.get_variable(dacycle.dasystem['background.ch4.fossil.flux']))
        term = np.array(file.get_variable(dacycle.dasystem['background.ch4.term.flux']))
        #mapped_parameters   = np.array(file.get_variable(dacycle.dasystem['final.param.mean.1x1']))
        file.close()

        next = ncf.inq_unlimlen()[0]

        vectorbio = statevector.grid2vector(griddata=bio * area, method='sum')
        vectorocn = statevector.grid2vector(griddata=ocean * area, method='sum')
        vectorfire = statevector.grid2vector(griddata=fire * area, method='sum')
        vectorfossil = statevector.grid2vector(griddata=fossil * area, method='sum')
        vectorterm = statevector.grid2vector(griddata=term * area, method='sum')


# Start adding datasets from here on, both prior and posterior datasets for bio and ocn

        for prior in [True, False]:
#
# Now fill the statevector with the prior values for this time step. Note that the prior value for this time step
# occurred nlag time steps ago, so we make a shift in the output directory, but only if we are more than nlag cycle away from the start date..
#
            if prior:
                qual_short = 'prior'
                for n in range(nlag, 0, -1):
                    priordate = enddate - timedelta(dt.days * n)
                    savedir = dacycle['dir.output'].replace(startdate.strftime('%Y%m%d'), priordate.strftime('%Y%m%d'))
                    filename = os.path.join(savedir,'savestate_%s.nc' % priordate.strftime('%Y%m%d'))
                    if os.path.exists(filename):
                        statevector.read_from_file(filename, qual=qual_short)
# Replace the mean statevector by all ones (assumed priors)
                        statemean = np.ones((statevector.nparams,))
                        choicelag = n
                        logging.debug('Read prior dataset from file %s, lag %d: ' % (filename, choicelag))
                        break
            else:
                qual_short = 'opt'
                savedir = dacycle['dir.output']
                filename = os.path.join(savedir, 'savestate_%s.nc' % startdate.strftime('%Y%m%d'))
                statevector.read_from_file(filename)
                choicelag = 1
                statemean = statevector.ensemble_members[choicelag - 1][0].param_values
                logging.debug('Read posterior dataset from file %s, lag %d: ' % (filename, choicelag))
#
# if prior, do not multiply fluxes with parameters, otherwise do
#
            data = statemean * vectorbio # units of mole region-1 s-1

            savedict = ncf.standard_var(varname='bio_flux_%s' % qual_short)
            savedict['values'] = data
            savedict['dims'] = dimdate + dimregs
            savedict['count'] = next
            ncf.add_data(savedict)

#
# Here comes a special provision for the posterior flux covariances: these are calculated relative to the prior flux covariance to
# ensure they are indeed smaller due to the data assimilation. If they would be calculated relative to the mean posterior flux, the 
# uncertainties would shift just because the mean flux had increased or decreased, which is not what we want.
#
# The implementation is done by multiplying the ensemble with the vectorbio only, and not with the statemean values
# which are assumed 1.0 in the prior always.
#

            members = statevector.ensemble_members[choicelag - 1] 
            deviations = np.array([mem.param_values * vectorbio for mem in members])
            deviations = deviations - deviations[0, :]

            savedict = ncf.standard_var(varname='bio_flux_%s_ensemble' % qual_short)

            savedict['values'] = deviations.tolist()
            savedict['dims'] = dimdate + dimmembers + dimregs
            savedict['comment'] = "This is the matrix square root, use (M x M^T)/(nmembers-1) to make covariance"
            savedict['units'] = "mol region-1 s-1"
            savedict['count'] = next
            ncf.add_data(savedict)

            savedict = ncf.standard_var('unknown')
            savedict['name'] = 'bio_flux_%s_std' % qual_short
            savedict['long_name'] = 'Biosphere flux standard deviation, %s' % qual_short
            savedict['values'] = deviations.std(axis=0)
            savedict['dims'] = dimdate + dimregs
            savedict['comment'] = "This is the standard deviation on each parameter"
            savedict['units'] = "mol region-1 s-1"
            savedict['count'] = next
            ncf.add_data(savedict)

            data = statemean * vectorfossil # units of mole region-1 s-1

            savedict = ncf.standard_var(varname='fossil_flux_%s' % qual_short)
            savedict['values'] = data
            savedict['dims'] = dimdate + dimregs
            savedict['count'] = next
            ncf.add_data(savedict)


#
# Here comes a special provision for the posterior flux covariances: these are calculated relative to the prior flux covariance to
# ensure they are indeed smaller due to the data assimilation. If they would be calculated relative to the mean posterior flux, the 
# uncertainties would shift just because the mean flux had increased or decreased, which is not what we want.
#
# The implementation is done by multiplying the ensemble with the vectorocn only, and not with the statemean values
# which are assumed 1.0 in the prior always.
#

            deviations = np.array([mem.param_values * vectorfossil for mem in members])
            deviations = deviations - deviations[0, :]

            savedict = ncf.standard_var(varname='fossil_flux_%s_ensemble' % qual_short)
            savedict['values'] = deviations.tolist()
            savedict['dims'] = dimdate + dimmembers + dimregs
            savedict['comment'] = "This is the matrix square root, use (M x M^T)/(nmembers-1) to make covariance"
            savedict['units'] = "mol region-1 s-1"
            savedict['count'] = next
            ncf.add_data(savedict)

            savedict = ncf.standard_var('unknown')
            savedict['name'] = 'fossil_flux_%s_std' % qual_short
            savedict['long_name'] = 'Fossil flux standard deviation, %s' % qual_short
            savedict['values'] = deviations.std(axis=0)
            savedict['dims'] = dimdate + dimregs
            savedict['comment'] = "This is the standard deviation on each parameter"
            savedict['units'] = "mol region-1 s-1"
            savedict['count'] = next
            ncf.add_data(savedict)

        data = vectorfire

        savedict = ncf.standard_var(varname='fire_flux_imp')
        savedict['values'] = data
        savedict['dims'] = dimdate + dimregs
        savedict['count'] = next
        ncf.add_data(savedict)

        data = vectorterm

        savedict = ncf.standard_var(varname='term_flux_imp')
        savedict['values'] = data
        savedict['dims'] = dimdate + dimregs
        savedict['count'] = next
        ncf.add_data(savedict)

        data = vectorocn

        savedict = ncf.standard_var(varname='ocn_flux_imp')
        savedict['values'] = data
        savedict['dims'] = dimdate + dimregs
        savedict['count'] = next
        ncf.add_data(savedict)

        savedict = ncf.standard_var(varname='date')
        savedict['values'] = ncfdate
        savedict['dims'] = dimdate
        savedict['count'] = next
        ncf.add_data(savedict)

        sys.stdout.write('.')
        sys.stdout.flush()
#
#   Done, close the new NetCDF file
#
    ncf.close()
#
#   Return the full name of the NetCDF file so it can be processed by the next routine
#
    logging.info("Vector weekly average fluxes now written")

    return saveas


