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

import numpy as np
import cPickle

def state_to_grid(values, regionmap, reverse=False, avg=False, mapname=None):
    """ 
    This method converts parameters from a CarbonTracker StateVector object to a gridded map of linear multiplication values. These
    can subsequently be used in the transport model code to multiply/manipulate fluxes

    """
    nregions = regionmap.max()
    try:
        if not mapname: 
            raise Exception 

        regionselect = cPickle.load(open('%s_regiondict.pickle' % mapname, 'rb'))
    except:

        # dictionary for region <-> map conversions
        regs = {}
        for r in np.arange(1, nregions + 1):
            sel = (regionmap.flat == r).nonzero()
            if len(sel[0]) > 0: 
                regs[r] = sel

        regionselect = regs
        
        cPickle.dump(regionselect, open('%s_regiondict.pickle' % mapname, 'wb'), -1)
        print 'Pickling region map'

    if reverse:
        """ project 1x1 degree map onto ecoregions """

        result = np.zeros(nregions, float)
        for k, v in regionselect.iteritems():
            if avg: 
                result[k - 1] = values.ravel().take(v).mean()
            else : 
                result[k - 1] = values.ravel().take(v).sum()
        return result

    else:
        """ project ecoregion properties onto 1x1 degree map """

        result = np.zeros((180, 360,), float)
        for k, v in regionselect.iteritems():
            result.put(v, values[k - 1])

        return result

def globarea(im=360, jm=180, silent=True):
    """ Function calculates the surface area according to TM5 definitions"""

    radius = 6.371e6  # the earth radius in meters
    deg2rad = np.pi / 180.
    g = 9.80665 

    dxx = 360.0 / im * deg2rad 
    dyy = 180.0 / jm * deg2rad 
    lat = np.arange(-90 * deg2rad, 90 * deg2rad, dyy)
    dxy = dxx * (np.sin(lat + dyy) - np.sin(lat)) * radius ** 2
    area = np.resize(np.repeat(dxy, im, axis=0) , [jm, im])
    if not silent:
        print 'total area of field = ', np.sum(area.flat)
        print 'total earth area    = ', 4 * np.pi * radius ** 2
    return area

