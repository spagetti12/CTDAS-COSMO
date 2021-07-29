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
# MolefractionToPermil.py

"""
Author : ivar 

Revision History:
File created on 11 May 2012.

"""


def molefraction_to_permil(filename,simulated):
    """ Converts 13C mole fractions to permil values"""
    import da.tools.io4 as io
    import logging
    import numpy as np

    pdb         = 0.011112

    trlength=len(simulated[0,:])
    memlength=trlength/2
    simulated=simulated #*1.e6 #convert to ppm
 #   np.set_printoptions(threshold=np.nan)
 #   msg='simulated shape',simulated.shape;logging.info(msg)
 #   msg='simulated',simulated[0,0],simulated[0,40];logging.info(msg)
 #   msg='simulated',simulated ;logging.info(msg)
    simulated=np.float64(simulated)
    simulated[:,memlength:trlength]=1000.*((simulated[:,memlength:trlength]/simulated[:,0:memlength]/pdb)-1.)
 #   msg='simulated',simulated[0,0],simulated[0,40],memlength,trlength,pdb;logging.info(msg)
 #   msg='simulated',simulated ;logging.info(msg)




    return simulated


if __name__ == "__main__":
    pass
