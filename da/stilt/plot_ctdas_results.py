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
"""
Purpose:
    To check the optimized results of CTDAS
    1. calculate the post state vector and  covriance;
    2. compare the prior and posterior covariances to illustrate the uncertainty reduction after optimization
    3. show H(Xa) is closer to observation than the prior, also the uncertainty in mole fraction is smaller
       after optimization

Date: Jan 26, 2015
Author: W. He
"""

#!/usr/bin/python
from pylab import *
import matplotlib.pyplot as plt
import glob
import netCDF4 as nc

import numpy as np
import datetime
import matplotlib.ticker as ticker

file='/Storage/CO2/wei/ctdas-stilt-proto/output/20100121/optimizer.20100121.nc'
ncf=nc.Dataset(file)

var1=[]
var1.append(ncf.variables['statevectormean_optimized'][:])

from numpy import array
var1=array(var1)
var1.shape

var2=[]
var2.append(ncf.variables['statevectordeviations_optimized'][:]) # for all ensenble members
var2=array(var2)
var2.shape

################################################
#calculate Xa using equation (5)
################################################
clf()

v1=var1[0,:]
plot(v1,ls='-',color='r')

v2=var2[0,:,9]   # only the 10th member
plot(v1+v2,ls='--',color='g')

plt.title('post state vector: 10th member')
plt.xlabel('state id')
plt.ylabel('value')
plt.legend(( 'X','Xa'), loc = 0)
savefig('poststatevector.pdf')

################################################
#calculate Pa using equation (6) and (7)
################################################
clf()
# define Matrix X using normalized ensemble of deviations
N=10
Xa=[]
Xa=1/sqrt(N-1)*var2

Pa=[]
Pa=Xa*np.transpose(Xa)

v3=Pa[:,8,:]   # only the 9th state

#a=plt.imshow(v3,extent=[1,1,10,10],interpolation='nearest')
a=plt.imshow(v3,interpolation='none')
plt.xlabel('ensemble id')
plt.ylabel('ensemble id')
plt.title('post covariance: 9th state')
plt.colorbar()

fig1 = plt.gcf()
fig1.savefig('mappostcovar.pdf',dpi=600)
plt.close()

################################################
#calculate Pb using equation (6) and (7)
################################################

var3=[]
var3.append(ncf.variables['statevectormean_prior'][:])

from numpy import array
var3=array(var3)
var3.shape

clf()
var4=[]
var4.append(ncf.variables['statevectordeviations_prior'][:]) # for all ensenble members
var4=array(var4)
var4.shape

v33=var3[0,:]
v44=var4[0,:,9]   # only the 10th member
plot(v33+v44,ls='-',color='r')
plot(v1+v2,ls='--',color='g')

plt.title('state vector: 10th member')
plt.xlabel('state id')
plt.ylabel('value')
plt.legend(( 'Xb','Xa'), loc = 0)
savefig('statevector.pdf')

#
clf()
Xb=[]
Xb=1/sqrt(N-1)*var4

Pb=[]
Pb=Xb*np.transpose(Xb)

v4=Pb[:,8,:]   # only the 9th state

a=plt.imshow(v4,interpolation='none')
plt.xlabel('ensemble id')
plt.ylabel('ensemble id')
plt.title('prior covariance: 9th state')
plt.colorbar()

fig1 = plt.gcf()
fig1.savefig('mappriorcovar.pdf',dpi=600)
plt.close()

#show figure together
v5=Pb[9,:,9]
v6=Pa[9,:,9]
plot(v5,ls='-',color='r')
plot(v6,ls='--',color='g')

plt.title('state covariance vector: 10th member')
plt.xlabel('state id')
plt.ylabel('value')
plt.legend(( 'Pb','Pa'), loc = 0)
savefig('statecovariance.pdf')

################################################
# figure H(Xa),H(Xb) and observations
################################################
clf()
var7=[]
var7.append(ncf.variables['modelsamplesmean_prior'][:]) # for all ensenble members
var7=array(var7)
var7.shape

var8=[]
var8.append(ncf.variables['modelsamplesmean_optimized'][:]) # for all ensenble members
var8=array(var8)
var8.shape

var9=[]
var9.append(ncf.variables['observed'][:]) # for all ensenble members
var9=array(var9)
var9.shape

v7=var7[0,:]*1e6
v8=var8[0,:]*1e6
v9=var9[0,:]*1e6

plot(v7,ls='-',color='r')
plot(v8,ls='-',color='g')
plot(v9,ls='-',color='b')
plt.title('model samples mean')
plt.xlabel('state id')
plt.ylabel('value')
plt.legend(( 'Prior','Posterior','Observed'), loc = 0)
savefig('modelsamplesmean.pdf')


