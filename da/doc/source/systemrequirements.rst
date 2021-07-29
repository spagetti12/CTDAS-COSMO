.. _systemrequirements:

System Requirements
======================

The CarbonTracker DAS is programmed in python and uses it built-in functionality for many of its tasks.
Users of CTDAS are required to have a python installation on their system, with a small set of 
open source add-ons. The minimum requirements are listed below, categorized by functionality.

Getting the CTDAS code
----------------------

* A working ``subversion`` (SVN) installation is needed, to check your system type::

     $ svn --version

  if the system returns an error, or a version < 1.5.0, please obtain svn from <http://subversion.tigris.org/>


Running CTDAS 
-------------

* A ``python2.4`` or later installation.  

.. note:: ``python3`` is not supported, one needs a python2.x version. To check your python version type::

    $ python --version

* The python module ``numpy``, not included with a standard installation. You can obtain numpy from <http://numpy.scipy.org>

* The python module ``netCDF4`` created by Jeff Whitaker. This package is freely available at <http://code.google.com/p/netcdf4-python/>

Atmospheric Transport
---------------------

An important component of the CTDAS system is the observation operator, usually in the form of an atmospheric
transport model. This model is a stand-alone piece of code that can subsample your state vector, and return a 
set of observations. **You therefore require an atmospheric transport model that can run indepently on your 
platform**. All examples on these pages assume you have access to the TM5 transport model.

.. important:: If you do not currently have access to such a transport model, or have no resources to run such a model
   for long periods of time and for many ensemble configurations, you have probably stumbled onto this page not
   understanding exactly what CTDAS is. We refer you to the :ref:`overview` for a more elaborate description.

The TM5 transport model requires:

    * Fortran 90
    * MPI, HDF4, HDF5, NetCDF4, LAPACK libraries
    * parallel computing capabilities (16-100 CPUs recommended)
    * Meteorological driver data (260 Gb of storage for each year of input data)

Analysis of results
-------------------

Although many packages can be used to perform analyses (IDL, MatLAB, GRADS, ...) we have included standard routines for python. Recommended system configuration for these is:
    * python 2.5 or higher <http://www.python.org/>
    * numpy <http://numpy.scipy.org>
    * matplotlib <http://matplotlib.sourceforge.net/>
    * basemap <http://matplotlib.sourceforge.net/basemap/doc/html/>
    * netCDF4 <http://code.google.com/p/netcdf4-python/>




