.. _installing:

Installing CTDAS
======================

.. note:: The current documentation assumes you will run CTDAS with the TM5 transport model. Although this is
          not an absolute requirement, it is advised to start with TM5 before modifying the system to use
          another transport model

To install a fully working CTDAS system on a new system is a two step process:

   1. Install and configure the TM5 transport model
   2. Install and configure the CTDAS shell

The first step is needed because the python based CTDAS shell relies on the FORTRAN based TM5 model to perform transport of tracers given a set of fluxes. Typically, the TM5 model needs to be compiled and tested only once, and is then available for subsequent data assimilation experiments. In principle, any transport model can take the place of TM5 to perform this service but no other transport models are currently supported. Although the current TM5 ctdas code only handles the transport and fluxes of CO2, even a novice TM5 user can make a new project to handle the transport of other tracers (CH4, SF6), or multiple tracers (CO2 + CO + 13CO2).

Installing TM5
--------------

The TM5 transport model is freely available to everyone who consents to the terms of use of the model, and thereby becomes a part of the TM5 user community.

.. important::  
    **Collaboration protocol TM5 model**

    *Rationale*

    The two-way nested TM5 global chemistry transport model has in the last decennia been developed by a consortium of persons 
    and organizations; namely IMAU, KNMI, JRC, and more recently also at NOAA. Whereas there is no formal property right on 
    TM5, and it is in general recognized that an increased user group will be beneficial for the scientific development and 
    credibility of TM5, there may also arise potential problems and conflicts of interest. These may refer to insufficient 
    communication, insufficient acknowledgement of intellectual efforts, duplication of efforts, the use of the same model 
    in competing projects, insufficient feed-back of users on potential problems, publishing scientific results without 
    clarifying what the difference with previous model versions used in previous calculations was, thereby confusing the 
    outside world as to which result to use/or to believe.

    *Rules of conduct*

    The user of TM5 agrees to:

    * To promptly report on bugs or problems in TM5 via e-mail or web-forum.
    * To acknowledge the principal authors of TM5 (e.g. Krol and Segers for the overall models, and several others for 
      sub-modules) in scientific papers or reports. In case that the contribution of a TM5 contributor was essential 
      for the publication co-authorship should be offered.
    * To include in a scientific publication a brief discussion of previous results on a certain topic and to explain 
      the differences.
    * To participate in the bi-annual TM meetings; or have at least one representative of the group to represent him/her.
    * To inform at the TM meeting (or before) the TM5 user group of projects in preparation, submitted and in progress, 
      in order to avoid internal competition.
    * To introduce new users into this protocol, 
    * To closely follow and give scientific and technical support to collaborations with external groups that involve 
      the remote use of TM5, in order to avoid damaging inappropriate use of the model.

To become a member of the TM5 user community, please contact `Wouter Peters at Wageningen University <http://www.maq.wur.nl/UK/Employees/WouterPeters/>`_. 

As a member of the TM5 community you gain access to the TM5 WiKi pages and the TM5 subversion server that allow you to complete an installation of the latest version of the model. The detailed instructions (members only!) are `here <https://www.surfgroepen.nl/sites/tm/Shared%20Wiki/Automatically%20installing%20TM5%20from%20subversion.aspx>`_. Note that in order to run the TM5 CTDAS project (i.e., the fortran code that performs the transport of the CO2 tracer within the overall CTDAS system), you also need to compile extra libraries supported in TM5. Specifically:

    * HDF4 (depends on JPEG, SZIP), http://www.hdfgroup.org/products/hdf4/
    * HDF5 (depends on JPEG, SZIP), http://www.hdfgroup.org/HDF5/
    * NetCDF4 (depends on HDF5, UDUNITS, SZIP), http://www.unidata.ucar.edu/software/netcdf/
    * MPI

Each of the first three libraries has to be built with FORTRAN support, and in addition to a regular build also have to be installed in combination with MPI to support parallel I/O (HDF5, NetCDF4). Also, make sure that HDF4 is compiled without netcdf support. The sets of flags to control this are:
    * for HDF4::

          ./configure \
              --disable-netcdf \
              --enable-fortran \
              --with-jpeg=${JPEG_HOME} \
              --with-szlib=${SZIP_HOME} \
              --disable-shared \
              --with-pic
        
    * for HDF5::

            ./configure \
                --enable-fortran \
                --with-szlib=${SZIP_HOME} \
                --disable-shared \
                --with-pic 

    * for NetCDF4::

            ./configure \
                --enable-netcdf-4 \
                --enable-f90 \
                --disable-dap \
                --with-hdf5=${HDF5_HOME} \
                --with-szlib=${SZIP_HOME} \
                --disable-shared \
                --with-pic 

And when compiling with the MPI compilers, also add the flag::

               --enable-parallel   (for HDF5 compiling)
               --enable-parallel-tests (for NetCDF4 compiling)

At this point, I advise you to continue to the :ref:`tutorial Chapter 0 <tut_chapter0>` for further instructions.


Installing the CTDAS shell
--------------------------

.. note::
    The CTDAS shell is currently not yet available in the open source domain, but is hosted on the password protected server at Wageningen University. Access can be granted to those interested in the project, by sending an email to  `Wouter Peters at Wageningen University <http://www.maq.wur.nl/UK/Employees/WouterPeters/>`_. In the future, we will release CTDAS as a typical python package, possibly through the `Cheeseshop <http://pypi.python.org>`. 

Once you have access to the subversion server, please check out the latest source code to a directory of your choice using::

   svn checkout https://maunaloa.wur.nl/subversion/das/ct/trunk ./

Accept any security certificate permanently (p) and allow the server to store your password unsecure if prompted.

*The CTDAS shell is now installed and ready to be configured!*

At this point, I advise you to continue to the :ref:`tutorial Chapter 1 <tut_chapter1>` for further instructions. 

