.. _tut_chapter0:

Chapter 0: Installing the transport model code
----------------------------------------------------

**In which you will learn to setup the TM5 model so it can be used within CTDAS**

0.1 Installing the TM5 base trunk with pycasso
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We recommend to test a base TM5 version first and to make sure that all libraries compile and work. The instructions for the current trunk of TM5 (note that there might be minor changes as the TM5 trunk evolves):

Make a clean directory on your machine to hold TM5, and later also CTDAS
from that directory, grab the TM5 install file ::

  svn export https://svn.knmi.nl/svn/TM5/install_tm5 

run the install script::

  ./install_tm5

you minimally need to check out from TM5::

  base/trunk
  grids/glb600x400/trunk
  levels/ml60/trunk
  user_output/trunk
  enkf/branches/ctdas

Once done, make a link to the pycasso runscript::

  ln -s  base/trunk/bin/pycasso_setup_tm5 setup_tm5

If needed, create a pycasso-machine-xxxx.rc file for your architecture::

  cd base/trunk/rc
  cp pycasso-machine-template.rc pycasso-machine-[yourmachinehere].rc

And edit it such that all keys have appropriate values

Return to the main directory and copy the pycasso-tm5.rc example file::

  cp base/trunk/rc/pycasso-tm5.rc .

Edit this file, specifically the keys::

  my.machine.rc       :  pycasso-machine-[yourmachinehere].rc
  my.tm5.define       :  slopes with_pycasso without_chemistry without_dry_deposition without_wet_deposition without_sedimentation

Now, try to compile and link the model, read the log file to work your way past possible error messages::

  ./setup_tm5 pycasso-tm5.rc

Once the model has successfully compiled, try to run it, again using the log messages to resolve issues::

  ./setup_tm5 -f -s pycasso-tm5.rc

After a successful completion on one processor, try the code on multiple processors by changing the pycasso-tm5.rc keys::

  par.mpi    :  T
  par.ntask  :  2

If all works fine, you're ready to try running the TM5 CTDAS project...
        
0.2 Installing the TM5 CTDAS code
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Make a copy of the tm5-ctdas.rc example file::

    cp proj/enkf/branches/ctdas/rc/tm5-ctdas.rc .

Edit the tm5-ctdas.rc file to reflect the following keys::

    my.machine.rc       :  pycasso-machine-[yourmachinehere].rc
    my.tm5.define       :  slopes with_pycasso without_chemistry without_dry_deposition without_wet_deposition without_sedimentation

And pay special attention to the locations of CTDAS input files that are needed (modify the path, and perhaps the filename)::

    ! This is the location of an observations file to read flask data from

    obsdir              : /Volumes/Storage/CO2/carbontracker/input/obsnc
    output.flask.infile : ${obsdir}/obs_forecast.nc

    ! This is the location of carbontracker a-priori fluxes for ocean, fores, and biosphere. 
    ! Set the params dir and file to none to simply use prior scaling factors of 1.0

    emis.input.dir          : /Volumes/Storage/CO2/carbontracker/input/ct09
    ct.params.input.dir     : None
    ct.params.input.file    : None
   
    ! This is the location and name of the fossil fuel input files

    ff_hdf.prefix           : ${emis.input.dir}/co2ff.country_cdiac_bp.seasnameuras.1x1.hdf
    ff_hdf.varname          : ff_country_new2009b.pro200907271745
 

If this project compiles successfully and completes a day without major warnings or errors, you are ready to proceed to the next step

.. warning::
    If you are planning to run the code on a machine with limited memory, or you are running it one only one processor, you might consider
    reducing the number of tracers (ensemble members) in the project. To do so, edit the file::
        proj/enkf/branches/ctdas/inv/src/chem_param.F90

This concludes the installation of the TM5 transport model code. Now continue installing the CTDAS shell code, see :ref:`installing <installing>`


