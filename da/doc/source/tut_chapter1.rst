.. _tut_chapter1:

Chapter 1: Configuring CTDAS to run an experiment
----------------------------------------------------

**In which you will learn to modify the three primary files that together control a CTDAS experiment**

Step 1: Copy the required job files and rc-files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In order to run CTDAS, you'll need: 

    1. A job script written in bash shell that works on your computer
    2. A control script written in python that initializes the needed objects
    3. Three rc-files to initialize differen python objects

Luckily, there are pre-made examples avaliable for each of these.

You can grab the example job script in ``da/examples/das.jb`` and copy it to your main directory (the one you used to check out CTDAS). This is an example of its contents: ::

    #$ das.jb 
    #$ /bin/sh 

     echo "All output piped to file das.out"
     module load python
     python das.py rc=da.rc $1  >& das.out

This job script simply prepares your environment and then starts a python executable. Note that the job script can have extra headers that allow you to submit this script to the queue on your computer, if relevant. Configuring this script is part of Step 5 in this chapter.

Next, grab the example python script in ``/da/examples/das.py`` and copy it to the same location. Configuring this script is part of Step 3 in this chapter.

Finally, grab the two example rc-files ``da/rc/da.rc`` and ``da/rc/carbontracker.rc`` and also copy them. We will modify these first. The third rc-file needed is actually the ``tm5-ctdas.rc`` file you created in :ref:`Chapter 0 <tut_chapter0>` so that part is done!

.. note:: The 4 files above must be located in the main directory of the CTDAS tree, i.e., ``${yourdir}/da/ct/trunk/``, only the TM5 rc-file can live somewhere else.

Step 2: Modify the two primary rc-files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Locations and settings that control the inner workings of the CTDAS system are specified in: 

   1. The ``da.rc`` file, which describes your CTDAS configuration with respect to experiment name, time period, and lag
   2. The ``carbontracker.rc`` file, which describes your CTDAS configuration with respect to observations and statevector

You can open these in any text editor and replace the values of each key with appropriate settings for your experiment. For
example: ::

    dir.da_run          : ${HOME}/tmp/test_da

can be replaced by: ::

    dir.da_run          : /scratch/${USER}/my_first_ctdas_run

Which, as you have likely guessed, will change the location where CTDAS creates a directory structure and places input/output files 
for your simulation. See :mod:`initexit` for more information on these settings. It is especially important to set the keys ::

    da.system.rc        : carbontracker.rc                   ! the settings needed in your inversion system
    da.obsoperator.rc   : ${HOME}/Modeling/TM5/tm5-ctdas.rc  ! the rc-file needed to run youobservation operator

correctly. The first one refers to the rc-file (2) described above, while the second one refers to the rc-file you used to compile the TM5 model in :ref:`Chapter 0 <tut_chapter0>`. 

.. note:: Files and paths specified in the two basic rc-files must exist, or the system will fail and alert you to the
         fact that they are missing.

Where the ``da.rc`` file is rather self-explanatory, the ``carbontracker.rc`` file has keys that refer to the inner workings of CTDAS as described in :mod:`dasystem`.

.. note:: The example files are found in the da/rc/ directory of your CTDAS tree. You are encouraged 
         to always create copies of these primary rc-files before modifying them. The rc filenames are specified to the
         system before running CTDAS and thus you can use a different copy of these 
         files for different experiments and purposes. You can even create a 
         sub-directory with the settings of all your experiments if you like.

With the three rc-files now in-place and modified, we'll continue to modify the python control script.

Step 3: Modify the python control script
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Open the das.py script (or whatever you called it) in an editor and take a look at the code. The python control script first initializes 
some python objects needed to log activity and to parse command line arguments, and is then followed by 
a block where all the modules that are needed in your experiment are imported. The example below shows the import of several classes that 
are needed to run the CarbonTracker CO2 system on a computer referred to as MaunaLoa, using TM5 as 
an :ref:`observationoperator <Observation Operator>` ::

    ###########################################################################################
    ### IMPORT THE APPLICATION SPECIFIC MODULES HERE, TO BE PASSED INTO THE MAIN PIPELINE!!! ##
    ###########################################################################################

    from da.platform.maunaloa import MaunaloaPlatForm 
    from da.ct.dasystem import CtDaSystem 
    from da.ct.statevector import CtStateVector 
    from da.ct.obs import CtObservations 
    from da.tm5.observationoperator import TM5ObservationOperator 
    from da.ct.optimizer import CtOptimizer
    
Once the classes are loaded successfully, the objects are created. ::

    PlatForm    = MaunaloaPlatForm()
    DaSystem    = CtDaSystem(DaCycle['da.system.rc'])
    ObsOperator = TM5ObservationOperator(DaCycle['da.obsoperator.rc'])
    Samples     = CtObservations()
    StateVector = CtStateVector()
    Optimizer   = CtOptimizer()

.. note:: See how the initilization of the DaSystem and ObservationOperator object make use of the keys specified in your primary rc-file ! 

Modification of these objects might be desirable for more advanced users, and in case of the :ref:`platform <Platform>` object, even 
necessary (see next section). Once the objects are created, they are simply passed to a pipeline for the CTDAS. In the first 
chapter of the tutorial, we will assume this pipeline is immutable.

The only thing you might want to alter for now is the initialization of the PlatForm object, which is computer specific. How to create your own PlatForm object is described next. After completing this task, make sure you import this object in the ``das.py`` script and initialize it similar to the example ::

    from da.platform.<yourplatform> import <yourplatform> 
    PlatForm    = <yourplatform>()

Step 4: Creating a PlatForm object for your system
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

From the description of the :ref:`platform <Platform>` object, you will understand that this object is (partly) unique for each user, 
or at least for the computing environment of each user. Information on the computing system is therefore coded into a specific python object. 

.. warning:: **This object will need to be created for your system by you.**

Luckily, part of the work is already done. In the ``da/baseclasses`` subdirectory you will find a baseclass :ref:`platform <Platform>` 
which serves as a blueprint for your own ``Platform`` object. This is done through class inheritance. As an example, you can open one 
of the files ``jet.py`` or ``maunaloa.py`` in the ``da/platform`` directory alongside the original ``da/baseclasses/platform.py``. One
of the first things to notice is the headers of the ``class PlatForm`` in the baseclass:  ::
    
    class PlatForm(object):

and in the derived class::

    from da.baseclasses.platform import PlatForm
    class PlatForm(PlatForm):

This tells you that the second object actually starts as a copy of the baseclass. But then, we see that the derived class has
a new implementation of the method ``GetJobTemplate`` from which the first set of lines are below: ::

    def GetJobTemplate(self,joboptions={},block=False):
        """ 
            Return the job template for a given computing system, 
            and fill it with options from the dictionary provided as argument
        """

        template = """#$ -N jobname \n"""+ \
                   """#$ -A jobaccount \n"""+ \
                   """#$ -pe jobnodes \n"""+ \
                   """#$ -l h_rt=jobtime \n"""+ \
                   """#$ -S jobshell \n"""+ \
                   """#$ -o joblog \n"""+ \
                   """#$ -cwd\n"""+ \
                   """#$ -r n\n"""+ \
                   """#$ -V\n"""+ \
                   """#$ -j y\n"""

While the baseclass did not have any functionality, a call to the ``GetJobTemplate`` method of the derived class will actually return a template for 
a job script on NOAA's "jet" supercomputer, so that we can submit jobs to its queue. **By modifying each of the methods in your own 
derived PlatForm class in the same way, you can make each method work on your system**. 

Once you have created your own PlatForm object, and you have successfully imported and instantiated it in your primary python run script, 
you are ready for the last step.

.. note:: Sometimes it is faster and easier to test your newly created class 'offline'. At the end of your module, following the __main__ 
         section you can add lines to test your PlatForm object before plugging it into the CTDAS. 

Step 5: Modifying the job script
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As a final step, open your job script (1) again and see whether it calls the right python control script, and whether it has the right rc filename specified. If so, you are ready for :ref:`Chapter 2: Running your first experiment <tut_chapter2>`.




