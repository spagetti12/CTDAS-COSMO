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
# tm5_tools.py

"""
Author : peters 

Revision History:
File created on 09 Feb 2009.
Major modifications to go to a class-based approach, July 2010.

This module holds specific functions needed to use the TM5 model within the data assimilation shell. It uses the information 
from the DA system in combination with the generic tm5.rc files. 

The TM5 model is now controlled by a python subprocess. This subprocess consists of an MPI wrapper (written in C) that spawns
a large number ( N= nmembers) of TM5 model instances under mpirun, and waits for them all to finish.

The design of the system assumes that the tm5.x (executable) was pre-compiled with the normal TM5 tools, and is residing in a 
directory specified by the ${RUNDIR} of a tm5 rc-file. This tm5 rc-file name is taken from the data assimilation rc-file. Thus,
this python shell does *not* compile the TM5 model for you!

"""

import os
import sys
import logging
import shutil
import datetime
import subprocess
from string import join
import glob
sys.path.append(os.getcwd())
sys.path.append("../../")

import da.tools.rc as rc
from da.tools.general import create_dirs, to_datetime
from da.baseclasses.observationoperator import ObservationOperator

identifier = 'TM5'
version = 'release 3.0'
mpi_shell_filename = 'tm5_mpi_wrapper'
mpi_shell_location = 'da/bin/'


################### Begin Class TM5 ###################



class TM5ObservationOperator(ObservationOperator):
    """ This class holds methods and variables that are needed to run the TM5 model. It is initiated with as only argument a TM5 rc-file
        location. This rc-file will be used to figure out the settings for the run. 
        
        *** This method of running TM5 assumes that a pre-compiled tm5.exe is present, and it will be run from time.start to time.final ***
        
        These settings can be modified later. To run a model version, simply compile the model using an existing TM5 rc-file, then 
        open python, and type:

           []> tm=TM5('/Users/peters/Modeling/TM5/tutorial.rc')
           []> tm.write_rc()
           []> tm.WriteRunRc()
           []> tm.run()

        To use this class inside a data assimilation cycle, a stand-alone method "setup()" is included which modifies the TM5
        settings according to an external dictionary of values to overwrite, and then runs the TM5 model.

    
    """

    def __init__(self, filename):
        """ The instance of an TMObservationOperator is application dependent """
        self.ID = identifier    # the identifier gives the model name
        self.version = version       # the model version used
        self.restart_filelist = []
        self.output_filelist = []
        
        self.outputdir = None # Needed for opening the samples.nc files created 

        self.load_rc(filename)   # load the specified rc-file
        self.validate_rc()         # validate the contents

        logging.info('Observation Operator initialized: %s (%s)' % (self.ID, self.version))
        
    def setup(self, dacycle):
        """ 
           Execute all steps needed to prepare the ObsOperator for use inside CTDAS, only done at the very first cycle normally
    
        """
        self.dacycle = dacycle
        
        if self.dacycle['time.restart'] == False or self.dacycle['transition'] == True:
            newitemsmeteo = {}
            if self.dacycle['transition']:
                if not (self.tm_settings[self.mysourcedirs]).find('ml91') >= 0:
                    newdirs = self.tm_settings[self.mysourcedirs] + self.addedsourcedirs
                else: newdirs = self.tm_settings[self.mysourcedirs]
                # Modify the rc-file to reflect directory structure defined by CTDAS
                # and change items for transition from 25 to 34 layers for od meteo in dec 2005
                newitems = {'my.basedir' : self.dacycle['dir.exec'],
                            'istart' : self.transitionvalue,
                             'my.source.dirs' : newdirs,
                             'build.copy.dirs' : newdirs}
                for k, v in self.tm_settings.items():
                    if type(v) is str:
                        if v.find('ml60') > -1:
                            newitemsmeteo[k] = v.replace('ml60','ml91')
                        if v.find('tropo25') > -1:
                            newitemsmeteo[k] = v.replace('tropo25','tropo34')
                        if v.find('ml25') > -1:
                            newitemsmeteo[k] = v.replace('ml25','tropo34')
                        if v.find('od_L60') > -1:
                            newitemsmeteo[k] = 'od_L91'
                        if k == 'my.meteo.nlev':
                            newitemsmeteo[k] = '91'
                logging.debug('Resetting TM5 to perform transition of od meteo from 25 to 34 levels')
            else:        
                logging.info('First time step, setting up and compiling the TM5 model before proceeding!')
                # Modify the rc-file to reflect directory structure defined by CTDAS
                newitems = {'my.basedir': self.dacycle['dir.exec']}

            self.modify_rc(newitems)
            self.modify_rc(newitemsmeteo)

            # Create the TM5 run directory to hold a copy of the modified rc-file

            tm5compiledir = self.tm_settings[self.rundirkey]
            create_dirs(tm5compiledir)

            rcfilename = os.path.join(tm5compiledir, 'tm5_setup_init.rc')
            self.write_rc(rcfilename)

            # Compile TM5 in the new directory, but only if this is a fresh start
            logging.debug('Original rc file: %s '%(self.dacycle['da.obsoperator.rc']))
            self.compile_tm5(rcfilename)
            
            newrcfilename = os.path.join(self.tm_settings['rundir'], self.tm_settings['install.rc'])

            #Use a TM5 restart file in the first cycle (instead of init file). Used now for the CO project.
            if 'da.obsoperator.restartfileinfirstcycle' in self.dacycle:
                restartfilename = self.dacycle['da.obsoperator.restartfileinfirstcycle']
                targetdir = self.tm_settings[self.savedirkey]
                create_dirs(targetdir)
                for file in glob.glob(restartfilename):
                    fname = os.path.split(file)[1]
                    logging.debug('Copied TM5 restart file to TM5 restart directory for first cycle: %s'%fname)
                    shutil.copy(file,os.path.join(targetdir,fname))

            # Replace the rc filename for TM5 with the newly created one in the new run directory

            
            logging.debug('Working copy of the tm5.rc file is in place (%s) ' % newrcfilename)

            # And also replace the path to the ObservationOperator in the dacycle object so we can work from the TM5 copy from here on

            self.dacycle['da.obsoperator.rc'] = newrcfilename

            logging.debug('...and set as the da.obsoperator.rc value in this dacycle ')


        logging.debug('Reloading the da.obsoperator.rc file for this dacycle')
        self.load_rc(self.dacycle['da.obsoperator.rc'])
        logging.debug('Note that the obsoperator is not recompiled if this is a recovery from a crash!!!')
        

    def compile_tm5(self, rcfilename):
        """
            Compile TM5 model using setup_tm5 and the modified rc-file
        """
        if 'da.obsoperator.home' in self.dacycle:
            os.chdir(self.dacycle['da.obsoperator.home'])
            logging.debug('Proceeding from TM5 root dir (%s)'% self.dacycle['da.obsoperator.home'])
        else:
            tm5_dir = os.path.split(self.dacycle['da.obsoperator.rc'])[0]
            logging.warning('Guessing your TM5 root dir from the rc filename')
            logging.warning('Try adding a key da.obsoperator.home to your da.rc')
            logging.warning('Proceeding from guessed TM5 root dir (%s) ' % tm5_dir)

            os.chdir(tm5_dir)
  
        if self.dacycle['transition']:
            cmd = ['python', 'setup_tm5', '-n', '--%s' % self.dacycle.daplatform.give_queue_type(), rcfilename]
        else: cmd = ['python', 'setup_tm5', '--%s' % self.dacycle.daplatform.give_queue_type(), rcfilename]

        logging.info('Starting the external TM5 setup script')
        logging.info('using command ... %s' % join(cmd))

        retcode = subprocess.call(cmd)
        os.chdir(self.dacycle['dir.da_submit'])

        if retcode != 0:
            logging.error('Compilation failed, quitting CTDAS')
            raise IOError
            sys.exit(2)
        else:
            logging.info('Compilation successful, continuing')

    def prepare_run(self):
        """ 
        Prepare a forward model TM5 run, this consists of:

          - reading the working copy TM5 rc-file, 
          - validating it, 
          - modifying the values,
          - Removing the existing tm5.ok file if present
    
        """

# Write a modified TM5 model rc-file in which run/break times are defined by our da system

        new_items = {
                    'submit.options': self.dacycle.daplatform.give_blocking_flag(),
                    self.timestartkey: self.dacycle['time.sample.start'],
                    self.timefinalkey: self.dacycle['time.sample.end'],
                    'jobstep.timerange.start': self.dacycle['time.sample.start'],
                    'jobstep.timerange.end': self.dacycle['time.sample.end'],
                    'jobstep.length': 'inf',
                    'ct.params.input.dir': self.dacycle['dir.input'],
                    'ct.params.input.file': os.path.join(self.dacycle['dir.input'], 'parameters'),
                    'output.flask.infile': self.dacycle['ObsOperator.inputfile'] ,
                    'output.flask': 'True'
                    }

        if self.dacycle['transition']:
            new_items[self.istartkey] = self.transitionvalue
            logging.debug('Resetting TM5 to perform transition of od meteo from 25 to 34 levels') 
        elif self.dacycle['time.restart']:  # If this is a restart from a previous cycle, the TM5 model should do a restart
            new_items[self.istartkey] = self.restartvalue
            logging.debug('Resetting TM5 to perform restart')
        else:
            if 'da.obsoperator.restartfileinfirstcycle' not in self.dacycle:
                new_items[self.istartkey] = self.coldstartvalue  # if not, start TM5 'cold'
                logging.debug('Resetting TM5 to perform cold start')
            else:
                new_items[self.istartkey] = self.restartvalue  # If restart file is specified, start TM5 with initial restartfile
                logging.debug('Resetting TM5 to start with restart file: %s'%self.dacycle['da.obsoperator.restartfileinfirstcycle'])

        if self.dacycle['time.sample.window'] != 0:  # If this is a restart from a previous time step within the filter lag, the TM5 model should do a restart
            new_items[self.istartkey] = self.restartvalue
            logging.debug('Resetting TM5 to perform restart')
        
        # If neither one is true, simply take the istart value from the tm5.rc file that was read

        self.modify_rc(new_items)
        self.write_rc(self.rc_filename)

	# Define the name of the file that will contain the modeled output of each observation

    	self.simulated_file = os.path.join(self.outputdir, 'flask_output.%s.nc' % self.dacycle['time.sample.stamp'])

    def load_rc(self, name):
        """ 
        This method loads a TM5 rc-file with settings for this simulation 
        """
        self.rcfile = rc.RcFile(name)
        self.tm_settings = self.rcfile.values
        self.rc_filename = name

        if 'my.source.dirs' in self.tm_settings.keys():
            self.rcfiletype = 'pycasso'
        else:
            self.rcfiletype = 'pre-pycasso'
        logging.debug('TM5 rc-file loaded successfully')

    def validate_rc(self):
        """
        Validate the contents of the tm_settings dictionary and add extra values. The required items for the TM5 rc-file
        are specified in the tm5_tools module, as dictionary variable "needed_rc_items".
        """

        if self.rcfiletype == 'pycasso':
            self.projectkey = 'my.project.dir'
            self.rundirkey = 'my.run.dir'
            self.outputdirkey = 'output.dir'
            self.savedirkey = 'restart.write.dir'
            self.timestartkey = 'timerange.start'
            self.timefinalkey = 'timerange.end'
            self.timelengthkey = 'jobstep.length'
            self.istartkey = 'istart'
            self.restartvalue = 33
            self.coldstartvalue = 9
            self.transitionvalue = 61
            self.mysourcedirs = 'my.source.dirs'
            self.addedsourcedirs = '  proj/levels/ml91/trunk proj/levels/ml91/tropo34/trunk'
        else:
            self.projectkey = 'runid'
            self.rundirkey = 'rundir'
            self.outputdirkey = 'outputdir'
            self.savedirkey = 'savedir'
            self.timestartkey = 'time.start'
            self.timefinalkey = 'time.final'
            self.timelengthkey = 'time.break.nday'
            self.istartkey = 'istart'
            self.restartvalue = 3
            self.coldstartvalue = 9
            self.transitionvalue = 61

        needed_rc_items = [
                            self.projectkey,
                            self.rundirkey,
                            self.outputdirkey,
                            self.savedirkey,
                            self.timestartkey,
                            self.timefinalkey,
                            self.timelengthkey,
                            self.istartkey
                          ]

        for k, v in self.tm_settings.items():
            if v == 'True' : self.tm_settings[k] = True
            if v == 'False': self.tm_settings[k] = False
            if 'date' in k : self.tm_settings[k] = to_datetime(v)
            if 'time.start' in k : 
                self.tm_settings[k] = to_datetime(v, fmt='TM5')
            if 'time.final' in k : 
                self.tm_settings[k] = to_datetime(v, fmt='TM5')
            if 'timerange.start' in k : 
                self.tm_settings[k] = to_datetime(v)
            if 'timerange.end' in k : 
                self.tm_settings[k] = to_datetime(v)

        for key in needed_rc_items:
            if key not in self.tm_settings:
                msg = 'Missing a required value in rc-file : %s' % key
                logging.error(msg)
                raise IOError(msg)
        logging.debug('rc-file has been validated succesfully')


    def modify_rc(self, newvalues):
        """ 
        Modify parts of the tm5 settings, for instance to give control of file locations to the DA shell
        instead of to the tm5.rc script. 

        Note that we replace these values in all {key,value} pairs of the tm5.rc file!

        """
    
        for k, v in newvalues.items():
            if key in self.tm_settings:
                # keep previous value
                v_orig = self.tm_settings[k]
                #replace with new
                self.tm_settings[k] = v
                #replace all instances of old with new, but only if it concerns a name of a path!!!
                if os.path.exists(str(v)): 
                    for k_old, v_old in self.tm_settings.items():
                        if not isinstance(v_old, str): 
                            continue
                        if str(v_orig) in str(v_old): 
                            v_new = str(v_old).replace(str(v_orig), str(v))
                            self.tm_settings[k_old] = v_new

                logging.debug('Replaced tm5 rc-item %s : %s' % (k,v))

            else:
                self.tm_settings[k] = v
                logging.debug('Added new tm5 rc-item %s : %s' % (k,v))


    def write_rc(self, tm5rcfilename):
        """ 
        Write the rc-file settings to a tm5.rc file in the rundir
        """
        rc.write(tm5rcfilename, self.tm_settings)
        logging.debug("Modified rc file for TM5 written (%s)" % tm5rcfilename)

    def validate_input(self):
        """
        Make sure that parameter files are written to the TM5 inputdir, and that observation lists are present
        """

        datadir = self.tm_settings['ct.params.input.dir']
        if not os.path.exists(datadir):
            msg = "The specified input directory for the TM5 model to read from does not exist (%s), exiting..." % datadir 
            logging.error(msg)
            raise IOError(msg)

        datafiles = os.listdir(datadir)

        obsfile = self.dacycle['ObsOperator.inputfile']

        if not os.path.exists(obsfile):
            msg = "The specified obs input file for the TM5 model to read from does not exist (%s), exiting..." % obsfile  
            logging.error(msg)
            if 'forward.savestate.dir' not in self.dacycle:
                raise IOError(msg)

        for n in range(int(self.dacycle['da.optimizer.nmembers'])):
            paramfile = 'parameters.%03d.nc' % n
            if paramfile not in datafiles:
                msg = "The specified parameter input file for the TM5 model to read from does not exist (%s), exiting..." % paramfile 
                logging.error(msg)
                raise IOError(msg)

        # Next, make sure there is an actual model version compiled and ready to execute

        targetdir = os.path.join(self.tm_settings[self.rundirkey])

        if self.rcfiletype == 'pycasso':
            self.tm5_exec = os.path.join(targetdir, self.tm_settings['my.basename'] + '.x')
        else:
            self.tm5_exec = os.path.join(targetdir, 'tm5.x')

        if not os.path.exists(self.tm5_exec):
            logging.error("Required TM5 executable was not found %s" % self.tm5_exec)
            logging.error("Please compile the model with the specified rc-file and the regular TM5 scripts first")
            raise IOError


    def get_initial_data(self):
        """ This method places all initial data needed by an ObservationOperator in the proper folder for the model.
            For TM5, this means copying the save_*.hdf* files to the dir.save directory from which TM5 will read initial
            concentrations for all tracers. 

            We get the input data from the restart.current directory at 2 times:
                (1) When the model starts the forecast over nlag cycles
                (2) When the model starts the advance step over 1 cycle


         """
        logging.debug("Moving TM5 model restart data from the restart directory to the TM5 save dir")

        # First get the restart data for TM5 from the current restart dir of the filter

        sourcedir = self.dacycle['dir.restart']
        targetdir = self.tm_settings[self.savedirkey]
        self.outputdir = self.tm_settings[self.outputdirkey]  # Needed further downstream to collect output data from TM5

        filterlist = '%s' % self.dacycle['time.start'].strftime('%Y%m%d')

        for f in os.listdir(sourcedir):
            fpath = os.path.join(sourcedir, f)
            if os.path.isdir(fpath): # skip dirs
                logging.debug("           [skip] .... %s " % fpath)
                continue    
            #if not f.startswith('save_'):
            if not f.startswith('TM5_restart'):
                logging.debug("           [skip] .... %s " % fpath)
                continue    
            if not filterlist in f:
                logging.debug("           [skip] .... %s " % fpath)
                continue    

            # all okay, copy file

            logging.debug("           [copy] .... %s " % fpath)
            shutil.copy(fpath, fpath.replace(sourcedir, targetdir))
        logging.debug("All restart data have been copied from the restart/current directory to the TM5 save dir")

    def run_forecast_model(self):
        self.prepare_run()
        self.validate_input()
        self.run()
        self.save_data()

    def run(self):
        """ 
         Start the TM5 executable. A new log file is started for the TM5 model IO, and then a subprocess is
         spawned with the tm5_mpi_wrapper and the tm5.x executable. The exit code of the model is caught and
         only if successfull on all processors will execution of the shell continue. 
         
        """
        cwd = os.getcwd()

        # From here on, several options should be implemented.

        #
        # (1) Where an mpi process is forked to do a TM5 instance with N tracers, each an ensemble member
        #
        # (2) Where N processes are spawned, each being one TM5 instance representing one member
        #
        # (3) Where N/m processes are spawned, each being a TM5 instance that handles m ensemble members
        #
        # In principle, it is best to make these processes produce scripts that can be executed stand-alone, or
        # be submitted to a queue.
        #

        # Open logfile and spawn model, wait for finish and return code

        # Code for Option (1)

        code = self.tm5_with_n_tracers()

        if code == 0:
            logging.info('Finished model executable succesfully (%s)' % code)
            self.Status = 'Success'
        else:
            logging.error('Error in model executable return code: %s ' % code)
            self.Status = 'Failed'
            raise OSError

        # Return to working directory

        os.chdir(cwd)

        return code

    def tm5_with_n_tracers(self):
        """ Method handles the case where one TM5 model instance with N tracers does the sampling of all ensemble members"""


        tm5submitdir = os.path.join(self.tm_settings[self.rundirkey])
        logging.info('tm5submitdir', tm5submitdir)

        # Go to executable directory and start the subprocess, using a new logfile

        os.chdir(tm5submitdir)
        logging.debug('Changing directory to %s ' % tm5submitdir)

        # Remove the tm5.ok file from a previous run, placed back only if a successful TM5 run is executed

        okfile = 'tm5.ok'
        if os.path.exists(okfile): 
            os.remove(okfile)

        # Prepare a job for the current platform, this job needs to account for all job parameters such as
        # runtime, queue request, number of processors, and other platform specific parameters

        # It is easiest if all of these options are 'pre-configured' through the tm5.rc file that is used, then
        # running the ObsOperator on all platforms simply becomes a matter of running the ./setup_tm5 script
        # that was also used to compile TM5 when not actually running CTDAS. The first run of CTDAS
        # then could be one where TM5 is actually compiled and run, the rest of the jobs would then *not* re-compile.
        # An option would need to be added to force a re-compile of the TM5 code, for debugging purposes.

        # file ID and names
        submitcommand = self.tm_settings['submit.command'] 
        logging.info('Submitting job at %s' % datetime.datetime.now())
        code = subprocess.call(submitcommand.split()) 
        logging.info('Resuming job at %s' % datetime.datetime.now())
        
        if not os.path.exists(okfile): 
            code = -1
        else:
            code = 0

        return code

    def save_data(self):
        """ Copy the TM5 recovery data from the outputdir to the TM5 savedir, also add the restart files to a list of names
            that is used by the dacycle object to collect restart data for the filter.

            WP Note: with the new pycasso restart files we no longer need to copy save files from outdir to savedir

            Note 2: also adding the weekly mean flux output to the output_filelist for later collection
         """

        sourcedir = os.path.join(self.tm_settings[self.savedirkey])
        filterlist = ['%s' % self.tm_settings[self.timefinalkey].strftime('%Y%m%d')]

        logging.debug("Creating a new list of TM5 restart data")
        logging.debug("           from directory: %s " % sourcedir)
        logging.debug("           with filter: %s " % filterlist)


        # Start from empty lists for each TM5 run. Note that these "private" lists from the obs operator are later on appended to the system
        # lists

        self.restart_filelist = []

        for fil in os.listdir(sourcedir):
            fil = os.path.join(sourcedir, fil)
            if os.path.isdir(fil): # skip dirs
                skip = True
            elif filterlist == []:      # copy all
                skip = False        
            else:                   # check filter
                skip = True         # default skip
                for f in filterlist:
                    if f in fil: 
                        skip = False # unless in filterlist
                        break
                
            if skip: 
                logging.debug("           [skip] .... %s " % fil)
                continue    

            self.restart_filelist.append(fil)
            logging.debug("           [added to restart list] .... %s " % fil)

        sourcedir = os.path.join(self.tm_settings[self.outputdirkey])
        sd_ed = self.dacycle['time.sample.stamp']
        filterlist = ['flask_output.%s' % sd_ed, 'flux1x1_%s' % sd_ed]

        logging.debug("Creating a new list of TM5 output data to collect")
        logging.debug("           from directory: %s " % sourcedir)
        logging.debug("           with filter: %s " % filterlist)


        # Start from empty lists for each TM5 run. Note that these "private" lists from the obs operator are later on appended to the system
        # lists

        self.output_filelist = []

        for fil in os.listdir(sourcedir):
            fil = os.path.join(sourcedir, fil)

            if os.path.isdir(fil): # skip dirs
                skip = True
            elif filterlist == []:      # copy all
                skip = False        
            else:                   # check filterlist
                skip = True         # default skip
                for f in filterlist:
                    if f in fil: 
                        skip = False # unless in filterlist
                        break
                
            if skip: 
                logging.debug("           [skip] .... %s " % fil)
                continue    

            self.output_filelist.append(fil)
            logging.debug("           [added to output list] .... %s " % fil)


################### End Class TM5 ###################
                   

if __name__ == "__main__":
    pass




