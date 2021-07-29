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
# da_initexit.py

"""
modified for module initext
Revision History:
File created on Apr. 2016 by Aki
File revised on Feb. 2017 by Aki

"""
import logging
import os
import sys
import numpy as np
from string import join

import da.tools.rc as rc
from da.tools.initexit import *

needed_da_items = [
    'time.start',
    'time.finish',
    'time.nlag',
    'time.cycle',
    'dir.da_run',
    'da.system',
    'da.system.rc',
    'da.obsoperator',
    'da.obsoperator.rc',
    'da.optimizer.nmembers',
    'da.suffixname']


#def validate_rc2(self):
#    """ 
#    Validate the contents of the rc-file given a dictionary of required keys. 
#    Currently required keys are :attr:`~da.tools.initexit.needed_da_items`
#    """
#
#    for k, v in self.iteritems():
#        if v in ['True', 'true', 't', 'T', 'y', 'yes']:
#            self[k] = True
#        if v in ['False', 'false', 'f', 'F', 'n', 'no']:
#            self[k] = False
#        if 'date' in k : 
#            self[k] = to_datetime(v)
#        if k in ['time.start', 'time.end', 'time.finish', 'da.restart.tstamp']:
#            self[k] = to_datetime(v)
#    for key in needed_da_items:
#        if not self.has_key(key):
#            msg = 'Missing a required value in rc-file : %s' % key
#            logging.error(msg)
#            logging.error('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ')
#            logging.error('Please note the update on Dec 02 2011 where rc-file names for DaSystem and ')
#            logging.error('are from now on specified in the main rc-file (see da/rc/da.rc for example)')
#            logging.error('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ')
#            raise IOError, msg
#    logging.debug('DA Cycle settings have been validated succesfully')
#CycleControl.validate_rc = validate_rc2

def submit_next_cycle_fmi(self):
    """ 
    Submit the next job of a DA cycle, this consists of 
        * Changing to the working directory from which the job was started initially
        * create a line to start the master script again with a newly created rc-file
        * Submitting the jobfile 

    If the end of the cycle series is reached, no new job is submitted.

    """
    

    if self['time.end'] < self['time.finish']:

        # file ID and names
        jobid = self['time.end'].strftime('%Y%m%d') 
        targetdir = os.path.join(self['dir.exec'])
        jobfile = os.path.join(targetdir, 'jb.%s.jb' % jobid)
        logfile = os.path.join(targetdir, 'jb.%s.log' % jobid)
        # Template and commands for job
        jobparams = {'jobname':"ctdas", 'jobtime':'05:00:00', 'logfile': logfile, 'errfile': logfile,
                      'jobpes':'120','jobnodes':'20'}
        template = self.daplatform.get_job_template(jobparams)
        execcommand = os.path.join(self['dir.da_submit'], sys.argv[0]) 
        if '-t' in self.opts:
            (self.opts).remove('-t') 
        #template += 'python %s rc=%s %s' % (execcommand, self['da.restart.fname'], join(self.opts, '')) 
        template += 'cd %s \n' %self['dir.da_submit']
        template += '%s rc=%s %s' % (sys.argv[0], self['da.restart.fname'], join(self.opts, '')) 
        # write and submit 
        self.daplatform.write_job(jobfile, template, jobid)
        jobid = self.daplatform.submit_job(jobfile, joblog=logfile) 
        job_file = 'ctdas.o%s' %jobid[0:-4]
        logging.info('Jobfile %s' %job_file )
        logging.info('DA cycle has been submitted.')

    else:
        logging.info('Final date reached, no new cycle started')

def setup_file_structure_fmi(self):

    # Aki: name of file structure changed
    # At FMI, job files are stored in executed folders, and cannot define where to put them.
    # So job file folder is not created.

    filtertime = self['time.start'].strftime('%Y%m%d')

    suf = self['da.suffixname']
    self['dir.exec'] = os.path.join(self['dir.da_run'], 'exec')
    self['dir.input'] = os.path.join(self['dir.da_run'], 'input_%s' %suf )
    self['dir.output'] = os.path.join(self['dir.da_run'], 'output_%s' %suf )
    self['dir.analysis'] = os.path.join(self['dir.da_run'], 'analysis_%s' %suf )
    #self['dir.jobs'] = os.path.join(self['dir.da_run'], 'jobs')
    self['dir.restart'] = os.path.join(self['dir.da_run'], 'restart_%s' %suf )

    logging.info("setup_file_structure %s" %self['dir.output'])
    create_dirs(self['dir.da_run'])
    create_dirs(os.path.join(self['dir.exec']))
    create_dirs(os.path.join(self['dir.input']))
    create_dirs(os.path.join(self['dir.output']))
    create_dirs(os.path.join(self['dir.analysis']))
    #create_dirs(os.path.join(self['dir.jobs']))
    create_dirs(os.path.join(self['dir.restart']))

    logging.info('Succesfully created the file structure for the assimilation job')

def setup_fmi(self):

    if self['transition']:
        logging.info("Transition of filter from previous step with od meteo from 25 to 34 levels")
        self.setup_file_structure()
        strippedname = os.path.split(self['jobrcfilename'])[-1]
        self['jobrcfilename'] = os.path.join(self['dir.exec'], strippedname)
        self.read_random_seed(False)

    elif self['time.restart']:
        logging.info("Restarting filter from previous step")
        self.setup_file_structure()
        strippedname = os.path.split(self['jobrcfilename'])[-1]
        self['jobrcfilename'] = os.path.join(self['dir.exec'], strippedname)
        self.read_random_seed(False)

    else: #assume that it is a fresh start, change this condition to more specific if crash recover added
        logging.info("First time step in filter sequence")
        self.setup_file_structure()

        # expand jobrcfilename to include exec dir from now on.
        # First strip current leading path from filename

        strippedname = os.path.split(self['jobrcfilename'])[-1]
        self['jobrcfilename'] = os.path.join(self['dir.exec'], strippedname)
        if 'copyregionsfile' in self.dasystem:
          shutil.copy(os.path.join(self.dasystem['regionsfile']),os.path.join(self['dir.exec'],'da','methane','analysis','copied_regions.nc'))
          logging.info('Copied regions file to the analysis directory: %s'%os.path.join(self.dasystem['regionsfile'])) 

          if 'extendedregionsfile' in self.dasystem:
              shutil.copy(os.path.join(self.dasystem['extendedregionsfile']),os.path.join(self['dir.exec'],'da','analysis','copied_regions_extended.nc')) 
              logging.info('Copied extended regions file to the analysis directory: %s'%os.path.join(self.dasystem['extendedregionsfile'])) 

        for filename in glob.glob(os.path.join(self['dir.exec'],'da','analysis','*.pickle')):
            logging.info('Deleting pickle file %s to make sure the correct regions are used'%os.path.split(filename)[1])
            os.remove(filename) 

        for filename in glob.glob(os.path.join(self['dir.exec'],'*.pickle')):
            logging.info('Deleting pickle file %s to make sure the correct regions are used'%os.path.split(filename)[1])
            os.remove(filename) 

        if 'random.seed.init' in self.dasystem:
            self.read_random_seed(True)

    self.parse_times()
       #self.write_rc(self['jobrcfilename'])


CycleControl.submit_next_cycle = submit_next_cycle_fmi
CycleControl.setup_file_structure = setup_file_structure_fmi
CycleControl.setup = setup_fmi

if __name__ == "__main__":
    pass

