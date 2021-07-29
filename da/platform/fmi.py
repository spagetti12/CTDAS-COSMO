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
# fmi.py

"""
Author : Aki 

Revision History:
File created on 25 Jun 2013.

"""

import logging
import subprocess

from da.baseclasses.platform import Platform
std_joboptions = {'jobaccount':'ch4'}

class FmiPlatform(Platform):
    def __init__(self):
        self.ID     = 'FMI voima'    # the identifier gives the platform name
        self.version = ''     # the platform version used

    #def give_blocking_flag(self):
    #    """
    #    Returns a blocking flag, which is important if tm5 is submitted in a queue system. The python ctdas code is forced to wait before tm5 run is finished
    #    """
    #    return "-s"

    #def give_queue_type(self):
    #    """
    #    Return a queue type depending whether your computer system has a queue system, or whether you prefer to run in the foreground. 
    #    On most large systems using the queue is mandatory if you run a large job.
    #    """
    #    return "queue"

    def get_job_template(self, joboptions={}, block=False):
        """ 
        Returns the job template for a given computing system, and fill it with options from the dictionary provided as argument.
        The job template should return the preamble of a job that can be submitted to a queue on your platform, 
        examples of popular queuing systems are:
            - SGE
            - MOAB
            - XGrid
            -

        A list of job options can be passed through a dictionary, which are then filled in on the proper line,
        an example is for instance passing the dictionary {'account':'co2'} which will be placed 
        after the ``-A`` flag in a ``qsub`` environment.

        An extra option ``block`` has been added that allows the job template to be configured to block the current
        job until the submitted job in this template has been completed fully.
        """

        template = """#!/bin/bash"""+ \
                   """ \n"""+ \
                   """ \n"""+ \
                   """#PBS -N jobname \n"""+ \
                   """#PBS -j oe \n"""+ \
                   """#PBS -l walltime=jobtime \n"""+ \
                   """#PBS -l mppwidth=jobpes \n"""+ \
                   """#PBS -l mppnppn=jobnodes \n"""+ \
                   """\n"""

        if 'depends' in joboptions:
            template += """#$ -hold_jid depends \n"""

        # First replace from passed dictionary
        for k,v in joboptions.items():
            while k in template:
                template = template.replace(k,v)

        # Fill remaining values with std_options
        for k,v in std_joboptions.items():
            while k in template:
                template = template.replace(k,v)

        return template


        #msg1  = 'Platform initialized: %s'%self.Identifier      ; logging.info(msg1)
        #msg2  = '%s version: %s'%(self.Identifier,self.Version) ; logging.info(msg2)


    def submit_job(self, jobfile, joblog=None, block=False):
        """ This method submits a jobfile to the queue, and returns the queue ID """

        if block:
            cmd = ["qsub", "-s", jobfile]
            logging.info("A new task will be started (%s)" % cmd)
            output = subprocess.Popen(cmd, stdout=subprocess.PIPE).communicate()[0]
            logging.info(output)
            print 'output', output
            jobid = output.split()[0]             
            print 'jobid', jobid
        else:
            cmd = ["qsub", jobfile]
            logging.info("A new task will be started (%s)" % cmd)
            output = subprocess.Popen(cmd, stdout=subprocess.PIPE).communicate()[0]  ; logging.info(output)
            jobid = output.split()[0]

        return jobid


if __name__ == "__main__":
    pass
