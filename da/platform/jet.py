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
# jet.py

"""
Author : peters 

Revision History:
File created on 06 Sep 2010.

"""
import os
import logging
import subprocess

from da.baseclasses.platform import Platform

std_joboptions = {'jobname':'test', 'jobaccount':'co2', 'jobnodes':'nserial 1', 'jobshell':'/bin/sh', 'depends':'', 'jobtime':'00:30:00', 'joblog':os.getcwd()}

class JetPlatform(Platform):
    def __init__(self):
        self.ID = 'NOAA jet'    # the identifier gives the platform name
        self.version = '1.0'     # the platform version used

        logging.debug('%s platform object initialized' % self.ID)
        logging.debug('%s version: %s' % (self.ID, self.version))


    def get_job_template(self, joboptions={}, block=False):
        """ Return the job template for a given computing system, and fill it with options from the dictionary provided as argument"""

        template = """#$ -N jobname \n""" + \
                   """#$ -A jobaccount \n""" + \
                   """#$ -pe jobnodes \n""" + \
                   """#$ -l h_rt=jobtime \n""" + \
                   """#$ -S jobshell \n""" + \
                   """#$ -o joblog \n""" + \
                   """#$ -cwd\n""" + \
                   """#$ -r n\n""" + \
                   """#$ -V\n""" + \
                   """#$ -j y\n"""

        if 'depends' in joboptions:
            template += """#$ -hold_jid depends \n"""

        if block:
            template += """#$ -sync y\n"""

        # First replace from passed dictionary
        for k, v in joboptions.items():
            while k in template:
                template = template.replace(k, v)

        # Fill remaining values with std_options
        for k, v in std_joboptions.items():
            while k in template:
                template = template.replace(k, v)

        return template

    def get_my_id(self):
        try:
            return os.environ['JOB_ID']
        except:
            return os.getpid()

    def submit_job(self, jobfile, joblog=None, block=False): 
        """ This method submits a jobfile to the queue, and returns the queue ID """
      
        cmd = ["qsub", jobfile]
        logging.info("A new task will be started (%s)" % cmd)
        output = subprocess.Popen(cmd, stdout=subprocess.PIPE).communicate()[0]
        logging.info(output)
        #jobid = output.split()[2]
        retcode = output.split()[-1]

        return retcode

    def kill_job(self, jobid):                   
        """ This method kills a running job """        
        output = subprocess.Popen(['qdel', jobid], stdout=subprocess.PIPE).communicate()[0]
        logging.info(output)
        return output

    def job_stat(self, jobid):                   
        """ This method gets the status of a running job """
        
        #output = subprocess.Popen(['sgestat'], stdout=subprocess.PIPE).communicate()[0]  ; logging.info(output)
        return ''

if __name__ == "__main__":
    pass
