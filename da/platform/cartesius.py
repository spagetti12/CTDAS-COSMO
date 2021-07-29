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
# cartesius.py

"""
Author : peters 

Revision History:
File created on 06 Sep 2010.

"""

import logging
import subprocess

from da.baseclasses.platform import Platform

std_joboptions = {'jobname':'test', 'jobaccount':'co2', 'jobtype':'serial', 'jobshell':'/bin/sh', 'depends':'', 'jobtime':'24:00:00', 'jobinput':'/dev/null', 'jobnodes':'1', 'jobtasks':'', 'modulenetcdf':'netcdf/4.1.2', 'networkMPI':'','jobqueue': 'normal'}


class CartesiusPlatform(Platform):
    def __init__(self):
        self.ID = 'cartesius'    # the identifier gives the platform name
        self.version = '1.0'     # the platform version used


    def give_blocking_flag(self):
        """
        Returns a blocking flag, which is important if tm5 is submitted in a queue system. The python ctdas code is forced to wait before tm5 run is finished

	    -on Huygens: return "-s"
            -on Maunaloa: return "" (no queue available)
            -on Jet/Zeus: return
        """
        return ""

    def give_queue_type(self):
        """
        Return a queue type depending whether your computer system has a queue system, or whether you prefer to run in the foreground. 
        On most large systems using the queue is mandatory if you run a large job.
            -on Huygens: return "queue"
            -on Maunaloa: return "foreground" (no queue available)
            -on Jet/Zeus: return  

        """
        return "foreground"

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
        
        #template = """## \n"""+ \
        #           """## This is a set of dummy names, to be replaced by values from the dictionary \n"""+ \
        #           """## Please make your own platform specific template with your own keys and place it in a subfolder of the da package.\n """+ \
        #           """## \n"""+ \
        #           """ \n"""+ \
        #           """#$ jobname \n"""+ \
        #           """#$ jobaccount \n"""+ \
        #           """#$ jobnodes \n"""+ \
        #           """#$ jobtime \n"""+ \
        #           """#$ jobshell \n"""+ \
        #           """\n"""+ \
        #           """source /usr/bin/sh\n"""+ \
        #           """module load python\n"""+ \
        #           """\n"""

       
        template = """#!/bin/bash \n""" + \
                   """## \n""" + \
                   """## This is a set of dummy names, to be replaced by values from the dictionary \n""" + \
                   """## Please make your own platform specific template with your own keys and place it in a subfolder of the da package.\n """ + \
                   """## \n""" + \
                   """#SBATCH -J jobname \n""" + \
                   """#SBATCH -p jobqueue \n""" + \
                   """#SBATCH -n jobnodes \n""" + \
                   """#SBATCH -t jobtime \n""" + \
                   """#SBATCH -o joblog \n""" + \
                   """module load python\n""" + \
                   """module load nco\n""" + \
		   """\n"""

        if 'depends' in joboptions:
            template += """#$ -hold_jid depends \n"""

        # First replace from passed dictionary
        for k, v in joboptions.items():
            while k in template:
                template = template.replace(k, v)

        # Fill remaining values with std_options
        for k, v in std_joboptions.items():
            while k in template:
                template = template.replace(k, v)

        return template


    def submit_job(self, jobfile, joblog=None, block=False):
        """ This method submits a jobfile to the queue, and returns the queue ID """


        #cmd     = ["llsubmit","-s",jobfile]
        #msg = "A new task will be started (%s)"%cmd  ; logging.info(msg)

        if block:
            cmd = ["salloc",'-n',std_joboptions['jobnodes'],'-t',std_joboptions['jobtime'], jobfile]
            logging.info("A new task will be started (%s)" % cmd)
            output = subprocess.Popen(cmd, stdout=subprocess.PIPE).communicate()[0]
            logging.info(output)
            print 'output', output
            jobid = output.split()[-1]             
            print 'jobid', jobid
        else:
            cmd = ["sbatch", jobfile]
            logging.info("A new job will be submitted (%s)" % cmd)
            output = subprocess.Popen(cmd, stdout=subprocess.PIPE).communicate()[0]  ; logging.info(output)
            jobid = output.split()[-1]
            
        return jobid







#        jobid   = output.split()[2]
#        retcode = output.split()[-1]
#	
#         #for huygens
#	print 'output', output
#	test   = output.split()[3]
#        dummy, jobid =test.split('nl.')
#	jobid='%s%s' %('"',jobid)
#	submitmsg ='%s%s%s'%(output.split()[4],output.split()[5],output.split()[6])
#        if submitmsg=='hasbeensubmitted.':
#	    retcode=2
#        print 'retcode',submitmsg,retcode	 
#        return retcode
#
#    def KillJob(self,jobid):
#        """ This method kills a running job """
#
#        output = subprocess.Popen(['qdel',jobid], stdout=subprocess.PIPE).communicate()[0]  ; logging.info(output)
#
#        return output
#
#    def StatJob(self,jobid):
#        """ This method gets the status of a running job """
#        import subprocess
#
#        #output = subprocess.Popen(['sgestat'], stdout=subprocess.PIPE).communicate()[0]  ; logging.info(output)
#
#        return ''
   


if __name__ == "__main__":
    pass
