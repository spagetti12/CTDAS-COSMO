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
# maunaloa.py

"""
Author : peters 

Revision History:
File created on 06 Sep 2010.

"""

import logging
import subprocess

from da.baseclasses.platform import Platform

std_joboptions = {'jobname':'test', 'jobaccount':'co2', 'jobtype':'serial', 'jobshell':'/bin/sh', 'depends':'', 'jobtime':'24:00:00', 'jobinput':'/dev/null', 'jobnode':'', 'jobtasks':'', 'modulenetcdf':'netcdf/4.1.2', 'networkMPI':''}


class HuygensPlatform(Platform):
    def __init__(self):
        self.ID = 'huygens'    # the identifier gives the platform name
        self.version = '1.0'     # the platform version used


    def give_blocking_flag(self):
        """
        Returns a blocking flag, which is important if tm5 is submitted in a queue system. The python ctdas code is forced to wait before tm5 run is finished

	    -on Huygens: return "-s"
            -on Maunaloa: return "" (no queue available)
            -on Jet/Zeus: return
        """
        return "-s"

    def give_queue_type(self):
        """
        Return a queue type depending whether your computer system has a queue system, or whether you prefer to run in the foreground. 
        On most large systems using the queue is mandatory if you run a large job.
            -on Huygens: return "queue"
            -on Maunaloa: return "foreground" (no queue available)
            -on Jet/Zeus: return  

        """
        return "queue"

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
		   """## @ node_usage = normal\n""" + \
                   """jobnode \n""" + \
                   """jobtasks \n""" + \
                   """networkMPI \n""" + \
                   """# @ notification = never\n""" + \
                   """# @ input = jobinput\n""" + \
                   """# @ output = logfile.$(jobid)\n""" + \
                   """# @ error = logfile.$(jobid)\n""" + \
                   """# @ wall_clock_limit = jobtime\n""" + \
                   """# @ job_type = jobtype \n""" + \
                   """# @ shell = /bin/bash\n""" + \
                   """# @ queue \n""" + \
                   """\n""" + \
                   """module load ctdas\n""" + \
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



   
      #  msg1 = 'Platform initialized: %s' % self.Identifier      ; logging.info(msg1)
#        #msg2  = '%s version: %s'%(self.Identifier,self.Version) ; logging.info(msg2)


#    def GetMyID(self):
#        try:
#            return os.environ['JOB_ID']
#        except:
#            return os.getpid()
#

#    def submit_job(self,jobfile,joblog=None,block=False):
#        """ This method submits a jobfile to the queue, and returns the queue ID """
#
#
#        cmd     = ["llsubmit","-s",jobfile]
#        msg = "A new task will be started (%s)"%cmd  ; logging.info(msg)
#
#        if block:
#            jobid = subprocess.call(cmd)
#	    print 'jobid',jobid
#        else:
#            jobid = subprocess.Popen(cmd).pid
#         #output  = subprocess.Popen(cmd, stdout=subprocess.PIPE).communicate()[0]  ; logging.info(output)
#        code=0




    def submit_job(self, jobfile, joblog=None, block=False):
        """ This method submits a jobfile to the queue, and returns the queue ID """


        #cmd     = ["llsubmit","-s",jobfile]
        #msg = "A new task will be started (%s)"%cmd  ; logging.info(msg)

        if block:
            cmd = ["llsubmit", "-s", jobfile]
            logging.info("A new task will be started (%s)" % cmd)
            output = subprocess.Popen(cmd, stdout=subprocess.PIPE).communicate()[0]
            logging.info(output)
            print 'output', output
            jobid = output.split()[3]             
            print 'jobid', jobid
        else:
            cmd = ["llsubmit", jobfile]
            logging.info("A new task will be started (%s)" % cmd)
            output = subprocess.Popen(cmd, stdout=subprocess.PIPE).communicate()[0]  ; logging.info(output)
            jobid = output.split()[3]
            
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
