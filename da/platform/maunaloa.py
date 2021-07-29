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

from da.baseclasses.platform import Platform, std_joboptions

class MaunaloaPlatform(Platform):
    def __init__(self):
        self.ID = 'WU maunaloa'    # the identifier gives the platform name
        self.version = '1.0'     # the platform version used

    def give_blocking_flag(self):
        return ""

    def give_queue_type(self):
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

        template = """## \n""" + \
                   """## This is a set of dummy names, to be replaced by values from the dictionary \n""" + \
                   """## Please make your own platform specific template with your own keys and place it in a subfolder of the da package.\n """ + \
                   """## \n""" + \
                   """ \n""" + \
                   """#$ jobname \n""" + \
                   """#$ jobaccount \n""" + \
                   """#$ jobnodes \n""" + \
                   """#$ jobtime \n""" + \
                   """#$ jobshell \n""" + \
                   """\n""" + \
                   """source /usr/local/Modules/3.2.8/init/sh\n""" + \
                   """module load python\n""" + \
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


if __name__ == "__main__":
    pass
