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
# ct_statevector_tools.py

"""
.. module:: statevector
.. moduleauthor:: Wouter Peters 

Revision History:
File created on 28 Jul 2010.

The module statevector implements the data structure and methods needed to work with state vectors (a set of unknown parameters to be optimized by a DA system) of different lengths, types, and configurations. Two baseclasses together form a generic framework:
    * :class:`~da.baseclasses.statevector.StateVector`
    * :class:`~da.baseclasses.statevector.EnsembleMember`

As usual, specific implementations of StateVector objects are done through inheritance form these baseclasses. An example of designing 
your own baseclass StateVector we refer to :ref:`tut_chapter5`.

.. autoclass:: da.baseclasses.statevector.StateVector 

.. autoclass:: da.baseclasses.statevector.EnsembleMember 

"""

import os
import logging
import numpy as np
from datetime import timedelta
import da.cosmo.io4 as io
from netCDF4 import Dataset

identifier = 'Baseclass Statevector '
version = '0.0'

################### Begin Class EnsembleMember ###################

class EnsembleMember(object):
    """ 
        An ensemble member object consists of:
           * a member number
           * parameter values
           * an observation object to hold sampled values for this member

        Ensemble members are initialized by passing only an ensemble member number, all data is added by methods 
        from the :class:`~da.baseclasses.statevector.StateVector`. Ensemble member objects have almost no functionality 
        except to write their data to file using method :meth:`~da.baseclasses.statevector.EnsembleMember.write_to_file`

        .. automethod:: da.baseclasses.statevector.EnsembleMember.__init__ 
        .. automethod:: da.baseclasses.statevector.EnsembleMember.write_to_file 
        .. automethod:: da.baseclasses.statevector.EnsembleMember.AddCustomFields 

    """

    def __init__(self, membernumber):
        """
           :param memberno: integer ensemble number
           :rtype: None

           An EnsembleMember object is initialized with only a number, and holds two attributes as containter for later
           data:
                * param_values, will hold the actual values of the parameters for this data
                * ModelSample, will hold an :class:`~da.baseclasses.obs.Observation` object and the model samples resulting from this members' data

        """
        self.membernumber = membernumber   # the member number
        self.param_values = None           # Parameter values of this member

################### End Class EnsembleMember ###################

################### Begin Class StateVector ###################


class StateVector(object):
    """ 
    The StateVector object first of all contains the data structure of a statevector, defined by 3 attributes that define the 
    dimensions of the problem in parameter space:
        * nlag
        * nparameters
        * nmembers

    The fourth important dimension `nobs` is not related to the StateVector directly but is initialized to 0, and later on 
    modified to be used in other parts of the pipeline:
        * nobs

    These values are set as soon as the :meth:`~da.baseclasses.statevector.StateVector.setup` is called from the :ref:`pipeline`. 
    Additionally, the value of attribute `isOptimized` is set to `False` indicating that the StateVector holds a-priori values 
    and has not been modified by the :ref:`optimizer`.

    StateVector objects can be filled with data in two ways
        1. By reading the data from file
        2. By creating the data through a set of method calls

    Option (1) is invoked using method :meth:`~da.baseclasses.statevector.StateVector.read_from_file`. 
    Option (2) consists of a call to method :meth:`~da.baseclasses.statevector.StateVector.make_new_ensemble`

    Once the StateVector object has been filled with data, it is used in the pipeline and a few more methods are
    invoked from there:
        * :meth:`~da.baseclasses.statevector.StateVector.propagate`, to advance the StateVector from t=t to t=t+1
        * :meth:`~da.baseclasses.statevector.StateVector.write_to_file`, to write the StateVector to a NetCDF file for later use

    The methods are described below:

    .. automethod:: da.baseclasses.statevector.StateVector.setup 
    .. automethod:: da.baseclasses.statevector.StateVector.read_from_file
    .. automethod:: da.baseclasses.statevector.StateVector.write_to_file
    .. automethod:: da.baseclasses.statevector.StateVector.make_new_ensemble
    .. automethod:: da.baseclasses.statevector.StateVector.propagate
    .. automethod:: da.baseclasses.statevector.StateVector.write_members_to_file

    Finally, the StateVector can be mapped to a gridded array, or to a vector of TransCom regions, using:

    .. automethod:: da.baseclasses.statevector.StateVector.grid2vector
    .. automethod:: da.baseclasses.statevector.StateVector.vector2grid
    .. automethod:: da.baseclasses.statevector.StateVector.vector2tc
    .. automethod:: da.baseclasses.statevector.StateVector.state2tc

    """

    def __init__(self):
        self.ID = identifier
        self.version = version

        # The following code allows the object to be initialized with a dacycle object already present. Otherwise, it can
        # be added at a later moment.

        logging.info('Statevector object initialized: %s' % self.ID)

    def setup(self, dacycle):
        """
        setup the object by specifying the dimensions. 
        There are two major requirements for each statvector that you want to build:
        
            (1) is that the statevector can map itself onto a regular grid
            (2) is that the statevector can map itself (mean+covariance) onto TransCom regions

        An example is given below.
        """

        self.nlag = int(dacycle['time.nlag'])
        self.nmembers = int(dacycle['da.optimizer.nmembers'])
#        self.nparams = int(dacycle.dasystem['nparameters'])
    #    self.nparams = 21
        self.nparams = 23
        self.nobs = 0
        
        self.obs_to_assimilate = ()  # empty containter to hold observations to assimilate later on

        # These list objects hold the data for each time step of lag in the system. Note that the ensembles for each time step consist 
        # of lists of EnsembleMember objects, we define member 0 as the mean of the distribution and n=1,...,nmembers as the spread.

        self.ensemble_members = list(range(self.nlag))

        for n in range(self.nlag):
            self.ensemble_members[n] = []


        # This specifies the file to read with the gridded mask at 1x1 degrees. Each gridbox holds a number that specifies the parametermember
        #  that maps onto it. From this map, a dictionary is created that allows a reverse look-up so that we can map parameters to a grid.

        mapfile = os.path.join(dacycle.dasystem['regionsfile'])
#        ncf = io.ct_read(mapfile, 'read')
#        self.gridmap = ncf.get_variable('regions')
#        self.tcmap = ncf.get_variable('transcom_regions')
#        ncf.close()
        regfile = Dataset(mapfile,mode='r')
        self.gridmap = np.squeeze(regfile.variables['regions'])
        self.tcmap = np.squeeze(regfile.variables['transcom_regions'])
        regfile.close()

        logging.debug("A TransCom  map on 1x1 degree was read from file %s" % dacycle.dasystem['regionsfile'])
        logging.debug("A parameter map on 1x1 degree was read from file %s" % dacycle.dasystem['regionsfile'])

        # Create a dictionary for state <-> gridded map conversions

        nparams = 23
        #nparams = self.gridmap.max()
        self.griddict = {}
#        for r in range(1, 10):
        for r in range(1, 11):
#        for r in range(1, int(nparams) + 1):
            sel = np.nonzero(self.gridmap.flat == r)
            if len(sel[0]) > 0: 
                self.griddict[r] = sel
                #self.griddict[r+10] = sel   # pavle - expand dictionary for nparam values because of RESP
                self.griddict[r+11] = sel   # pavle - expand dictionary for nparam values because of RESP
        # pavle: sel sorts out regions by PFT

        logging.debug("A dictionary to map grids to states and vice versa was created")

        # Create a matrix for state <-> TransCom conversions

#        self.tcmatrix = np.zeros((self.nparams, 21), 'float') 
        self.tcmatrix = np.zeros((self.nparams, 23), 'float') 

        for r in range(1, self.nparams + 1):
            sel = np.nonzero(self.gridmap.flat == r)
            if len(sel[0]) < 1: 
                continue
            else:
                n_tc = set(self.tcmap.flatten().take(sel[0]))
                if len(n_tc) > 1: 
                    logging.error("Parameter %d seems to map to multiple TransCom regions (%s), I do not know how to handle this" % (r, n_tc))
                    raise ValueError
                self.tcmatrix[r - 1, int(n_tc.pop()) - 1] = 1.0

        logging.debug("A matrix to map states to TransCom regions and vice versa was created")

        # Create a mask for species/unknowns

        self.make_species_mask()

    def make_species_mask(self):

        """

        This method creates a dictionary with as key the name of a tracer, and as values an array of 0.0/1.0 values 
        specifying which StateVector elements are constrained by this tracer. This mask can be used in 
        the optimization to ensure that certain types of osbervations only update certain unknowns.

        An example would be that the tracer '14CO2' can be allowed to only map onto fossil fuel emissions in the state

        The form of the mask is:

        {'co2': np.ones(self.nparams), 'co2c14', np.zeros(self.nparams)  }

        so that 'co2' maps onto all parameters, and 'co2c14' on none at all. These arrays are used in the Class 
        optimizer when state updates are actually performed

        """
        self.speciesdict = {'co2': np.ones(self.nparams)}
        logging.debug("A species mask was created, only the following species are recognized in this system:")
        for k in self.speciesdict.keys(): 
            logging.debug("   ->    %s" % k)


    def make_new_ensemble(self, lag, covariancematrix=None):
        """ 
        :param lag: an integer indicating the time step in the lag order
        :param covariancematrix: a matrix to draw random values from
        :rtype: None
    
        Make a new ensemble, the attribute lag refers to the position in the state vector. 
        Note that lag=1 means an index of 0 in python, hence the notation lag-1 in the indexing below.
        The argument is thus referring to the lagged state vector as [1,2,3,4,5,..., nlag]

        The optional covariance object to be passed holds a matrix of dimensions [nparams, nparams] which is
        used to draw ensemblemembers from. If this argument is not passed it will ne substituted with an 
        identity matrix of the same dimensions.

        """    

#        if covariancematrix == None: 
 #           covariancematrix = np.identity(self.nparams)

        # Make a cholesky decomposition of the covariance matrix


        try:
            _, s, _ = np.linalg.svd(covariancematrix)
        except:
            s = np.linalg.svd(covariancematrix, full_matrices=1, compute_uv=0) #Cartesius fix
        dof = np.sum(s) ** 2 / sum(s ** 2)
#        covariancematrix = np.identity(self.nparams)   # pavle - for testing, final matrix will be 23x23
        C = np.linalg.cholesky(covariancematrix)

        logging.debug('Cholesky decomposition has succeeded ')
        logging.info('Appr. degrees of freedom in covariance matrix is %s' % (int(dof)))


        # Create mean values 

        newmean = np.ones(self.nparams, float) # standard value for a new time step is 1.0

        # If this is not the start of the filter, average previous two optimized steps into the mix
        if lag == self.nlag - 1 and self.nlag >= 3:
            newmean += self.ensemble_members[lag - 1][0].param_values + \
                                           self.ensemble_members[lag - 2][0].param_values 
            newmean = newmean / 3.0

        # Create the first ensemble member with a deviation of 0.0 and add to list

        newmember = EnsembleMember(0)
        newmember.param_values = newmean.flatten()  # no deviations
        self.ensemble_members[lag].append(newmember)

        # Create members 1:nmembers and add to ensemble_members list
        for member in range(1, self.nmembers):
            #rands = np.random.uniform(low=-0.5, high=0.5, size=self.nparams-1)
            rands = np.random.uniform(low=-1., high=1., size=self.nparams-1)
            rands_bg = np.random.uniform(low=-0.05, high=0.05, size=1)

#            if member == 1 and lag == 0:
 #               rands=np.array([-0.5402309, 0.3021888, 0.346095, 0.4942706, 0.6276228, 0.7450897, 0.5155229, -0.3006688, 0.8386819, -0.1197997, -1, -0.9025682, -1,    0.5018581, 0.923963, -0.2253573, 0.4716174, 0.4637589, 0.7411034,     0.6526616, -0.3188075, -1, 0.01435228])
  #          if member == 2 and lag == 0:
   #             rands=np.array([0.1876527, -0.4323478, 0.2562822, -0.7032343, 0.83417, 0.1378682, 0.3969375, 0.2101818, 0.5018877, 0.1390382, -1, 0.9267841, 0.1679351,     0.1019922, 0.4781876, -0.08182608, 0.9295725, 0.05694205, -0.119032,     -0.8025782, -0.8610154, -1, -0.03559723])
    #        if member == 3 and lag == 0:
     #           rands=np.array([-0.8030154, 0.08172598, -0.7622297, -0.4516532, 0.4969012, 0.8023733,     -0.08237687, -0.8089882, -0.8048904, -0.4082271, -1, -0.1105984,     -0.4859152, 0.3398715, -0.09020101, 0.462075, 0.5386102, -0.865896,     0.5063729, -0.01304542, -0.2682804, -1, 0.04435889])
      #      if member == 4 and lag == 0:
       #         rands=np.array([-0.555479, -0.6089986, -0.664399, -0.03338498, -1, -1, 0.1820136,     0.1071932, -0.9634671, -1, -1, -0.8547595, 0.273797, -0.802806,     0.1928162, -0.5642408, 0.6471559, 0.2174408, -0.05968733, -0.004075899,     0.9075829, -1, -0.02335107])
        #    if member == 5 and lag == 0:
         #       rands=np.array([-0.674361, 0.3454686, -0.4433188, -0.8486685, -0.7747395, 0.203481, -1,     -0.09257206, -0.3892125, -0.02775251, -1, -0.6137948, -0.8518382,     -0.8779979, -0.6405048, -0.8992591, -0.6312428, 0.4434872, -0.05931631,     -0.2671082, -0.7183358, -1, -0.02273439])
          #  if member == 6 and lag == 0:
           #     rands=np.array([-0.4753621, -1, 0.6610132, -0.3525974, -0.987463, -0.2421917, -0.8256455,     0.1535834, 0.6603956, 0.2776085, -1, -0.6808311, -0.6428144, -1,     -0.518185, -1, -1, 0.4222396, 0.7343495, 0.6562976, -0.2627398, -1,     0.00420737])
            #if member == 7 and lag == 0:
             #   rands=np.array([0.9235092, -0.08976277, -0.693751, 0.7479168, -0.5282985, -0.325265,     0.5994204, 0.1087474, -0.8495662, 0.5450585, -1, -0.4888758, -0.6259907,     0.5597795, 0.7957715, 0.04843958, 0.2360165, -0.7671149, 0.2514055,     -0.872794, -0.5365307, -1, -0.02834922])
#            if member == 8 and lag == 0:
 #               rands=np.array([-0.2851942, 0.2424788, 0.6036716, 0.1861704, 0.2724828, -0.1842902,     0.4256651, 0.9268349, -0.5386767, -0.282816, -1, 0.1387041, 0.6110855,     -0.5958148, -0.9803032, -0.5043366, -0.6089351, 0.2659722, -0.7406561,     0.409389, -0.07663137, -1, -0.04040922])
  #          if member == 9 and lag == 0:
   #             rands=np.array([0.7698311, 0.2813426, 0.4236566, -0.5569826, 0.06671345, 0.1696645,     -0.5119113, 0.5269542, 0.4109541, 0.5580103, -1, -0.4272517, -0.556591,     0.1682963, 0.02460092, -0.2566992, -0.9167355, -0.7835014, 0.7101904,     -0.5695945, -0.1746321, -1, 0.04223957])
    #        if member == 10 and lag == 0:
     #           rands=np.array([-0.6941028, -0.9384637, -0.9250847, -1, -0.959329, -1, 0.04787683,     -0.896835, -0.4911715, -0.01541466, -1, 0.9626108, 0.6731427, 0.5944214,     1.05593, 1.0805, -0.1627321, -0.154229, 0.6146008, 0.420991, 0.7912573,     -1, -0.02473733])
      #      if member == 11 and lag == 0:
       #         rands=np.array([0.96286, 0.9415601, -0.7406601, 0.119507, 0.5506695, 0.5813184, 0.9465097,     -0.658381, -0.5169885, 0.3601404, -1, -0.1221153, -0.911568, 0.4350488,     0.983259, 0.4603367, 0.4221282, -0.3391019, -0.192896, 0.3793332,     0.8009401, -1, 0.03354285])
        #    if member == 12 and lag == 0:
         #       rands=np.array([0.4491677, -0.5364176, 0.1777764, 0.3457868, 0.4886276, 0.6140389,     -0.9525508, 0.7758317, 0.7040095, 0.8019212, -1, 0.2732197, 0.994552,     -0.199925, 0.7970242, 0.5645052, -0.6018554, 0.01864898, 0.0417577,     0.3765279, -0.2793771, -1, -0.008472145])
          #  if member == 13 and lag == 0:
           #     rands=np.array([-0.8257746, 0.5533802, -0.7910783, 0.4864411, -0.5154006, -0.5978034,     -0.202793, -0.9381504, 0.2701632, 0.3198348, -1, 0.3341525, 1.001176,     -0.2380723, 0.1933595, -0.6892359, -0.4271372, 0.6341484, -0.5026917,     -0.3053397, 0.750175, -1, -0.02964199])
            #if member == 14 and lag == 0:
             #   rands=np.array([0.701148, -0.1619855, -0.063532, -0.4197367, -0.6209376, -1, -0.8872383,     -0.4309417, -1, -1, -1, 0.1452325, -0.2525761, 0.1687515, -0.3141938,     -0.1218985, -0.4126042, 0.1641715, 0.3630136, -0.5947977, 0.7974821, -1,     -0.004756214])
#            if member == 15 and lag == 0:
 #               rands=np.array([-0.2590746, -0.5643597, -0.2787321, -0.2195021, -0.7023828, 0.5482118,     0.6232124, 0.5892509, 0.6535864, -0.5035124, -1, -0.6891443, 0.3709087,     -0.1523427, -0.9473059, 0.1686965, 0.04820674, 0.07875252, 0.5189905,     0.2113448, -0.708707, -1, -0.01257968])
  #          if member == 16 and lag == 0:
   #             rands=np.array([0.4056181, 0.04536337, -0.6866099, -0.6401204, -0.2879847, -0.7785423,     0.6233089, 0.2740796, -0.8484439, 0.9647028, -1, 0.789863, 1.163193,     -0.2164443, 0.08588136, 0.5577631, -0.3278782, -0.1751566, -0.7633625,     0.7689492, -0.819483, -1, -0.03903025])
    #        if member == 17 and lag == 0:
     #           rands=np.array([0.2155399, 0.1660601, 0.4649725, 0.5862293, -0.4452171, -0.7389292,     -0.07185336, 0.9012561, -0.2824533, -0.7226023, -1, 0.4194106,     -0.2461459, 0.109305, 0.2449136, -0.09785557, 0.7447503, -0.1461458,     0.3223822, 0.4597907, 0.9065241, -1, -0.0007600483])
      #      if member == 18 and lag == 0:
       #         rands=np.array([0.04116038, 0.7824114, 0.4636176, 0.1557975, 0.4527292, -0.4133459,     0.994531, -0.3580682, 0.5178631, 0.2745949, -1, 0.5335826, 0.265704,     0.1912572, 0.536295, -0.7093778, 0.07351156, -0.2118817, -0.7678121,     0.43079, 0.512942, -1, 0.03622137])
        #    if member == 19 and lag == 0:
         #       rands=np.array([0.3184969, 0.2478539, 0.1378153, -0.6565253, 0.6140021, -0.2936001, 0.6872,     0.8872706, 0.4869822, -0.2930091, -1, 0.2853186, 0.9346947, 0.7716588,     -0.1080804, 0.3278734, -0.2541114, -0.7900179, -0.259804, 0.6708161,     -0.8029982, -1, 0.000798168])
          #  if member == 20 and lag == 0:
           #     rands=np.array([-0.7732743, -0.3873762, -0.8547782, -0.07859742, -1, -0.4397394, 0.2048925,     -0.3029643, 0.7081758, -0.06017716, -1, -0.466458, -0.2800255, 0.1250267,     0.632107, -0.7519485, 0.4637258, 0.5653031, 0.4806009, -0.01624972,     0.05281375, -1, 0.01537185])
            #if member == 1 and lag == 1:
             #   rands=np.array([0.4061777, -0.5013675, -0.2862317, 0.6116837, 0.05803703, -0.3041813,     0.4205766, -0.04346657, -0.4199566, -0.2088348, -1, 0.8685411,     -0.1089062, 0.8366363, -0.278753, -0.4832861, -0.147602, 0.3227859,     0.2312384, 0.357722, -0.03119155, -1, -0.03759673])
#            if member == 2 and lag == 1:
 #               rands=np.array([-0.3528758, 0.2302615, 0.8945007, -0.2294943, -0.7223361, 0.4254291,     0.2018437, 0.555488, -0.03991819, 0.2198012, -1, 0.918696, 0.3752462,     1.003467, -0.3229041, 0.6925161, 1.036227, 0.7551609, 1.008353,     0.6003665, -0.6694191, -1, 0.03411322])
  #          if member == 3 and lag == 1:
   #             rands=np.array([0.732138, 0.0511689, 0.5308715, -0.3620436, 0.1138427, 0.3123921,     0.7225392, 0.8311706, -0.7413127, 0.02754115, -1, 0.2934451, -0.7651439,     -0.574313, -0.9932647, -0.5833654, -0.5632083, 0.2126374, -0.4117518,     -0.6983544, -0.5544686, -1, 0.03195827])
    #        if member == 4 and lag == 1:
     #           rands=np.array([-0.5030815, -0.02779967, -0.9290833, 0.1386678, -0.510511, 0.297296,     -0.6832802, -0.1534745, -0.7899511, 0.4862387, -1, -0.991854, -0.2682375,     -0.7712117, -0.8905275, 0.3489431, 0.6761383, 0.1792146, 0.8675478,     -0.7689523, -0.8875247, -1, 0.01273509])
      #      if member == 5 and lag == 1:
       #         rands=np.array([0.3466186, 0.6708195, -0.4234087, 0.791083, 0.5170874, -0.4746583,    0.9374641, 0.9230136, -0.5804763, -0.02489379, -1, 0.9292546, 0.3582832,     -0.2736852, -0.4632118, 1.018919, 0.09223775, 0.3311332, 0.2550537,     0.7294372, 0.4026192, -1, 0.03795173])
        #    if member == 6 and lag == 1:
         #       rands=np.array([0.3490212, 0.2809639, 0.08980897, -0.6552433, 0.3108646, 0.9600278,     -0.7445735, -0.929427, 0.6871993, -0.1122855, -1, 0.2775505, 0.7744503,     0.5353921, 0.8993617, 0.5068132, 0.6897388, -0.5868602, -0.5443915,     0.2293534, -0.880945, -1, -0.04769203])
          #  if member == 7 and lag == 1:
           #     rands=np.array([0.7321754, 0.06772549, 0.2179533, 0.2919183, 1.062581, 0.7554626,     0.5785612, 1.011077, 0.6109563, -0.350915, -1, -0.579572, -0.07610841,     -0.3409401, 0.7316601, -0.904253, -1, -0.01014853, 0.5529011, 0.534495,     -0.3398449, -1, -0.007357216])
            #if member == 8 and lag == 1:
             #   rands=np.array([-0.4717625, -0.881941, -0.1655097, 0.2257761, 0.8264135, -0.7091568,     0.2949457, -0.296845, 0.794414, -0.743104, -1, 0.9034545, -0.2101449,     0.3021894, -0.4184734, -0.4638716, 0.6389658, 0.9005892, -0.5002754,     0.6826805, -0.3747995, -1, 0.03224342])
#            if member == 9 and lag == 1:
 #               rands=np.array([-0.5808461, -0.8636444, 0.7052112, 0.7909082, -0.8140187, -0.02347049,     0.5638158, -0.7550521, 0.5693654, 0.07545193, -1, -0.2417251, 0.6826159,     0.9508564, -0.3956808, 0.8499226, 0.09175592, 0.5388746, 0.8177778,     0.4029195, 0.4111119, -1, -0.01329621])
  #          if member == 10 and lag == 1:
   #             rands=np.array([0.2166863, -0.437274, -0.7773396, 0.3759225, -0.7887743, -0.5783228,     0.2769741, -0.8951628, 0.5877202, -0.7473742, -1, 0.8568208, -0.3967594,     -0.02381929, 0.911634, 0.1854728, 0.347165, -0.7240146, 0.6684847,     -0.05459603, -0.2197302, -1, -0.04627957])
    #        if member == 11 and lag == 1:
     #           rands=np.array([-0.5553841, 0.1981317, -0.01706193, 0.3010842, 0.3364336, 0.6311802,     -0.4822653, -1, -0.2353798, -0.4258329, -1, -0.3820901, -0.1350132,     0.1199797, -0.4477318, -0.03276956, 0.2006144, 0.7612231, -0.4382411,     0.3140147, -0.7996647, -1, -0.0351413])
      #      if member == 12 and lag == 1:
       #         rands=np.array([0.1795571, 0.7117739, 0.3948795, 0.01603693, 0.4342803, 0.3915263,     -0.1187998, 0.3347789, -0.1160258, -0.6199519, -1, 0.9331298, 0.2128856,     0.1029978, 0.5622161, -0.03326744, 0.1113134, 0.718336, -0.7506785,     0.8058456, -0.165978, -1, 0.03335277])
        #    if member == 13 and lag == 1:
         #       rands=np.array([0.1133219, 0.1868301, 0.6434304, 0.06698409, 0.7134196, 0.6549824,     -0.3687601, 0.2397162, 0.3933871, -0.8959517, -1, -0.9083285, 0.1290343,     -0.8501528, -0.1640191, -0.04551715, 0.2021651, 0.1527653, 0.4235286,     0.7385602, -0.1721378, -1, 0.008510991])
          #  if member == 14 and lag == 1:
           #     rands=np.array([0.5128248, -0.025835, -0.3690266, 0.5315114, -0.3562845, -0.8088912,     0.6070702, 0.8863306, 0.127363, 0.4851948, -1, 0.2463859, 0.5403649,     0.1009353, -0.3769973, 0.8082039, 0.642301, -0.3246575, -0.8921521,     0.485802, 0.06472221, -1, 0.042198])
            #if member == 15 and lag == 1:
             #   rands=np.array([-0.6925536, -0.1888134, 0.7955419, -0.1718084, -0.4925962, -0.63966,     0.712366, 0.7542216, 0.6485603, 0.9584286, -1, 0.2036391, -0.547542,     -0.1605409, 0.5345823, -0.8189099, -0.4620013, 0.2116281, 0.5273049,     -0.5836464, -0.7646005, -1, 0.03012722])
#            if member == 16 and lag == 1:
 #               rands=np.array([-0.9031112, -0.9747161, 0.6262707, 0.1299344, -0.656873, -1, 0.7236382,     0.9888599, 0.599555, -0.3344736, -1, 0.1587437, 0.9665449, 0.7327783,     0.2025174, 1.045059, 0.7680195, 0.202516, 0.3561991, -0.08263882,     -0.7926366, -1, -0.0340635])
  #          if member == 17 and lag == 1:
   #             rands=np.array([-0.4468101, -0.6611752, 0.02045458, -0.818266, -0.06696166, -0.2870016,     -0.8372102, -1, 0.6159038, 0.1100874, -1, -0.2941307, -0.8388062,     0.7502818, -0.02831579, -0.4171109, 0.3452238, 0.3286322, -0.1746102,     0.741174, -0.1962416, -1, -0.03235182])
    #        if member == 18 and lag == 1:
     #           rands=np.array([-0.117545, 0.6358077, -0.8399339, -0.4401485, -0.4653139, -0.1412946,     0.368727, 0.5386308, 0.07112807, 0.2995664, -1, -0.3923609, -0.8697481,     0.283159, -0.1124206, -0.2029478, -1, -0.9128151, 0.2968773, -0.4699179,     -0.1645328, -1, -0.04468708])
      #      if member == 19 and lag == 1:
       #         rands=np.array([-0.3586554, -0.8109571, -0.4370517, 0.1238612, -0.996357, -0.6696278,     0.5044641, 0.6546122, -0.0149131, 0.8158504, -1, 0.8035823, 0.565383,     0.2036764, -0.6208936, 0.4221191, 0.4371626, -0.6936532, -0.6731901, -1,     0.1557615, -1, -0.02435511])
        #    if member == 20 and lag == 1:
         #       rands=np.array([-0.6186742, 0.1736774, 0.3401965, 0.8706512, -0.6920347, -0.2222075,     0.1594849, 0.3773506, 1.039722, 0.3430154, -1, 0.8352696, 1.13694,     -0.20392, 0.1833831, 1.04152, 0.1299644, 0.07491443, -0.6444398,     0.2261413, 0.01262438, -1, -0.03116111])

            newmember = EnsembleMember(member)
         #   newmember.param_values = rands + newmean
       #     newmember.param_values = (np.hstack((np.dot(C, rands[0:10]),np.dot(C, rands[10:20]), rands_bg)) + newmean).ravel()
            newmember.param_values = (np.hstack((np.dot(C, rands[0:11]),np.dot(C, rands[11:22]), rands_bg)) + newmean).ravel()
         #   newmember.param_values[10] = 0.
          #  newmember.param_values[21] = 0.
            newmember.param_values[newmember.param_values<0.] = 0.
            newmember.param_values[newmember.param_values>2.] = 2.
            self.ensemble_members[lag].append(newmember)

        logging.debug('%d new ensemble members were added to the state vector # %d' % (self.nmembers, (lag + 1)))


    def propagate(self, dacycle):
        """ 
        :rtype: None

        Propagate the parameter values in the StateVector to the next cycle. This means a shift by one cycle 
        step for all states that will
        be optimized once more, and the creation of a new ensemble for the time step that just 
        comes in for the first time (step=nlag). 
        In the future, this routine can incorporate a formal propagation of the statevector.

        """
        
        # Remove State Vector n=1 by simply "popping" it from the list and appending a new empty list at the front. This empty list will
        # hold the new ensemble for the new cycle 

        self.ensemble_members.pop(0)
        self.ensemble_members.append([])

        # And now create a new time step of mean + members for n=nlag
        date = dacycle['time.start'] + timedelta(days=(self.nlag - 0.5) * int(dacycle['time.cycle']))
        cov = self.get_covariance(date, dacycle)
        self.make_new_ensemble(self.nlag - 1, cov)

        logging.info('The state vector has been propagated by one cycle')


    def write_to_file(self, filename, qual):
        """
        :param filename: the full filename for the output NetCDF file
        :rtype: None

        Write the StateVector information to a NetCDF file for later use. 
        In principle the output file will have only one two datasets inside 
        called:
            * `meanstate`, dimensions [nlag, nparamaters]
            * `ensemblestate`, dimensions [nlag,nmembers, nparameters]

        This NetCDF information can be read back into a StateVector object using 
        :meth:`~da.baseclasses.statevector.StateVector.read_from_file`

        """
        #import da.tools.io4 as io
        #import da.tools.io as io

        if qual == 'prior':
            f = io.CT_CDF(filename, method='create')
            logging.debug('Creating new StateVector output file (%s)' % filename)
            #qual = 'prior'
        else:
            f = io.CT_CDF(filename, method='write')
            logging.debug('Opening existing StateVector output file (%s)' % filename)
            #qual = 'opt'

        dimparams = f.add_params_dim(self.nparams)
        dimmembers = f.add_members_dim(self.nmembers)
        dimlag = f.add_lag_dim(self.nlag, unlimited=True)

        for n in range(self.nlag):
            members = self.ensemble_members[n]
            mean_state = members[0].param_values

            savedict = f.standard_var(varname='meanstate_%s' % qual)
            savedict['dims'] = dimlag + dimparams 
            savedict['values'] = mean_state
            savedict['count'] = n
            savedict['comment'] = 'this represents the mean of the ensemble'
            f.add_data(savedict)

            members = self.ensemble_members[n]
            devs = np.asarray([m.param_values.flatten() for m in members])
            data = devs - np.asarray(mean_state)

            savedict = f.standard_var(varname='ensemblestate_%s' % qual)
            savedict['dims'] = dimlag + dimmembers + dimparams 
            savedict['values'] = data
            savedict['count'] = n
            savedict['comment'] = 'this represents deviations from the mean of the ensemble'
            f.add_data(savedict)
        f.close()

        logging.info('Successfully wrote the State Vector to file (%s) ' % filename)

    def read_from_file(self, filename, qual='opt'):
        """ 
        :param filename: the full filename for the input NetCDF file
        :param qual: a string indicating whether to read the 'prior' or 'opt'(imized) StateVector from file
        :rtype: None

        Read the StateVector information from a NetCDF file and put in a StateVector object
        In principle the input file will have only one four datasets inside 
        called:
            * `meanstate_prior`, dimensions [nlag, nparamaters]
            * `ensemblestate_prior`, dimensions [nlag,nmembers, nparameters]
            * `meanstate_opt`, dimensions [nlag, nparamaters]
            * `ensemblestate_opt`, dimensions [nlag,nmembers, nparameters]

        This NetCDF information can be written to file using 
        :meth:`~da.baseclasses.statevector.StateVector.write_to_file`

        """

        #import da.tools.io as io
        f = io.ct_read(filename, 'read')
        meanstate = f.get_variable('statevectormean_' + qual)
        ensmembers = f.get_variable('statevectorensemble_' + qual)
        f.close()

        for n in range(self.nlag):
            if not self.ensemble_members[n] == []:
                self.ensemble_members[n] = []
                logging.warning('Existing ensemble for lag=%d was removed to make place for newly read data' % (n + 1))

            for m in range(self.nmembers):
                newmember = EnsembleMember(m)
                newmember.param_values = ensmembers[n, m, :].flatten() + meanstate[n]  # add the mean to the deviations to hold the full parameter values
                self.ensemble_members[n].append(newmember)

        logging.info('Successfully read the State Vector from file (%s) ' % filename)

    def write_members_to_file(self, lag, outdir,endswith='.nc'):
        """ 
           :param: lag: Which lag step of the filter to write, must lie in range [1,...,nlag]
           :param: outdir: Directory where to write files
           :param: endswith: Optional label to add to the filename, default is simply .nc
           :rtype: None

           Write ensemble member information to a NetCDF file for later use. The standard output filename is 
           *parameters.DDD.nc* where *DDD* is the number of the ensemble member. Standard output file location 
           is the `dir.input` of the dacycle object. In principle the output file will have only two datasets inside 
           called `parametervalues` which is of dimensions `nparameters` and `parametermap` which is of dimensions (180,360). 
           This dataset can be read and used by a :class:`~da.baseclasses.observationoperator.ObservationOperator` object. 

           .. note:: if more, or other information is needed to complete the sampling of the ObservationOperator you
                     can simply inherit from the StateVector baseclass and overwrite this write_members_to_file function.

        """

        # These import statements caused a crash in netCDF4 on MacOSX. No problems on Jet though. Solution was
        # to do the import already at the start of the module, not just in this method.
           
        #import da.tools.io as io
        #import da.tools.io4 as io

        members = self.ensemble_members[lag]

        for mem in members:
#            filename = os.path.join(outdir, 'parameters.%03d%s' % (mem.membernumber, endswith))
            filename = os.path.join(outdir, 'parameters_lag'+str(lag)+'.%03d%s' % (mem.membernumber, endswith))
            ncf = io.CT_CDF(filename, method='create')
            dimparams = ncf.add_params_dim(self.nparams)
            dimgrid = ncf.add_latlon_dim()

            data = mem.param_values

            savedict = io.std_savedict.copy()
            savedict['name'] = "parametervalues"
            savedict['long_name'] = "parameter_values_for_member_%d" % mem.membernumber
            savedict['units'] = "unitless"
            savedict['dims'] = dimparams 
            savedict['values'] = data
            savedict['comment'] = 'These are parameter values to use for member %d' % mem.membernumber
            ncf.add_data(savedict)

            griddata = self.vector2grid(vectordata=data)

            savedict = io.std_savedict.copy()
            savedict['name'] = "parametermap"
            savedict['long_name'] = "parametermap_for_member_%d" % mem.membernumber
            savedict['units'] = "unitless"
            savedict['dims'] = dimgrid 
            savedict['values'] = griddata.tolist()
            savedict['comment'] = 'These are gridded parameter values to use for member %d' % mem.membernumber
            ncf.add_data(savedict)

            ncf.close()

            logging.debug('Successfully wrote data from ensemble member %d to file (%s) ' % (mem.membernumber, filename))

    def write_members_for_cosmo(self, lag, outdir,endswith='.nc'):
        members = self.ensemble_members[lag]

        self.nparams=int(self.nparams)
        for mem in members:

# GPP
            filename_gpp = os.path.join(outdir, 'parameters_gpp_lag'+str(lag)+'.%03d%s' % (mem.membernumber, endswith))
            ncf = io.CT_CDF(filename_gpp, method='create')
            dimparams = ncf.add_params_dim(11)
#            dimparams = ncf.add_params_dim(10)
            #dimparams = ncf.add_params_dim((self.nparams-1)//2)
            dimgrid = ncf.add_latlon_dim()

        #    data = np.append(mem.param_values[0:10],1.)
            data = mem.param_values[0:11]
            #data = 2.*np.ones(11)-mem.param_values[0:11]       # unblock this maybe
#            print(str(lag),data)
            #data = mem.param_values[0:(self.nparams-1)//2]

            savedict = io.std_savedict.copy()
            savedict['name'] = "parametervalues"
            savedict['long_name'] = "parameter_values_for_member_%d" % mem.membernumber
            savedict['units'] = "unitless"
            savedict['dims'] = dimparams 
            savedict['values'] = data
            savedict['comment'] = 'These are parameter values to use for member %d' % mem.membernumber
            ncf.add_data(savedict)

            griddata = self.vector2grid(vectordata=data)

            savedict = io.std_savedict.copy()
            savedict['name'] = "parametermap"
            savedict['long_name'] = "parametermap_for_member_%d" % mem.membernumber
            savedict['units'] = "unitless"
            savedict['dims'] = dimgrid 
            savedict['values'] = griddata.tolist()
            savedict['comment'] = 'These are gridded parameter values to use for member %d' % mem.membernumber
            ncf.add_data(savedict)

            ncf.close()
# RESP
            filename_resp = os.path.join(outdir, 'parameters_resp_lag'+str(lag)+'.%03d%s' % (mem.membernumber, endswith))
            ncf = io.CT_CDF(filename_resp, method='create')
            #dimparams = ncf.add_params_dim(10)
            dimparams = ncf.add_params_dim(11)
            #dimparams = ncf.add_params_dim((self.nparams-1)//2)
            dimgrid = ncf.add_latlon_dim()

            #data = mem.param_values[10:20]
            data = mem.param_values[11:22]
            #data = np.append(mem.param_values[10:20],1.)
            #data = mem.param_values[(self.nparams-1)//2:self.nparams-1]

            savedict = io.std_savedict.copy()
            savedict['name'] = "parametervalues"
            savedict['long_name'] = "parameter_values_for_member_%d" % mem.membernumber
            savedict['units'] = "unitless"
            savedict['dims'] = dimparams 
            savedict['values'] = data
            savedict['comment'] = 'These are parameter values to use for member %d' % mem.membernumber
            ncf.add_data(savedict)

            griddata = self.vector2grid(vectordata=data)

            savedict = io.std_savedict.copy()
            savedict['name'] = "parametermap"
            savedict['long_name'] = "parametermap_for_member_%d" % mem.membernumber
            savedict['units'] = "unitless"
            savedict['dims'] = dimgrid 
            savedict['values'] = griddata.tolist()
            savedict['comment'] = 'These are gridded parameter values to use for member %d' % mem.membernumber
            ncf.add_data(savedict)

            ncf.close()

# BACKGROUND CO2
#            filename_bg = os.path.join(outdir, 'parameters_bg_lag'+str(lag)+'.%03d%s' % (mem.membernumber, endswith))
 #           ncf = io.CT_CDF(filename_bg, method='create')
  #          dimparams = ncf.add_params_dim(11)
   #         dimgrid = ncf.add_latlon_dim()
#
 #           data = mem.param_values[-1]
  #          data = np.hstack((data, data, data, data, data, data, data, data, data, data, data)).ravel()
#
 #           savedict = io.std_savedict.copy()
  #          savedict['name'] = "parametervalues"
   #         savedict['long_name'] = "parameter_values_for_member_%d" % mem.membernumber
    #        savedict['units'] = "unitless"
     #       savedict['dims'] = dimparams 
      #      savedict['values'] = data
       #     savedict['comment'] = 'These are parameter values to use for member %d' % mem.membernumber
        #    ncf.add_data(savedict)
#
 #           griddata = self.vector2grid(vectordata=data)
#
 #           savedict = io.std_savedict.copy()
  #          savedict['name'] = "parametermap"
   #         savedict['long_name'] = "parametermap_for_member_%d" % mem.membernumber
    #        savedict['units'] = "unitless"
     #       savedict['dims'] = dimgrid 
      #      savedict['values'] = griddata.tolist()
       #     savedict['comment'] = 'These are gridded parameter values to use for member %d' % mem.membernumber
        #    ncf.add_data(savedict)
#
 #           ncf.close()
#
#
    def grid2vector(self, griddata=None, method='avg'):
    # not used --pavle
        """ 
            Map gridded data onto a vector of length (nparams,)

           :param griddata: a gridded dataset to use. This dataset is mapped onto a vector of length `nparams`
           :param method: a string that specifies the method to combine grid boxes in case reverse=True. Must be either ['avg','sum','minval']
           :rtype: ndarray: size (nparameters,)

           This method makes use of a dictionary that links every parameter number [1,...,nparams] to a series of gridindices. These 
           indices specify a location on a 360x180 array, stretched into a vector using `array.flat`. There are multiple ways of calling 
           this method::

               values       = self.grid2vector(griddata=mygriddeddata,method='minval') # 
                                                                    using the minimum value of all datapoints covered by that parameter index

               values       = self.grid2vector(griddata=mygriddeddata,method='avg') # 
                                                                    using the average value of all datapoints covered by that parameter index

               values       = self.grid2vector(griddata=mygriddeddata,method='sum') # 
                                                                    using the sum of values of all datapoints covered by that parameter index

           .. note:: This method uses a DaSystem object that must be initialized with a proper parameter map. See :class:`~da.baseclasses.dasystem` for details

        """

        methods = ['avg', 'sum', 'minval']
        if method not in methods:
            logging.error("To put data from a map into the statevector, please specify the method to use (%s)" % methods)
            raise ValueError

        result = np.zeros((self.nparams,), float)
        for k, v in self.griddict.items():
#            print(k,k-1,result.shape, v)
            if method == "avg": 
                result[k - 1] = griddata.take(v).mean()
            elif method == "sum" : 
                result[k - 1] = griddata.take(v).sum()
            elif method == "minval" : 
                result[k - 1] = griddata.take(v).min()
        return result # Note that the result is returned, but not yet placed in the member.param_values attrtibute!


    def vector2grid(self, vectordata=None):
        """ 
            Map vector elements to a map or vice cersa

           :param vectordata: a vector dataset to use in case `reverse = False`. This dataset is mapped onto a 1x1 grid and must be of length `nparams`
           :rtype: ndarray: an array of size (360,180,) 

           This method makes use of a dictionary that links every parameter number [1,...,nparams] to a series of gridindices. These 
           indices specify a location on a 360x180 array, stretched into a vector using `array.flat`. There are multiple ways of calling 
           this method::

               griddedarray = self.vector2grid(vectordata=param_values) # simply puts the param_values onto a (180,360,) array

           .. note:: This method uses a DaSystem object that must be initialzied with a proper parameter map. See :class:`~da.baseclasses.dasystem` for details

        """
        result = np.zeros(self.gridmap.shape, float)
        for k, v in self.griddict.items():
#            print(k,v)
#            if k<=10:
            if k<=11:
                result.put(v, vectordata[k - 1])
        return result         

    def vector2tc(self, vectordata, cov=False):
        """ 
            project Vector onto TransCom regions 

           :param vectordata: a vector dataset to use, must be of length `nparams`
           :param cov: a Boolean to specify whether the input dataset is a vector (mean), or a matrix (covariance)
           :rtype: ndarray: an array of size (23,) (cov:F) or of size (23,23,) (cov:T)
        """

        M = self.tcmatrix
        if cov:
            return np.dot(np.transpose(M), np.dot(vectordata, M))
        else:
            return np.dot(vectordata.squeeze(), M)

    def state_to_grid(self, fluxvector=None, lag=1):
        """ 
            Transforms the StateVector information (mean + covariance) to a 1x1 degree grid.

            :param: fluxvector: a vector of length (nparams,) that holds the fluxes associated with each parameter in the StateVector
            :param: lag: the lag at which to evaluate the StateVector
            :rtype: a tuple of two arrays (gridmean,gridvariance) with dimensions (180,360,)

            If the attribute `fluxvector` is not passed, the function will return the mean parameter value and its variance on a 1x1 map.
            
            ..note:: Although we can return the variance information for each gridbox, the covariance information contained in the original ensemble is lost when mapping to 1x1 degree!

        """

        if fluxvector == None:
            fluxvector = np.ones(self.nparams)

        ensemble = self.ensemble_members[lag - 1]
        ensemblemean = ensemble[0].param_values

        # First transform the mean
        gridmean = self.vector2grid(vectordata=ensemblemean * fluxvector)

        # And now the covariance, first create covariance matrix (!), and then multiply
        deviations = np.array([mem.param_values * fluxvector - ensemblemean for mem in ensemble])
        ensemble = []
        for mem in deviations:
            ensemble.append(self.vector2grid(mem))

        return (gridmean, np.array(ensemble))

    def state2tc(self, fluxvector=None, lag=1):
        """ 
            Transforms the StateVector information (mean + covariance) to the TransCom regions.

            :param: fluxvector: a vector of length (nparams,) that holds the fluxes associated with each parameter in the StateVector
            :param: lag: the lag at which to evaluate the StateVector
            :rtype: a tuple of two arrays (mean,covariance) with dimensions ((23,), (23,23,) )

        """
        ensemble = self.ensemble_members[lag - 1]
        ensemblemean = ensemble[0].param_values

        # First transform the mean

        mean = self.vector2tc(vectordata=ensemble[0].param_values * fluxvector)

        # And now the covariance, first create covariance matrix (!), and then multiply

        deviations = np.array([mem.param_values * fluxvector - ensemblemean for mem in ensemble])
        covariance = np.dot(np.transpose(deviations), deviations) / (self.nmembers - 1)
        cov = self.vector2tc(covariance, cov=True)

        return (mean, cov)

    def get_covariance(self, date, cycleparams):
        pass
    
################### End Class StateVector ###################

if __name__ == "__main__":
    pass

