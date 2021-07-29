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

            if member == 1 and lag == 0:
                rands=np.array([-0.6865332,0.3558595,0.7385287,0.8752266,-0.6774718,-0.1993868,0.8089826,0.754124,0.8117286,-0.2621968,0,-0.277585,-0.2579767,0.1577864,0.2642049,0.901803,1,-0.2622429,0.4879664,0.2333322,0.1841461,0,0.005501025])
            if member == 2 and lag == 0:
                rands=np.array([0.854864,0.6554427,0.4478732,0.9646293,-0.5737644,-0.568122,-0.4400042,0.1481517,0.3478025,0.7362454,0,-0.8154956,-0.4146178,0.4900621,0.3009739,-0.3057978,-0.293597,0.09253449,0.08424578,0.8465578,0.9420493,0,-0.0219593])
            if member == 3 and lag == 0:
                rands=np.array([-0.6324199,-0.002833775,0.1350392,0.3754209,-0.06226669,0.3743071,0.0659928,-0.7189865,0.3424602,0.397964,0,-0.7215626,-0.007501739,0.3145521,0.3148233,0.1010255,0.6437694,0.9282547,0.03162892,1,0.447092,0,0.04623652])
            if member == 4 and lag == 0:
                rands=np.array([-0.04206607,-0.2271331,-0.3525559,-0.805801,0.4870083,0.3036494,0.5990416,-0.5578776,-0.4825368,0.5263366,0,-0.58166,-1,0.2483742,-0.1055822,0.5741059,0.8543911,-0.3701363,0.3532563,0.4482012,-0.1065976,0,0.006667121])
            if member == 5 and lag == 0:
                rands=np.array([-0.7474335,-0.2877136,-0.9696097,-0.4469923,-1,0.07805356,-0.5317366,0.2719007,-0.6825466,-0.03831927,0,-0.8331428,-0.4398521,-0.4037699,-0.7918211,0.6683037,0.673253,-0.2467666,-0.3684593,-0.290229,-0.2196416,0,-0.004629975])
            if member == 6 and lag == 0:
                rands=np.array([0.6362485,1,-0.3317139,-0.5523597,0.5481536,-0.157634,0.7913657,0.1318698,-0.2012699,0.4593242,0,0.8619153,0.5313798,-0.2796631,0.2096716,0.4844697,1,0.7664097,-0.5802394,-0.2993027,0.8965644,0,-0.01908342])
            if member == 7 and lag == 0:
                rands=np.array([0.7788572,-0.394289,0.8726047,0.3016202,0.8938656,-0.0611987,-0.503676,-0.8730192,0.552582,0.4025556,0,-0.5100973,-0.9004463,0.6479801,-0.7253889,-0.05091742,0.7477584,0.05360383,0.3514898,-0.1618125,-0.1032103,0,0.00449502])
            if member == 8 and lag == 0:
                rands=np.array([-0.385563,-0.7897603,-1,-0.1659514,-0.04117743,-1,0.4338062,-0.3634367,0.5738267,-0.2297338,0,-0.2141016,0.1454927,0.7582307,-0.1823771,0.4570185,0.5109985,-0.8624253,-0.9817374,0.3879869,0.3235793,0,0.007567421])
            if member == 9 and lag == 0:
                rands=np.array([0.02401263,-0.8005567,-0.2080475,0.518029,-0.4546398,0.740793,0.4948097,1,0.6126591,-0.06489301,0,-0.3209458,-0.9661453,-0.5611605,-0.1535542,0.2005962,-0.3554438,0.5859257,0.1857199,-0.3683482,-0.6147193,0,-0.01431568])
            if member == 10 and lag == 0:
                rands=np.array([-0.6043704,0.663038,-0.7966857,0.09454904,0.001207841,-0.1590567,0.1638788,-0.6406232,-0.6433946,0.4149962,0,0.0364261,0.7709227,1,0.4498931,0.005313124,-0.6495739,-0.05076865,0.1495458,0.7621301,0.9752402,0,-0.04433159])
            if member == 11 and lag == 0:
                rands=np.array([-0.8747755,-0.8251697,-0.7172815,0.3534077,0.4100716,-0.2430783,-0.8679456,0.3135208,0.7610837,0.6283105,0,-0.6373556,-0.09406921,-0.6684051,0.3434846,-0.6966717,-1,0.4377855,0.8649635,0.07060629,0.2723278,0,-0.02111423])
            if member == 12 and lag == 0:
                rands=np.array([-0.3798023,-0.08313549,-0.3811138,-0.1380353,0.6731774,-0.1338478,0.4247293,0.958501,0.161129,0.5277465,0,0.8514146,0.307235,-0.5815715,0.3430769,0.08776622,-0.3751208,0.6033686,-0.8277515,0.2015533,-0.6236387,0,-0.03661139])
            if member == 13 and lag == 0:
                rands=np.array([0.002277945,0.6822468,0.4901579,-0.4742505,-0.004097998,0.1202785,0.5017828,-0.4101433,-0.5463102,-0.8954731,0,-0.9228076,0.5828303,0.2402924,-0.4813775,0.1176145,-0.7799372,0.5447594,-0.3260884,0.8725086,-0.3888961,0,0.01957087])
            if member == 14 and lag == 0:
                rands=np.array([-0.5426614,0.3118683,-0.6139171,0.1220636,0.2437155,-0.4154192,0.9160472,0.9443241,1,1,0,-0.8581675,0.3852032,-0.6203499,-0.1785067,-0.05744135,0.374786,0.3990103,-0.1792489,-0.7705899,-0.2908233,0,0.006832954])
            if member == 15 and lag == 0:
                rands=np.array([0.9016066,-0.08090271,0.00404155,0.491021,0.9821634,0.535128,-0.5960406,-0.7466612,-0.2362669,-0.2788562,0,-0.04636844,-0.7322,-0.3879876,-0.6992238,-0.4355617,0.657073,0.7211763,0.6343221,-0.6410404,0.5267332,0,-0.003195141])
            if member == 16 and lag == 0:
                rands=np.array([-0.1508753,0.5437378,0.3320757,0.7468662,0.4373512,0.1174251,-0.9389582,-0.7015296,-0.8911921,-0.02896086,0,-0.09633388,-0.7430251,-0.3673911,-0.6709117,-0.3223665,-0.9080444,-0.7327068,-0.5885021,-0.4674551,-0.4288323,0,-0.04847308])
            if member == 17 and lag == 0:
                rands=np.array([0.4761747,0.8022127,-0.653662,-0.4827527,-0.8192345,0.008643345,-0.5563272,-0.4879796,-0.03064474,-0.5957319,0,-0.9211775,-0.7920011,-1,-0.4938926,-0.8254521,0.1540332,-0.05616575,-0.8989872,-0.191664,0.5984833,0,-0.04690293])
            if member == 18 and lag == 0:
                rands=np.array([0.5920676,-0.1176148,0.981807,0.5902563,-0.5075156,-0.7445332,0.3569477,-0.7992576,0.4009235,-0.5261548,0,-0.7176116,0.3382274,-0.3050335,-0.3559639,0.2135999,0.687503,0.1296655,0.7078268,0.0410387,0.7803431,0,-0.02164244])
            if member == 19 and lag == 0:
                rands=np.array([-0.936323,-0.004591943,-0.14192,0.01133621,0.8316824,0.137822,-0.2853822,0.9385482,0.5062207,0.231412,0,0.05361838,-0.7293404,0.3445812,-0.1972966,-0.1143838,0.4922296,0.1472069,-0.1203217,-0.3042997,-0.6008288,0,-0.0113157])
            if member == 20 and lag == 0:
                rands=np.array([0.9894645,-0.3027939,0.01672673,0.9668484,0.2120544,-0.3184669,-0.3186892,-0.618659,0.1277333,0.4045099,0,-0.7060359,-0.8333552,-0.4819083,0.5372346,-0.7353554,0.1417296,0.2413463,0.2322551,0.5773139,0.7776622,0,0.03813688])

            if member == 1 and lag == 1:
                rands=np.array([-0.4110718,-0.88892,0.2737297,-0.1388108,0.535023,-0.259322,0.8953033,0.6445236,-0.1747866,0.1983998,0,0.2436623,0.003266148,-0.5901974,0.6207106,-0.06291912,-0.06488661,0.9013689,0.03002323,0.7768148,0.7369732,0,0.0432417])
            if member == 2 and lag == 1:
                rands=np.array([-0.4944594,-0.9457298,-1,-0.1354843,-0.5887785,-0.15213,0.5869877,0.2237128,0.3282241,-0.7923626,0,-0.4187486,-0.9977051,0.4647258,0.4916891,-0.9991935,0.4891383,0.4701071,-0.5617338,-0.05081338,0.4587131,0,0.0129059])
            if member == 3 and lag == 1:
                rands=np.array([0.8555521,-0.3292904,0.2105805,-0.2740403,0.3676851,-0.713368,0.6796349,0.6576366,0.9874129,1,0,-0.5120251,0.08154307,0.06941926,0.6559066,-0.170521,-0.7472118,0.4263452,0.4463684,-0.8195287,-0.5793816,0,0.03172341])
            if member == 4 and lag == 1:
                rands=np.array([0.355482,0.8289347,0.01088188,-0.3384805,0.6268374,1,-0.2082716,0.6971263,-0.3328478,-0.9343829,0,0.4667932,0.02340866,-0.003623952,-0.684778,0.5278456,-0.6662421,-0.1499434,-0.4975282,0.6486844,-0.1274243,0,0.005940544])
            if member == 5 and lag == 1:
                rands=np.array([-0.8808182,0.5226131,-0.8998929,-0.2023875,-0.6950441,-0.2061953,0.367494,-0.6245959,0.9017839,-0.7293146,0,-0.2431296,0.2236341,-0.1754546,0.5178533,-0.4121083,0.8061285,0.1651406,-0.9084586,-1,0.5162783,0,0.01680873])
            if member == 6 and lag == 1:
                rands=np.array([0.06329523,-0.8564324,-0.8164288,-0.06235388,-0.142993,-0.3798361,0.7687777,-0.3145005,0.7740852,-0.6808619,0,-0.6647659,-0.683647,-0.1215333,0.7125296,-0.6953543,-0.9192593,-0.2288445,-0.2399431,0.09424525,0.2002725,0,-0.04112585])
            if member == 7 and lag == 1:
                rands=np.array([0.5613453,0.6810185,-0.05334619,0.9520684,0.5774254,-0.361872,0.5680153,0.6171691,0.3382697,1,0,0.5451796,0.7761162,0.7679793,0.2387103,0.8507853,-0.4915903,-0.1660851,-0.4182827,-0.3919242,-0.5486335,0,0.04137752])
            if member == 8 and lag == 1:
                rands=np.array([-0.3455416,-0.8222767,0.3813077,-0.2371309,0.4406765,-0.3945573,0.07800323,0.3014143,0.5139478,0.395775,0,0.3030554,0.3006979,-0.2408492,-0.5032662,-0.3779641,-0.7752672,-0.6641146,-0.3464549,-0.9769148,-0.268956,0,0.04573557])
            if member == 9 and lag == 1:
                rands=np.array([-0.6286283,-0.5469043,-0.1558352,-0.3414178,0.6410176,0.1428595,0.260674,-0.5694779,0.4645444,0.02658742,0,-0.8926882,-0.1016959,0.3635801,0.6021298,-0.5903702,-0.5853319,0.4431821,0.3508914,0.2912112,0.3294599,0,0.009273289])
            if member == 10 and lag == 1:
                rands=np.array([0.3917164,0.01293407,-0.1093186,0.9098022,-0.2894615,0.7123989,-0.8058071,0.3248312,-0.04407874,0.1965292,0,0.6192394,-0.5915784,0.7952496,0.6682949,-0.1606273,-0.518069,-0.4748175,0.1653854,-0.01001874,0.002010202,0,-0.02243915])
            if member == 11 and lag == 1:
                rands=np.array([-0.9141791,-0.241822,0.6085451,-0.3847938,-0.3815159,0.5437576,0.8222843,0.07603163,0.8769036,-0.3614713,0,-0.37404,-0.4086722,-0.1993674,-0.4439494,-0.2315528,-0.7542738,-0.2330906,0.1764065,-0.1569091,-0.3884712,0,-0.04503115])
            if member == 12 and lag == 1:
                rands=np.array([0.3941141,-0.7811217,-0.5984285,0.1827517,0.03472878,0.01964883,-0.8944655,-0.5264197,-1,-0.7233462,0,0.9766217,0.2840011,-0.5602549,-0.1843678,0.7192143,-0.1638302,-0.3837061,0.6651133,0.2903731,-0.1519475,0,-0.006677282])
            if member == 13 and lag == 1:
                rands=np.array([0.4172789,0.1763749,-0.2614514,-0.9539597,0.3494455,0.09771311,0.2078965,-0.5552518,0.8163491,0.1633967,0,0.18388,0.1524801,0.8148156,1,-0.5217537,-0.4205278,0.007124796,-0.345917,-0.3345222,-0.7838224,0,-0.02678765])
            if member == 14 and lag == 1:
                rands=np.array([0.8430449,0.1706411,0.1663785,0.2937862,0.1398123,0.9573265,0.3895729,-0.5957338,0.7086564,-0.4036448,0,0.3960939,0.5560951,-0.1537433,-0.4057478,0.1520463,0.3624675,0.1260598,0.4589263,-0.3366966,0.6250498,0,0.001293249])
            if member == 15 and lag == 1:
                rands=np.array([-0.7933657,0.03643509,0.779883,1,-0.3035239,-0.1296805,0.2410311,0.1341446,0.686308,0.4610172,0,-0.2255432,0.1047536,-0.1770354,-0.08493959,-0.1798149,0.8312726,0.5267493,0.2262887,-0.7638802,0.9652162,0,0.02293831])
            if member == 16 and lag == 1:
                rands=np.array([-0.9171996,-0.4681723,-0.3022334,-0.3627997,-0.7253119,-1,-0.2520878,-0.07770864,-0.7227243,0.6474443,0,-0.7581744,0.2004271,-0.5477204,-0.5136528,-0.8772402,0.1973221,-0.74787,-0.5328937,-0.9429494,0.5461707,0,-0.02661308])
            if member == 17 and lag == 1:
                rands=np.array([-0.3845477,-0.5994932,-0.2322402,0.1683227,-0.2870263,-0.921927,0.47649,0.101328,-0.09137523,-0.7140735,0,-0.9559138,0.4263357,-0.8882524,-0.1593087,0.5885355,0.271398,-0.6340341,0.3255026,-0.1300831,-0.5140868,0,0.01148178])
            if member == 18 and lag == 1:
                rands=np.array([0.07484807,0.4029489,0.116799,-0.02001172,-0.2764485,-0.6566548,0.5756098,-0.7976735,0.7077577,-0.2850362,0,-0.3977833,0.1738415,-0.7056463,0.4181726,-0.02520084,-0.9478775,-0.5392118,0.2711374,-0.9957126,0.597542,0,-0.008706607])
            if member == 19 and lag == 1:
                rands=np.array([-0.9447101,-1,-1,-1,0.1799199,-0.7011237,0.6601787,0.741608,0.8797131,1,0,-0.0915675,0.4198072,-0.6981921,0.6845254,-0.6768668,0.1707276,-0.2262995,0.4536436,0.400492,-0.439103,0,-6.52E-05])
            if member == 20 and lag == 1:
                rands=np.array([0.4629868,-0.6783764,-0.4005693,0.6058718,0.0389765,-0.5943372,-0.5393985,-0.8952301,0.1574879,0.3106839,0,0.5308425,0.2482044,0.8471618,-0.1107048,0.05897403,0.2191929,0.02109409,-0.5249093,0.6646884,-0.4554588,0,0.04243601])

            newmember = EnsembleMember(member)
            newmember.param_values = rands + newmean
       #####     newmember.param_values = (np.hstack((np.dot(C, rands[0:10]),np.dot(C, rands[10:20]), rands_bg)) + newmean).ravel() #not this one
#            newmember.param_values = (np.hstack((np.dot(C, rands[0:11]),np.dot(C, rands[11:22]), rands_bg)) + newmean).ravel()    # THIS ONE
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
            dimgrid = ncf.add_latlon_dim()

            data = mem.param_values[0:11]

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
            dimparams = ncf.add_params_dim(11)
            dimgrid = ncf.add_latlon_dim()

            data = mem.param_values[11:22]

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

#
    def write_members_for_cosmo_non_frac(self, lag, outdir,endswith='.nc'):
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

