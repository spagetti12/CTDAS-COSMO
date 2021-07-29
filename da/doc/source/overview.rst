.. _overview:

.. index:: overview, philosophy, development centers

Overview
========

This page gives a brief overview of CTDAS, with the intention to provide background information for (potential) new users. Afer reading the information below we recommend contacting one of the CTDAS development centers (Wageningen University and NOAA ESRL), or proceed to the :doc:`tutorial`.

What is CTDAS?
--------------

CTDAS is short for "CarbonTracker Data Assimilation Shell".  
This is the implementation of an extendible data assimilation framework for CarbonTracker,
developed by NOAA ESRL & Wageningen University, in close cooperation with many partners around the world.

The aim of the CTDAS system is to facilitate the use of CarbonTracker, and to 
foster its development by its many international partners. 


How can you use CTDAS?
----------------------
At its most basic, the CTDAS is a simple 
control system for CarbonTracker that deals with the running, optimization, analysis, and time stepping of the system. 
For advanced users, the CTDAS provides an easy way to extend or modify the system by introducing new configurations 
for the state vector, for the observations, or even for the transport model or optimization method.

CTDAS Design philosophy
-----------------------

CTDAS is implemented in python. Each component of the data assimilation system is written as a separate python class, which are 
combined in a pipeline. Separate pipelines exist for an inverse simulation, or a forward run. Four classes compose the core of 
the system: Observations, StateVector, ObservationOperator, and Optimizer. These are controlled by a CycleControl object, which 
itself holds information from class PlatForm and class DaSystem. For each of these seven components, a "baseclass" exists that 
describes the basic layout, and mandatory methods that each class needs to have. A specific implementation of a baseclass 
"inherits" this basic behavior, and then extends or overwrites the methods. 

As a typical example, in the module :mod:`~da.baseclasses.observationoperator`, a baseclass :class:`~da.baseclasses.observationoperator.ObservationOperator` is a nearly empty object that contains methods to :meth:`~da.baseclasses.observationoperator.ObservationOperator.Initialize()` the object, to :meth:`~da.baseclasses.observationoperator.ObservationOperator.Validate()` the input for it, to :meth:`~da.baseclasses.observationoperator.ObservationOperator.Run()` the operator, and to :meth:`~da.baseclasses.observationoperator.ObservationOperator.SaveData()` that is needed for a next cycle. The class :class:`~da.tm5.observationoperator.TM5ObservationOperator` is derived from this baseclass, and has a new implementation of each of these methods. When they are called in the pipeline, a TM5 model run is prepared and executed by the methods of the class, returning a set of samples for the optimizer to use in the 
minimum least squares method.

This design philosophy makes the pipeline completely independent of the specific implementation. It simply calls methods from the clasees it receives in a given order, and ensures that the flow of time and data works. The pipeline does not know (or need to know) whether it is running a gridded inversion with satellite data and the TM5 model without zoom, or whether it is running a methane inversion with TM3. The implementation is thus in the hands of the user.

Extending CTDAS
---------------
To extend CTDAS, one needs to write new classes that inherit from one of the baseclasses, and have specific functionality implemented under the methods called from the pipeline. This can be a very simple task, or a very hard task depending on which functionality you want to add, or alter.

For instance, to make a new StateVector for your assimilation system with a different number of unknowns and different covariance between them, it would suffice to make a new class "MyStateVector" that inherits everything from the :class:`~da.baseclasses.statevector.StateVector` class, and then write one specific method to replace the standard :meth:`~da.baseclasses.statevector.StateVector.GetCovariance()` method.

For more specific instructions and examples, see the :ref:`Tutorial`.



