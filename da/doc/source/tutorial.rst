.. _tutorial:

Tutorial
======================

The CTDAS tutorial is written to provide some guidance for common tasks when running, extending, or modifying CTDAS for your own purpopes. 

.. warning:: It is not a basic course in data assimilation techniques, nor in the use of python, object-oriented programming, or UNIX. 

The descriptions assume that you are familiar with these. It also assumes that you have successfully compiled and 
run a transport model (usually TM5), that can serve as observation operator for your CTDAS runs. Instructions 
on how to obtain, compile, or write such a project are not included here. We refer to the :ref:`installation <installing>` section for 
further help with those steps. 

.. toctree::
   :maxdepth: 1

   Chapter 0: Configuring TM5 to run as transport model<tut_chapter0>
   Chapter 1: Configuring CTDAS to run an experiment <tut_chapter1>
   Chapter 2: Running your first experiment  <tut_chapter2>
   Chapter 3: Controlling your experiment: Cold start, Restart, or Recover from crash<tut_chapter3>
   Chapter 4: Adding more observations <tut_chapter4>
   Chapter 5: Modifying the state vector<tut_chapter5>
   Chapter 6: Changing the covariance structure<tut_chapter6>
   Chapter 7: Adding a new type of observations <tut_chapter7>



