.. _tut_chapter2:

Chapter 2: Running your first experiment
----------------------------------------------------

**In which you will run a small sample CTDAS experiment**

Configure a small test run
^^^^^^^^^^^^^^^^^^^^^^^^^^

Modify the primary da.rc file to: ::  

    time.start   : 2000-01-01 00:00:00
    time.finish  : 2000-01-03 00:00:00
    time.cycle   : 1
    time.nlag    : 3

If these settings do not mean anything to you yet, it is probably good to go back and read again how CTDAS works.

Next, modify the das.py script to: ::

    (a) have the TM5ObservationOperator object point to the correct (compiled and tested) tm5.rc file
    (b) Initialize to the correct DaSystem rc-file + object
    (c) Initialize to the correct PlatForm object

I refer to the tutorial chapter 2 is you need a reminder on how to do this.


