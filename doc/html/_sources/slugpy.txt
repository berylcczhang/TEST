.. highlight:: rest

.. _sec-slugpy:

slugpy -- The Python Helper Library
===================================

Basic Usage
-----------

SLUG comes with the python module slugpy, which contains an extensive set of routines for reading, writing, and manipulating SLUG outputs. The most common task is to read a set of SLUG outputs into memory so that they can be processed. To read the data from a SLUG run using slugpy, one can simply do the following::

   from slugpy import *
   idata = read_integrated('SLUG_MODEL_NAME')
   cdata = read_cluster('SLUG_MODEL_NAME')

The ``read_integrated`` function reads all the integrated-light data (i.e., the data stored in the ``_integrated_*`` files -- see :ref:`sec-output`) for a SLUG output whose name is given as the argument. This is the base name specified by the ``model_name`` keyword (see :ref:`ssec-basic-keywords`), without any extensions; the slugpy library will automatically determine which outputs are available and in what format, and read the appropriate files. It returns a ``namedtuple`` containing all the output data available for that simulation. The ``read_cluster`` function is analogous, except that instead of reading the whole-galaxy data, it reads data on the individual star clusters, as stored in the ``_cluster_*`` output files.

Full Documentation
------------------

.. automodule:: slugpy
   :members:
