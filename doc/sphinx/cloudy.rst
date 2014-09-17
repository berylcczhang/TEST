.. highlight:: rest

.. _sec-cloudy-slug:

cloudy-slug: An Automated Interface to cloudy
=============================================

SLUG stochastically generates stellar spectra, but it does not compute
the nebular lines produced when those photons interact with the
interstellar medium. To perform such calculations, SLUG includes an
automated interface to `cloudy <http://nublado.org/>`_ (`Ferland et
al., 2013, RMxAA, 49, 137
<http://adsabs.harvard.edu/abs/2013RMxAA..49..137F>`_). This can be
used to post-process the output of a SLUG run in order to compute
nebular emission.

cloudy-slug Basics
------------------

The basic steps (described in greater detail below) are as follows:

1. Get cloudy installed and compiled, following the directions on the
   `cloudy website <http://nublado.org/>`_.

2. Set the environment variable ``$CLOUDY_DIR`` to the directory where
   the cloudy executable ``cloudy.exe`` is located.  If you are using
   a ``bash``-like shell, the syntax for this is::

      export CLOUDY_DIR = /path/to/cloudy

   while for a ``csh``-like shell, it is::

      setenv CLOUDY_DIR /path/to/cloudy

3. If you desire, edit the cloudy input template
   ``cloudy_slug/cloudy.in_template`` and the line list
   ``cloudy_slug/LineList_HII.dat``. There are the template input files
   that will be used for all the cloudy runs, and their syntax follows
   the standard cloudy syntax. They control things like the density and
   element abundances in the nebula -- see :ref:`ssec-cloudy-template`
   for more details.

4. Invoke the cloudy-slug interface script via::

     python cloudy_slug/cloudy_slug.py SLUG_MODEL_NAME

   where ``SLUG_MODEL_NAME`` is the name of the SLUG run to be
   processed. See :ref:`ssec-cloudy-slug-options` for more details on
   the python script.

5. The output will be stored as a series of additional output files of
   with names of the form SLUG_MODEL_NAME_*cloudy*.ext, where the
   extension is .txt, .bin, or .fits, depending on the format in which
   the orignal SLUG output was stored. These files can be processed
   automatically by the slugpy helper routines (see
   :ref:`sec-slugpy`). See :ref:`ssec-cloudy-output` for a full
   description of the outputs.

.. _ssec-cloudy-template:

The cloudy-slug Input Template
-------------------------

.. _ssec-cloudy-slug-options:

The cloudy-slug Interface Script
--------------------------------

.. _ssec-cloudy-output:

Full Description of cloudy-slug Output
--------------------------------------
