.. highlight:: rest

Running a SLUG simulation
=========================

Once SLUG is compiled, running a simulation is extremely simple. The first step, which is not required but makes life a lot simpler, is to set the environment variable ``SLUG_DIR`` to the directory where you have installed SLUG. If you are using a ``bash``-like shell, the syntax for this is::

   export SLUG_DIR = /path/to/slug

while for a ``csh``-like shell, it is::

   setenv SLUG_DIR /path/to/slug

This is helpful because SLUG needs a lot of input data, and if you don't set this variable, you will have to manually specify where to find it.

Next, to run on a single processor, just do::

   ./bin/slug param/filename.param

where ``filename.param`` is the name of a parameter file, formatted as specified in :ref:`sec-parameters`. The code will write a series of output files as described in :ref:`sec-output`.

If you have more than one core at your disposal, you can also run SLUG in parallel, using the command line::

   python ./bin/slug.py param/filename.param

This called a python script that automatically divides up the Monte Carlo trials you have requested between the available processors, then consolidates the output so that it looks the same as if you had run a single-processor job. The python script allows fairly fine-grained control of the parallelism. It accepts the following command line arguments:

* ``-n NPROC, --nproc NPROC``: this parameter specifies the number of simultaneous SLUG processes to run. It defaults to the number of cores present on the machine where the code is running
* ``-b BATCHSIZE, --batchsize BATCHSIZE``: this specifies how to many trials to do per SLUG process. It defaults to the total number of trials requested divided by the total number of processes, rounded up, so that only one SLUG process is run per processor. *Rationale*: The default behavior is optimal from the standpoint of minimizing the overhead associated with reading data from disk, etc. However, if you are doing a very large number of runs that are going to require hours, days, or weeks to complete, and you probably want the code to checkpoint along the way. In that case it is probably wise to set this to a value smaller than the default in order to force output to be dumped periodically.
* ``-nc, --noconsolidate``: by default the ``slug.py`` script will take all the outputs produced by the parallel runs and consolidate them into single output files, matching what would have been produced had the code been run in serial mode. If set, this flag suppresses that behavior, and instead leaves the output as a series of files whose root names match the model name given in the parameter file, plus the extension ``_pPPPPP_nNNNNN``, where the digits ``PPPPP`` give the number of the processor that produces that file, and the digits ``NNNNN`` give the run number on that processor. *Rationale*: normally consolidation is convenient. However, if the output is very large, this may produce undesirably bulky files. Furthermore, if one is doing a very large number of simulations over an extended period, and the ``slug.py`` script is going to be run multiple times (e.g.\ due to wall clock limits on a cluster), it may be preferable to leave the files unconsolidated until all runs have been completed.

