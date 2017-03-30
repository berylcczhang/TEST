.. highlight:: rest

Running a SLUG simulation
=========================

Basic Serial Runs
-----------------

Once SLUG is compiled, running a simulation is extremely simple. The first step, which is not required but makes life a lot simpler, is to set the environment variable ``SLUG_DIR`` to the directory where you have installed SLUG. If you are using a ``bash``-like shell, the syntax for this is::

   export SLUG_DIR = /path/to/slug

while for a ``csh``-like shell, it is::

   setenv SLUG_DIR /path/to/slug

This is helpful because SLUG needs a lot of input data, and if you don't set this variable, you will have to manually specify where to find it.

Next, to run on a single processor, just do::

   ./bin/slug param/filename.param

where ``filename.param`` is the name of a parameter file, formatted as specified in :ref:`sec-parameters`. The code will write a series of output files as described in :ref:`sec-output`.

Thread-Based Parallelism
------------------------

If you have more than one core at your disposal, you can also run SLUG in parallel using threads, via the command line::

   python ./bin/slug.py param/filename.param

This called a python script that automatically divides up the Monte Carlo trials you have requested between the available processors, then consolidates the output so that it looks the same as if you had run a single-processor job. The python script allows fairly fine-grained control of the parallelism. It accepts the following command line arguments (not an exhaustive list -- do ``python ./bin/slug.py --help`` for the full list):

* ``-n NPROC, --nproc NPROC``: this parameter specifies the number of simultaneous SLUG processes to run. It defaults to the number of cores present on the machine where the code is running.
* ``-b BATCHSIZE, --batchsize BATCHSIZE``: this specifies how to many trials to do per SLUG process. It defaults to the total number of trials requested divided by the total number of processes, rounded up, so that only one SLUG process is run per processor. *Rationale*: The default behavior is optimal from the standpoint of minimizing the overhead associated with reading data from disk, etc. However, if you are doing a very large number of runs that are going to require hours, days, or weeks to complete, and you probably want the code to checkpoint along the way. In that case it is probably wise to set this to a value smaller than the default in order to force output to be dumped periodically.
* ``-nc, --noconsolidate``: by default the ``slug.py`` script will take all the outputs produced by the parallel runs and consolidate them into single output files, matching what would have been produced had the code been run in serial mode. If set, this flag suppresses that behavior, and instead leaves the output as a series of files whose root names match the model name given in the parameter file, plus the extension ``_pPPPPP_nNNNNN``, where the digits ``PPPPP`` give the number of the processor that produces that file, and the digits ``NNNNN`` give the run number on that processor. *Rationale*: normally consolidation is convenient. However, if the output is very large, this may produce undesirably bulky files. Furthermore, if one is doing a very large number of simulations over an extended period, and the ``slug.py`` script is going to be run multiple times (e.g., due to wall clock limits on a cluster), it may be preferable to leave the files unconsolidated until all runs have been completed.

MPI-Based Parallelism
---------------------

SLUG can also run in parallel on distributed-memory architectures using MPI. To use MPI, you must first compile the code with MPI support -- see :ref:`ssec-compiling`. Then to start an MPI-parallel computation, do::

  mpirun -np N bin/slug param/filename.param

where `N` is the number of parallel processes to run. In this mode each MPI process will write its own output files, which will be named as `MODELNAME_XXXX_FILETYPE.EXT` where `MODELNAME` is the model name specified in the parameter file (see :ref:`sec-parameters`), `XXXX` is the process number of the process that wrote the file, `FILETYPE` is the type of output file (see :ref:`sec-output`), and `EXT` is the extension specifying the file format (see :ref:`sec-output`).

If it is desirable to do so, the output files produced by an MPI run can be combined into a single output file using the ``consolidate.py`` script in the ``tools`` subdirectory.

Note that full parallel computation is only available under MPI implementations that support the MPI 3.0 standard or later. Earlier versions of MPI allow MPI functionality for SLUG in library mode (see :ref:`sec-library-mode`), but do not allow MPI parallel runs of the slug executable.


Checkpointing and Restarting
----------------------------

When running a large number of trials, it is often desirable to checkpoint the calculation, i.e., to write intermediate outputs rather than waiting until the entire calculation is done to write. SLUG can checkpoint after a specified number of trials; this number is controlled by the `checkpoint_interval` parameter (see :ref:`sec-parameters`). Checkpoint files are are named as `MODELNAME_chkYYYY_FILETYPE.EXT` (or `MODELNAME_XXXX_chkYYYY_FILETYPE.EXT` for MPI runs) where `YYYY` is the number of the checkpoint, starting at 0. Checkpoints are valid output files with some added information -- see :ref:`ssec-checkpoint-files` for details.

To restart a run from checkpoints, just give the command line option `--restart`, for example::

    mpirun -np N bin/slug param/filename.param --restart

SLUG will automatically search for checkpoint files (using the file names specified in `filename.param`), determine how many trials they contain, and resume the run to complete any remaining trials neede to reach the target number specified in the parameter file.

As with MPI runs, the output checkpoint files run can be combined into a single output file using the ``consolidate.py`` script in the ``tools`` subdirectory.
