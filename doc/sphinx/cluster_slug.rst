.. highlight:: rest

.. _sec-cluster-slug:

cluster_slug: Bayesian Inference of Star Cluster Properties
===========================================================

The slugpy.cluster_slug module computes posterior probabilities for the mass, age, and extinction of star clusters from a set of input photometry.  It is implemented as a wrapper around :ref:`sec-bayesphot`, so for details on how the calculation is performed see the bayesphot documentation.

Getting the Default Library
---------------------------

The cluster_slug module requires a pre-computed library of slug simulations to use as a "training set" for its calculations. Due to its size, the default library *is not* included in the slug git repository. Instead, it is provided for download from the `SLUG data products website <http://www.slugsps.com/data>`_. Download the two files ``clusterslug_mw_cluster_phot.fits`` and ``clusterslug_mw_cluster_prop.fits`` and save them in the ``cluster_slug`` directory of the main respository. If you do not do so, and do not provide your own library when you attempt to use cluster_slug, you will be prompted to download the default library.


Basic Usage
-----------

For an example of how to use cluster_slug, see the file ``cluster_slug/cluster_slug_example.py`` in the repository. All funtionality is provided through the cluster_slug class. The basic steps are as follows:

1. Import the library and instantiate an ``sfr_slug`` object (see :ref:`ssec-cluster-slug-full` for full details)::

     from slugpy.cluster_slug import cluster_slug
     cs = cluster_slug(photsystem=photsystem)

This creates a cluster_slug object, using the default simulation library. If you have another library of simulations you'd rather use, you can use the ``libname`` keyword to the ``cluster_slug`` constructor to select it. The optional argument ``photsystem`` specifies the photometric system you will be using for your data. Possible values are ``L_nu`` (flux per unit frequency, in erg/s/Hz), ``L_lambda`` (flux per unit wavelength, in erg/s/Angstrom), ``AB`` (AB magnitudes), ``STMAG`` (ST magnitudes), and ``Vega`` (Vega magnitudes). If left unspecified, the photometric system will be whatever the library was written in; the default library is in the ``L_nu`` system. Finally, if you have already read a library into memory using ``read_cluster``, you can set the keyword ``lib`` in the ``cluster_slug`` constructor to specify that library should be used.

2. Specify your filter(s), for example::

     cs.add_filters(['WFC3_UVIS_F336W', 'WFC3_UVIS_F438W', 'WFC3_UVIS_F555W',
                     'WFC3_UVIS_F814W', 'WFC3_UVIS_F657N'])

The ``add_filter`` method takes as an argument a string or list of strings specifying which filters were used for the observations you're going to analyze. You can have more than one set of filters active at a time (just by calling ``add_filters`` more than once), and then specify which set of filters you're using for any given calculation.

3. Specify your priors, for example::

     # Set priors to be flat in log T and A_V, but vary with log M as
     # p(log M) ~ 1/M
     def priorfunc(physprop):
        # Note: physprop is an array of shape (N, 3) where physprop[:,0] =
	# log M, physprop[:,1] = log T, physprop[:,2] = A_V
        return 1.0/exp(physprop[:,0])
     cs.priors = prorfunc

The ``priors`` property specifies the assumed prior probability distribution on the physical properties of star clusters. It can be either ``None`` (in which case all simulations in the library are given equal prior probability), an array with as many elements as there are simulations in the library giving the prior for each one, or a callable that takes a vector of physical properties as input and returns the prior for it.

4. Generate a marginal posterior probability distribuiton via::

     logm, pdf = cs.mpdf(idx, phot, photerr = photerr)

The first argument ``idx`` is an index for which posterior distribution should be computed -- a value of 0 generates the posterior in log mass, a value of 1 generates the posterion on log age, and a value of generates the posterior in A_V. The second argument ``phot`` is an array giving the photometric values in the filters specified in step 2; make sure you're using the same photometric system you used in step 1. For the array ``phot``, the trailing dimension must match the number of filters, and the marginal posterior-finding exercise is repeated over every value in the leading dimensions. If you have added two or more filter sets, you need to specify which one you want to use via the ``filters`` keyword. The optional argument ``photerr`` can be used to provide errors on the photometric values. The shape rules on it are the same as on ``phot``, and the two leading dimensions of the two arrays will be broadcast together using normal broadcasting rules.

The ``cluster_slug.mpdf`` method returns a tuple of two quantities. The first is a grid of values for log M, log T, or A_V, depending on the value of ``idx``. The second is the posterior probability distribution at each value of of the grid. Posteriors are normalized to have unit integral. If the input consisted of multiple sets of photometric values, the output will contains marginal posterior probabilities for each input. The output grid will be created automatically be default, but all aspects of it (shape, size, placement of grid points) can be controlled by keywords -- see :ref:`ssec-cluster-slug-full`.


Using cluster_slug in Parallel
------------------------------

The ``cluster_slug`` module has full support for threaded computation using the python `multiprocessing module <https://docs.python.org/2.7/library/multiprocessing.html>`_. This allows efficient use of multiple cores on a shared memory machine, without the need for every project to read a large simulation library or store it in memory. See :ref:`ssec-bayesphot-threading` for full details on the recommended paradigm for parallel computing. The full list of thread-safe ``cluster_slug`` methods is:

* ``cluster_slug.logL``
* ``cluster_slug.mpdf``
* ``cluster_slug.mcmc``
* ``cluster_slug.bestmatch``
* ``cluster_slug.make_approx_phot``
* ``cluster_slug.make_approx_phys``
* ``cluster_slug.squeeze_rep``
* ``cluster_slug.mpdf_approx``


Making Your Own Library
-----------------------

You can generate your own library by running slug; you might want to do this, for example, to have a library that works at different metallicity or for a different set of stellar tracks. An example parameter file (the one that was used to generate the default clusterslug_mw library) is included in the ``cluster_slug`` directory. This file uses slug's capability to pick the output time and the cluster mass from specified PDFs.

One subtle post-processing step you should take once you've generated your library is to read it in using :ref:`sec-slugpy` and then write the photometry back out using the ``slugpy.write_cluster_phot`` routine with the format set to ``fits2``. This uses an alternative FITS format that is faster to search when you want to load only a few filters out of a large library. For large data sets, this can reduce cluster_slug load times by an order of magnitude. (To be precise: the default format for FITS outputs to put all filters into a single binary table HDU, while the ``fits2`` format puts each filter in its own HDU. This puts all data for a single filter into a contiguous block, rather than all the data for a single cluster into a contiguous block, and is therefore faster to load when one wants to load the data filter by filter.)


Variable Mode IMF
-----------------

If your library was run with variable IMF parameters, these can also be used in ``cluster_slug``. When creating a ``cluster_slug`` object, you can pass the array ``vp_list`` as an argument. This list should have an element for each variable parameter in your library. Each element should then be either ``True`` or ``False`` depending on whether you wish to include this parameter in the analysis.
For example, for a library with four variable parameters you could have::
      
      vp_list=[True,False,True,True]


.. _ssec-cluster-slug-full:

Full Documentation of slugpy.cluster_slug
-----------------------------------------

.. autoclass:: slugpy.cluster_slug.cluster_slug
   :members:
   :special-members:
