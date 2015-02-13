.. highlight:: rest

.. _sec-sfr-slug:

sfr_slug: Bayesian Inference of Star Formation Rates
====================================================

The slugy.sfr_slug module computes posterior probabilities on star formation rates given a set of star formation rates estimated using the "point mass estimate" (i.e., the estimate you would get for a fully sampled stellar population) for the SFR based on the ionizing, FUV, or bolometric luminosity. It is implemented as a wrapper around :ref:`sec-bayesphot`, so for details on how the calculation is performed see the bayesphot documentation.

Getting the Default Library
---------------------------

The sfr_slug module requires a pre-computed library of slug simulations to use as a "training set" for its calculations. Due to its size, the default library *is not* included in the slug git repository. Instead, it is provided for download from the `SLUG data products website <http://www.slugsps.com/data>`_. Download the two files ``SFR_SLUG_integrated_phot.fits`` and ``SFR_SLUG_integrated_prop.fits`` and save them in the ``sfr_slug`` directory of the main respository. If you do not do so, and do not provide your own library when you attempt to use sfr_slug, you will be prompted to download the default library.


Basic Usage
-----------

The ``sfr_slug/sfr_slug_example.py`` file in the repository provides an example of how to use sfr_slug. Usage of is simple, as the functionality is all implemented through a single class, sfr_slug. The required steps are as follows:

1. Import the library and instantiate an ``sfr_slug`` object (see :ref:`sec-sfr-slug-full` for full details)::

     from slugpy.sfr_slug import sfr_slug
     sfr_estimator = sfr_slug()

This creates an sfr_slug object, using the default simulation library, $SLUG_DIR/sfr_slug/SFR_SLUG. If you have another library of simulations you'd rather use, you can use the ``libname`` keyword to the ``sfr_slug`` constructor to select it.

2. Specify your filter(s), for example::

     sfr_estimator.add_filters('QH0')

The ``add_filter`` method takes as an argument a string or list of strings specifying which filters you're going to point mass SFRs based on. You can have more than one set of filters active at a time (just by calling ``add_filters`` more than once), and then specify which set of filters you're using for any given calculation.

3. Specify your priors, for example::

     sfr_estimator.priors = 'schechter'

The ``priors`` property specifies the assumed prior probability distribution on the star formation rate. It can be either ``None`` (in which case all simulations in the library are given equal prior probability), an array with as many elements as there are simulations in the library giving the prior for each one, a callable that takes a star formation rate as input and returns the prior for it, or a string whose value is either "flat" or "prior". The two strings specify, respectively, a prior distribution that is either flat in log SFR or follows the Schechter function SFR distribution from `Bothwell et al. (2011) <http://adsabs.harvard.edu/abs/2011MNRAS.415.1815B>`_:

.. math:: p(\log\mathrm{SFR}) \propto \mathrm{SFR}^{\alpha} \exp(-\mathrm{SFR}/\mathrm{SFR}_*)

with :math:`\alpha = -0.51` and :math:`\mathrm{SFR}_* = 9.2\,M_\odot\,\mathrm{yr}^{-1}`.

4. Generate the posterior probability distribuiton of SFR via::

     logSFR, pdf = sfr_estimator.mpdf(logSFR_in, logSFRphoterr = logSFR_err)

The argument ``logSFR_in`` can be a float or an array specifying one or more point mass estimates of the SFR in your chosen filter. For a case with two or more filters, then ``logSFR_in`` must be an array whose trailing dimension matches the number of filters. If you have added two or more filter sets, you need to specify which one you want to use via the ``filters`` keyword. The optional argument ``logSFRphoterr`` can be used to provide errors on the photometric SFRs. Like ``logSFR_in``, it can be a float or an array.

The ``sfr_slug.mpdf`` method returns a tuple of two quantities. The first is a grid of log SFR values, and the second is the posterior probability distribution at each value of log SFR. If the input consisted of multiple photometric SFRs, the output will contains posterior probabilities for each input. The output grid will be created automatically be default, but all aspects of it (shape, size, placement of grid points) can be controlled by keywords -- see :ref:`sec-sfr-slug-full`.


.. _sec-sfr-slug-full:

Full Documentation of slugpy.sfr_slug
-------------------------------------

.. autoclass:: slugpy.sfr_slug.sfr_slug
   :members:
   :special-members:
