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

4. Perform the desired SLUG simulation. The SLUG simulation outputs
   must include spectra and photometry, and one of the photometric
   bands output must be ``QH0`` (see :ref:`ssec-phot-keywords`). If
   running in integrated mode (the default -- see
   :ref:`ssec-cloudy-cluster`), integrated specta and photometry are
   required, and if running in cluster mode, cluster spectra and
   photometry are required.

5. Invoke the cloudy-slug interface script via::

     python cloudy_slug/cloudy_slug.py SLUG_MODEL_NAME

   where ``SLUG_MODEL_NAME`` is the name of the SLUG run to be
   processed. See :ref:`ssec-cloudy-cluster` for more information on
   the underlying physical model assumed in the calculation, and
   :ref:`ssec-cloudy-slug-options` for more details on the python
   script and its options.

6. The output will be stored as a series of additional output files of
   with names of the form SLUG_MODEL_NAME_*cloudy*.ext, where the
   extension is .txt, .bin, or .fits, depending on the format in which
   the orignal SLUG output was stored. These files can be processed
   automatically by the slugpy helper routines (see
   :ref:`sec-slugpy`). See :ref:`ssec-cloudy-output` for a full
   description of the outputs.

.. _ssec-cloudy-cluster:

The cloudy-slug Physical Model: Integrated Mode Versus Cluster Mode
-------------------------------------------------------------------

Associating nebular emission with the stellar populations produced by
SLUG requires some assumptions about geometry, and some choices about
what quantities one is interested in computed. SLUG outputs both
integrated spectra for all the stars in a galaxy, and spectra for
individual clusters. One on hand, one could make the exteme assumption
that all the star clusters are spatially close enough to one another
that one can think of the entire galaxy as a single giant HII region,
and compute the nebular emission for the galaxy as a whole. This may
be a reasonable assumption for galaxies where the star formation is
highly spatially-concentrated. At the other extreme, one may assume that
there is no overlap whatsoever between the HII regions surrounding
different star clusters, so that nebular emission should be computed
for each one independently. This may be a reasonable assumption for
extended, slowly star-forming systems like the outer disk of the Milky
Way. Each of these assumptions entails somewhat different choices
about how to set the inner radius and inner density the HII region, as
required by cloudy. The cloudy-slug interface can compute nebular
emission under either of these scenarios; we refer to the former as
integrated mode, and to the latter as cluster mode.

Integrated Mode
^^^^^^^^^^^^^^^

In integrated mode, cloudy-slug will read all the spectra contained in
the SLUG integrated_spec output file, and for each stellar spectrum it
will perform a cloudy run to produce a calculation of the nebular
emission produced by that stellar spectrum interacting with a
surrounding HII region. The density in the first zone of the HII
region will be as specified by the standard ``hden`` keyword in the
cloudy input template (see :ref:`ssec-cloudy-template`). The inner
radius of the HII region will be computed automatically, and will be
set to :math:`10^{-3}` of the Stromgren radius for that density, where

.. math:: r_{\mathrm{St}} = \left(\frac{3.0 Q(\mathrm{H}^0)}{4\pi
	  \alpha_B n_{\mathrm{H}}^2}\right)^{1/3}

where :math:`Q(\mathrm{H}^0)` is the ionizing luminosity computed by
SLUG, :math:`n_{\mathrm{H}}` is the hydrogen number density stored in
the cloudy input template, and :math:`\alpha_B` is the case B
recombination coefficient, which is taken to have a value of
:math:`2.59\times 10^{-13}\;\mathrm{cm}^3\;\mathrm{s}^{-1}`.


Cluster Mode
^^^^^^^^^^^^

In cluster mode, cloudy-slug will read all the individual cluster
spectra contained in the SLUG cluster_spec file, and for each one it
will perform a cloudy calculation to determine the corresponding
nebular emission. The density and radius are handled somewhat
differently in this case, since, for a mono-age stellar population, it
is possible to compute the time evolution of the HII region radius and
density.

In cluster mode, the hydrogen number density :math:`n_{\mathrm{H}}`
stored in the cloudy input template (see :ref:`ssec-cloudy-template`)
is taken to specify the density of the *neutral* gas around the HII
region, not the density of the gas inside the HII region. The outer
radius of the HII region is then computed using the approximate
analytic solution for the expansion of an HII region into a uniform
medium, including the effects of radiation presssure and stellar wind
momentum deposition, given by `Krumholz & Matzner (2009, ApJ,
703, 1352) <http://adsabs.harvard.edu/abs/2009ApJ...703.1352K>`_. The
radius is computed from the ionizing luminosity
:math:`Q(\mathrm{H}^0)`, hydrogen number density
:math:`n_{\mathrm{H}}`, and star cluster age :math:`t` as

.. math::

   r_{\mathrm{II}} & = r_{\mathrm{ch}}
   \left(x_{\mathrm{II,rad}}^{7/2} +
   x_{\mathrm{II,gas}}^{7/2}\right)^{2/7} \\

   x_{\mathrm{II,rad}} &= (2\tau^2)^{1/4} \\

   x_{\mathrm{II,gas}} &= (49\tau^2/36)^{2/7} \\

   \tau &= t/t_{\mathrm{ch}} \\

   r_{\mathrm{ch}} & = \frac{\alpha_B}{12\pi\phi}
   \left(\frac{\epsilon_0}{2.2 k_B T_{\mathrm{II}}}\right)^2
   f_{\mathrm{trap}}^2 \frac{\psi^2 Q(\mathrm{H}^0)}{c^2} \\

   t_{\mathrm{ch}} & = \left(\frac{4\pi \mu m_{\mathrm{H}}
   n_{\mathrm{H}} c r_{\mathrm{ch}}^4}{3 f_{\mathrm{trap}}
   Q(\mathrm{H}^0) \psi \epsilon_0}\right)^{1/2}

where :math:`\alpha_B = 2.59\times
10^{-13}\;\mathrm{cm}^3\;\mathrm{s}^{-1}` is the case B recombination
coefficient, :math:`\phi = 0.73` is the fraction of ionizing photons absorbed
by hydrogen atoms rather than dust, :math:`epsilon_0 =
13.6\;\mathrm{eV}` is the hydrogen ionization potential,
:math:`T_{\mathrm{II}} = 10^4\;\mathrm{K}` is the temperature inside
the HII region, :math:`f_{\mathrm{trap}} = 2` is the trapping factor
that accounts for stellar wind and trapped infrared radiation
pressure, :math:`\psi = 3.2` is the mean photon energy in Rydberg for
a fully sampled IMF at zero age, and :math:`\mu = 1.33` is the mean
mass per hydrogen nucleus for gas of the standard cosmic
composition. See `Krumholz & Matzner (2009)
<http://adsabs.harvard.edu/abs/2009ApJ...703.1352K>`_ for a discussion
of the fiducial choices of these factors.

Once the outer radius is known, cloudy-slug sets the starting radius
for the cloudy calculation to :math:`10^{-3} r_{\mathrm{II}}`, and
sets the starting density to the value expected for photoionization
equilibrum in a uniform HII region,

.. math:: n_{\mathrm{II}} = \left(\frac{3
	  Q(\mathrm{H}^0)}{4\pi\alpha_B
	  r_{\mathrm{II}}^3}\right)^{1/2}

Note that this approximation will be highly inaccurate if
:math:`r_{\mathrm{II}} \ll r_{\mathrm{ch}}`, but no better analytic
approximation is available, and this phase should be very short-lived
for most clusters.


.. _ssec-cloudy-template:

The cloudy-slug Input Template
------------------------------

.. _ssec-cloudy-slug-options:

The cloudy-slug Interface Script
--------------------------------

.. _ssec-cloudy-output:

Full Description of cloudy-slug Output
--------------------------------------
