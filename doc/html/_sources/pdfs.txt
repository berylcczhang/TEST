.. highlight:: rest

.. _sec-pdfs:

Probability Distribution Functions
==================================

The SLUG code regards the IMF, the CMF, the CLF, the SFH, and the extinction :math:`A_V` as probability distribution functions -- see :ref:`ssec-slugpdfs`. The code provides a generic file format through which PDFs can be specified. Examples can be found in the ``lib/imf``, ``lib/cmf``, ``lib/clf``, and ``lib/sfh`` directories of the SLUG distribution.

PDFs in SLUG are generically written as functions

.. math:: \frac{dp}{dx} = n_1 f_1(x; x_{1,a}, x_{1,b}) + n_2 f_2(x; x_{2,a}, x_{2,b}) + n_3 f_3(x; x_{3,a}, x_{3,b}) + \cdots,

where :math:`f_i(x; x_{i,a}, x_{i,b})` is non-zero only for :math:`x \in [x_{i,a}, x_{i,b}]`. The functions :math:`f_i` are simple continuous functional forms, which we refer to as *segments*. Functions in this form can be specified in SLUG in two ways.

Basic Mode
----------

The most common way of specifying a PDF is in basic mode. Basic mode describes a PDF that has the properties that

#. the segments are contiguous with one another, i.e., :math:`x_{i,b} = x_{i+1,a}`
#. :math:`n_i f_i(x_{i,b}; x_{i,a}, x_{i,b}) = n_{i+1} f_{i+1}(x_{i+1,a}; x_{i+1,a}, x_{i+1,b})`
#. the overall PDF is normalized such that :math:`\int (dp/dx)\, dx = 1`

Given these constraints, the PDF can be specified fully simply by giving the :math:`x` values that define the edges of the segments and the functional forms :math:`f` of each segment; the normalizations can be computed from the constraint equations. Note that SFH PDFs cannot be described using basic mode, because they are not normalized to unity. Specifying a non-constant SFH requires advanced mode.

An example of a basic mode PDF file is as follows::

   ###############################################################
   # This is an IMF definition file for SLUG v2.
   # This file defines the Chabrier (2005) IMF          
   ###############################################################

   # Breakpoints: mass values where the functional form changes
   # The first and last breakpoint will define the minimum and
   # maximum mass
   breakpoints 0.08 1 120

   # Definitions of segments between the breakpoints

   # This segment is a lognormal with a mean of log_10 (0.2 Msun) 
   # and dispersion 0.55; the dispersion is in log base 10, not 
   # log base e
   segment
   type lognormal
   mean 0.2
   disp 0.55

   # This segment is a powerlaw of slope -2.35
   segment
   type powerlaw
   slope -2.35

This example represents a `Chabrier (2005) <http://adsabs.harvard.edu/abs/2005ASSL..327...41C>`_ IMF from :math:`0.08 - 120` :math:`M_\odot`, which is of the functional form

.. math:: \frac{dp}{dm} \propto \left\{\begin{array}{ll} \exp[-\log(m/m_0)^2/(2\sigma^2)] (m/m_b)^{-1} , & m < m_b \\ \exp[-\log(m_b/m_0)^2/(2\sigma^2)] (m/m_b)^{-2.35}, & m \geq m_b \end{array} \right.,

where :math:`m_0 = 0.2` :math:`M_\odot`, :math:`\sigma = 0.55`, and :math:`m_b = 1` :math:`M_\odot`.

Formally, the format of a basic mode file is as follows. Any line beginning with ``#`` is a comment and is ignored. The first non-empty, non-comment line in a basic mode PDF file must be of the form::

   breakpoints x1 x2 x3 ...

where ``x1``, ``x2``, ``x3``, ``...`` are a non-decreasing series of real numbers. These represent the breakpoints that define the edges of the segment, in units of :math:`M_\odot`. In the example given above, the breakpoints are are :math:`0.08`, :math:`1`, and :math:`120`, indicating that the first segment goes from :math:`0.08 - 1` :math:`M_\odot`, and the second from :math:`1 - 120` :math:`M_\odot`.

After the ``breakpoints`` line, there must be a series of entries of the form::

   segment
   type TYPE
   key1 VAL1
   key2 VAL2
   ...

where ``TYPE`` specifies what functional form describes the segment, and ``key1 VAL1``, ``key2 VAL2``, etc. are a series of (key, value) pairs the define the free parameters for that segment. In the example above, the first segment is described as having a ``lognormal`` functional form, and the keywords ``mean`` and ``disp`` specify that the lognormal has a mean of 0.2 :math:`M_\odot` and a dispersion of 0.55 in :math:`\log_{10}`. The second segment is of type ``powerlaw``, and it has a slope of :math:`-2.35`. The full list of allowed segment types and the keywords that must be specified with them are listed in the :ref:`tab-segtypes` Table. Keywords and segment types are case-insensitive. Where more than one keyword is required, the order is arbitrary.

The total number of segments must be equal to one less than the number of breakpoints, so that each segment is described. Note that it is not necessary to specify a normalization for each segment, as the segments will be normalized relative to one another automatically so as to guarantee that the overall function is continuous.

.. _tab-segtypes:

.. table:: Segment Types

   +-----------------+--------------------------------------------------------+-----------+---------------------------+-----------+-------------------------------------------------+
   | Name            | Functional form                                        | Keyword   | Meaning                   | Keyword   | Meaning                                         |
   +=================+========================================================+===========+===========================+===========+=================================================+
   | ``delta``       | :math:`\delta(x-x_a)`                                  |           |                           |           |                                                 |
   +-----------------+--------------------------------------------------------+-----------+---------------------------+-----------+-------------------------------------------------+
   | ``exponential`` | :math:`\exp(-x/x_*)`                                   | ``scale`` | Scale length, :math:`x_*` |           |                                                 |
   +-----------------+--------------------------------------------------------+-----------+---------------------------+-----------+-------------------------------------------------+
   | ``lognormal``   | :math:`x^{-1} \exp\{-[\log_{10}(x/x_0)]^2/2\sigma^2\}` | ``mean``  | Mean, :math:`x_0`         | ``disp``  | Dispersion in :math:`\log_{10}`, :math:`\sigma` |
   +-----------------+--------------------------------------------------------+-----------+---------------------------+-----------+-------------------------------------------------+
   | ``normal``      | :math:`\exp[-(x-x_0)^2/2\sigma^2]`                     | ``mean``  | Mean, :math:`x_0`         | ``disp``  | Dispersion, :math:`\sigma`                      |
   +-----------------+--------------------------------------------------------+-----------+---------------------------+-----------+-------------------------------------------------+
   | ``powerlaw``    | :math:`x^p`                                            | ``slope`` | Slope, :math:`p`          |           |                                                 |
   +-----------------+--------------------------------------------------------+-----------+---------------------------+-----------+-------------------------------------------------+
   | ``schechter``   | :math:`x^p \exp(-x/x_*)`                               | ``slope`` | Slope, :math:`p`          | ``xstar`` | Cutoff, :math:`x_*`                             |
   +-----------------+--------------------------------------------------------+-----------+---------------------------+-----------+-------------------------------------------------+

Advanced Mode
-------------

In advanced mode, one has complete freedom to set all the parameters describing the PDF: the endpoints of each segment :math:`x_{i,a}` and :math:`x_{i,b}`, the normalization of each segment :math:`n_i`, and the functional forms of each segment :math:`f_i`. This can be used to defined PDFs that are non-continuous, or that are overlapping; the latter option can be used to construct segments with nearly arbitrary functional forms, by constructing a Taylor series approximation to the desired functional form and then using a series of overlapping ``powerlaw`` segments to implement that series.

An example of an advanced mode PDF file is as follows::

   ###############################################################
   # This is a SFH definition file for SLUG v2.
   # This defines a SF history consisting of a series of
   # exponentially-decaying bursts with a period of 100 Myr and
   # a decay timescale of 10 Myr, with an amplitude chosen to
   # give a mean SFR of 10^-3 Msun/yr.
   ###############################################################

   # Declare that this is an advanced mode file
   advanced

   # First exponential burst
   segment
   type exponential
   min      0.0
   max      1.0e8         # Go to 100 Myr
   weight   1.0e5         # Form 10^5 Msun of stars over 100 Myr
   scale    1.0e7         # Decay time 10 Myr

   # Next 4 bursts
   segment
   type exponential
   min      1.0e8
   max      2.0e8
   weight   1.0e5
   scale    1.0e7

   segment
   type exponential
   min      2.0e8
   max      3.0e8
   weight   1.0e5
   scale    1.0e7

   segment
   type exponential
   min      3.0e8
   max      4.0e8
   weight   1.0e5
   scale    1.0e7

   segment
   type exponential
   min      4.0e8
   max      5.0e8
   weight   1.0e5
   scale    1.0e7

This represents a star formation history that is a series of exponential bursts, separated by 100 Myr, with decay times of 10 Myr. Formally, this SFH follows the functional form

.. math:: \dot{M}_* = n e^{-(t\,\mathrm{mod}\, P)/t_{\rm dec}},

where :math:`P = 100` Myr is the period and :math:`t_{\rm dec} = 10` Myr is the decay time, from times :math:`0-500` Myr. The normalization constant :math:`n` is set by the condition that :math:`(1/P) \int_0^P \dot{M}_* \,dt = 0.001` :math:`M_\odot\;\mathrm{yr}^{-1}`, i.e., that the mean SFR averaged over a single burst period is 0.001 :math:`M_\odot\;\mathrm{yr}^{-1}`.

Formally, the format of an advanced mode file is as follows. First, all advanced mode files must start with the line::

   advanced

to declare that the file is in advanced mode. After that, there must be a series of entries of the form::

   segment
   type TYPE
   min MIN
   max MAX
   weight WEIGHT
   key1 VAL1
   key2 VAL2
   ...

The ``type`` keyword is exactly the same as in basic mode, as are the segment-specific parameter keywords ``key1``, ``key2``, :math:`\ldots`. The same functional forms, listed in the :ref:`tab-segtypes` Table, are available as in basic mode. The additional keywords that must be supplied in advanced mode are ``min``, ``max``, and ``weight``. The ``min`` and ``max`` keywords give the upper and lower limits :math:`x_{i,a}` and :math:`x_{i,b}` for the segment; the probability is zero outside these limits. The keyword ``weight`` specifies the integral under the segment, i.e., the weight :math:`w_i` given for segment :math:`i` is used to set the normalization :math:`n_i` via the equation

.. math:: w_i = n_i \int_{x_{i,a}}^{x_{i,b}} f_i(x) \, dx.

In the case of a star formation history, as in the example above, the weight :math:`w_i` of a segment is simply the total mass of stars formed in that segment. In the example given above, the first segment declaration sets up a PDF that with a minimum at 0 Myr, a maximum at 100 Myr, following an exponential functional form with a decay time of :math:`10^7` yr. During this time, a total mass of :math:`10^5` :math:`M_\odot` of stars is formed.

Note that, for the IMF, CMF, and CLF, the absolute values of the weights to not matter, only their relative values. On the other hand, for the SFH, the absolute weight does matter.

.. _sampling_metod_label:

Sampling Methods
----------------

A final option allowed in both basic and advanced mode is a specification of the sampling method. The sampling method is a description of how to draw a population of objects from the PDF, when the population is specified as having a total sum :math:`M_{\rm target}` (usually but not necessarily a total mass) rather than a total number of members :math:`N`; there are a number of ways to do this, which do not necessarily yield identical distributions, even for the same underlying PDF. To specify a sampling method, simply add the line::

   method METHOD

to the PDF file. This line can appear anywhere except inside a ``segment`` specification, or before the ``breakpoints`` or ``advanced`` line that begins the file. The following values are allowed for ``METHOD`` (case-insensitive, as always):

* ``stop_nearest``: this is the default option: draw until the total mass of the population exceeds :math:`M_{\rm target}`. Either keep or exclude the final star drawn depending on which choice brings the total mass closer to the target value.
* ``stop_before``: same as ``stop_nearest``, but the final object drawn is always excluded.
* ``stop_after``: same as ``stop_nearest``, but the final object drawn is always kept.
* ``stop_50``: same as ``stop_nearest``, but keep or exclude the final object with 50% probability regardless of which choice gets closer to the target.
* ``number``: draw exactly :math:`N = M_{\rm target}/\langle M\rangle` object, where :math:`\langle M\rangle` is the expectation value for a single draw.
* ``poisson``: draw exactly :math:`N` objects, where the value of :math:`N` is chosen from a Poisson distribution with expectation value :math:`\langle N \rangle = M_{\rm target}/\langle M\rangle`
* ``sorted_sampling``: this method was introduced by `Weidner & Kroupa (2006, MNRAS. 365, 1333) <http://adsabs.harvard.edu/abs/2006MNRAS.365.1333W>`_, and proceeds in steps. One first draws exactly :math:`N= M_{\rm target}/\langle M\rangle` as in the ``number`` method. If the resulting total mass :math:`M_{\rm pop}` is less than :math:`M_{\rm target}`, the procedure is repeated recursively using a target mass :math:`M_{\rm target} - M_{\rm pop}` until :math:`M_{\rm pop} > M_{\rm target}`. Finally, one sorts the resulting stellar list from least to most massive, and then keeps or removes the final, most massive star using a ``stop_nearest`` policy. 

See the file ``lib/imf/wk06.imf`` for an example of a PDF file with a ``method`` specification.
