.. highlight:: rest

.. _sec-cloudy-slug:

cloudy_slug: An Automated Interface to cloudy
=============================================

SLUG stochastically generates stellar spectra, and it includes an
approximate computation of the nebular lines produced when those
photons interact with the interstellar medium. However, this
approximation ignores a number of potentially important effects, and
does not properly account for the stochastic nature of the stellar
spectra. To perform a much more accurate calculation, SLUG includes an
automated interface to `cloudy <http://nublado.org/>`_ (`Ferland et
al., 2013, RMxAA, 49, 137
<http://adsabs.harvard.edu/abs/2013RMxAA..49..137F>`_). This can be
used to post-process the output of a SLUG run in order to compute
nebular emission.

cloudy_slug Basics
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
   ``cloudy_slug/LineList_HII.dat``. These are the template input files
   that will be used for all the cloudy runs, and their syntax follows
   the standard cloudy syntax. They control things like the density and
   element abundances in the nebula -- see :ref:`ssec-cloudy-template`
   for more details.

4. Perform the desired SLUG simulation. The SLUG simulation outputs
   must include spectra and photometry, and one of the photometric
   bands output must be ``QH0`` (see
   :ref:`ssec-phot-keywords`). Depending on whether one is running in
   integrated or cluster mode (see
   :ref:`sssec-cloudy-integrated-cluster`), either integrated specta
   and photometry or cluster spectra and photometry are
   required.

5. Invoke the cloudy_slug interface script via::

     python cloudy_slug/cloudy_slug.py SLUG_MODEL_NAME

   where ``SLUG_MODEL_NAME`` is the name of the SLUG run to be
   processed. See :ref:`ssec-cloudy-model` for more information on
   the underlying physical model assumed in the calculation, and
   :ref:`ssec-cloudy-slug-options` for more details on the python
   script and its options.

6. The output will be stored as a series of additional output files of
   with names of the form SLUG_MODEL_NAME_*cloudy*.ext, where the
   extension is .txt, .bin, or .fits, depending on the format in which
   the orignal SLUG output was stored. These files can be processed
   automatically by the slugpy helper routines (see
   :ref:`sec-slugpy`). See :ref:`ssec-cloudy-output` for a description
   of the outputs.

Note that some care is required in selecting the conditions passed to
cloudy to ensure that the results are physically sensible. Users are
strongly encouraged to read :ref:`ssec-cloudy-model` to understand
exactly what physical assumptions are being made, and to ensure that
they are reasonable.

   
.. _ssec-cloudy-model:

The cloudy_slug Physical Model
------------------------------

The cloudy_slug code computes emission from a spherical HII region
surrounding a stellar population. The stellar population comes from
SLUG, and the emission calculation is performed with cloudy. Combining
the two requires some physical assumptions and inputs, which are
explained in this section.

.. _sssec-cloudy-integrated-cluster:

Integrated versus Cluster Spectra
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

SLUG outputs both integrated spectra for all the stars in a galaxy,
and spectra for individual clusters. Both the integrated spectra and
the individual cluster spectra can be processed by cloudy. However, it
is important to understand the implicit physical assumptions that one
is making while doing so. If one has a galaxy where all stars are in
clusters (i.e., cluster formation fraction is unity and there is no
cluster disruption), then the integrated starlight spectrum is just
the sum of the individual cluster spectra. For nebular emission,
however, this is not the case: nebular emission does not, in general,
add linearly.

For this reason, if one processes the integrated spectrum through
cloudy, the implicit physical assumption is that the entire galaxy is
a single giant HII region being ionized by the starlight of all the
clusters present. If one processes the individual cluster spectra
instead, the implicit physical picture is that there is no overlap
whatsoever between the HII regions surrounding different star
clusters. Reality almost certainly lies somewhere between these two
extremes, but it is important to understand physically what assumption
one is making by adopting one or the other. We refer to processing the
integrated spectrum as integrated mode, and to processing the
individual cluster spectra as cluster mode. Note that cluster mode can
be very computationally intensive if there are many clusters present,
and that in cluster mode there is no processing of nebular emission
produced by field stars.

In either mode, the spectrum that is used to compute the nebular
emission will be the *unextincted, non-redshifted* spectrum computed
by SLUG.

.. _sssec-cloudy-nebular properties:

Nebular Properties
^^^^^^^^^^^^^^^^^^

Computing the nebular emission requires specifying the physical
properties of the interstellar gas into which the
photons propagate. Codes like cloudy require that the HII region be
described by an inner radius :math:`r_0` and a number density
:math:`n_{0}` of hydrogen nuclei at that radius. One option for
cloudy_slug is that these parameters can be set in the cloudy inputs
as they would be for a normal cloudy run. However, these parameters
are not necessarily the most convenient or descriptive ones with which
to characterize HII regions. For this reason, cloudy_slug allows users
to specify HII region properties in a number of other more convenient
ways.

The basic assumptions made in cloudy_slug's parameterization is
that the HII region is isobaric and isothermal, at all points hydrogen
is fully ionized and helium is singly ionized, and that radiation
pressure is negligible. (Important note: these are the assumptions
used in cloudy_slug's way of writing out the parameters, and they are
approximately true for most HII regions. However, they are *not*
exactly true for the final cloudy calculation, where in general the
temperature is not constant, the ionization states of hydrogen and
helium vary through the nebula, and radiation pressure may or may not
be important.) The HII region occupies a spherical shell bounded by an
inner radius :math:`r_0` and an outer radius :math:`r_1`. The inner
radius is set by the presence of a bubble of shocked stellar wind
material at a temperature :math:`\sim 10^6` K, which is assumed to be
optically thin to ionizing photons. The outer radius is set by the
location where all the ionizing photons have been absorbed.

Under these assumptions, the inner density :math:`n_0` is simply the
(uniform) density :math:`n_{\mathrm{II}}` throughout the ionized
region, and the ionizing photon luminosity passing through a shell of
material at a distance :math:`r` from the stars is 

.. math:: Q(r) = Q(\mathrm{H}^0)
	  \left[1 - \left(\frac{r}{r_S}\right)^3 +
	  \left(\frac{r_0}{r_S}\right)^3\right],

where :math:`Q(\mathrm{H}^0)` is the hydrogen-ionizing luminosity of
the source and :math:`r_S` is the Stromgren radius, given by

.. math:: r_S = \left(\frac{3 Q(\mathrm{H}^0)}{4\pi
	  \alpha_B f_e n_{\mathrm{II}}^2}\right)^{1/3}.

Here :math:`\alpha_B` is the case B recombination coefficient and
:math:`f_e` is the abundance of electrons per H nucleus. For the
purposes of cloudy_slug, we take these two quantities to have
the fixed values :math:`\alpha_B = 2.59\times
10^{-13}\;\mathrm{cm}^3\;\mathrm{s}^{-1}`, appropriate for a
temperature of :math:`10^4` K, and :math:`f_e = 1.1`, appropriate for
a region where He is singly ionized.

From this setup one can define some useful dimensionless numbers. One
is the wind parameter :math:`\Omega` introduced by `Yeh &
Matnzer (2012, ApJ, 757, 108)
<http://adsabs.harvard.edu/abs/2012ApJ...757..108Y>`_, which under the
simple assumptions made in cloudy_slug is given by

.. math:: \Omega = \frac{r_0^3}{r_1^3-r_0^3}

i.e., it is just the ratio of the volume occupied by the wind gas to
that occupied by the photoionized gas. The value of :math:`\Omega`
determines whether winds are important (:math:`\Omega \gg 1`) or
unimportant (:math:`\Omega \ll 1`) for the dynamics of the HII
region. The second dimensionless parameter is the volume-averaged
ionization parameter

.. math:: \mathcal{U} = \frac{3}{4\pi (r_1^3-r_0^3)} \int_{r_0}^{r_1}
	  \left(\frac{Q(r)}{4\pi r^2 c f_i n_{\mathrm{II}}}\right)
	  4\pi r^2 \, dr.

Here :math:`f_i` is the number of free ions per H nucleus, and is
equal to :math:`f_i = 1.1` under the assumption that He is singly
ionized. The quantity in parentheses is the ratio of the ionizing
photon to ion number densities at radius :math:`r`. The value of
:math:`\mathcal{U}` is, together with :math:`n_{\mathrm{II}}`, the
most important factor in determining the output spectrum. A third
useful dimensionless parameter is the ionization parameter at the
inner radius,

.. math:: \mathcal{U}_0 = \frac{Q(\mathrm{H}^0)}
	  {4\pi r_0^2 f_i n_{\mathrm{II}} c}.

The various quantities are not unrelated. It is straightforward to
show that they are constrained by the following relationships:

.. math:: r_0 & = \Omega^{1/3} r_S \\

	  r_1 & = \left(1 + \Omega\right)^{1/3} r_S \\

	  \mathcal{U} & = \left[\frac{81 \alpha_B^2 n_{\mathrm{II}}
	  Q(\mathrm{H}^0)}{256 \pi c^3 f_e}\right]^{1/3}
	  \left[\left(1 + \Omega\right)^{4/3} 
	  - \Omega^{1/3} \left(\frac{4}{3}+\Omega\right)\right] \\

	  & = \left[\frac{81 \alpha_B Q(\mathrm{H}^0)}
	  {64 \pi c^2 f_e r_S}\right]^{1/2}
	  \left[\left(1 + \Omega\right)^{4/3} 
	  - \Omega^{1/3} \left(\frac{4}{3}+\Omega\right)\right] \\

	  \mathcal{U}_0 &= \left[
	  \frac{\alpha_B^2 n_{\mathrm{II}} Q(\mathrm{H}^0)}
	  {36 \pi c^3 f_e}\right]^{1/3} \frac{1}{\Omega^{2/3}} \\

	  &= \frac{4}{9}\Omega^{-2/3} \left[(1+\Omega)^{4/3} -
	  \Omega^{1/3}\left(\frac{4}{3}+\Omega\right)\right]^{-1}
	  \mathcal{U} \\

These relations may be used to compute any four of the quantities
:math:`n_{\mathrm{II}}`, :math:`r_0`, :math:`r_1`, :math:`\mathcal{U}`,
:math:`\mathcal{U}_0` and :math:`\Omega` given the other two.
:ref:`sec-slugpy` provides a class ``hiiregparam`` that can be used
to perform such a computation.

Given this background, cloudy_slug allows the user to specify the
physical properties of the HII region by setting any two of the
following six quantities:

#. The photoionized gas density :math:`n_{\mathrm{II}}`.
#. The inner radius :math:`r_0`.
#. The outer radius :math:`r_1`.
#. The volume-averaged ionization parameter :math:`\mathcal{U}`.
#. The inner radius ionization parameter :math:`\mathcal{U}_0`.
#. The wind parameter :math:`\Omega`.

The two quantities chosen can be specified exactly, or can be drawn
from a specified PDF. One final option, which is only available in
cluster mode, is to obtain the required quantities from a dynamical
model -- see :ref:`sssec-cloudy-dynamical-cluster-mode`.

A few caveats are in order at this point.

#. Not all combinations of values are realizable. In addition to the
   obvious constraints (e.g., :math:`r_1 > r_0`), there are some
   subtle ones. For example, for any given ionizing luminosity
   :math:`Q(\mathrm{H}^0)` and density :math:`n_{\mathrm{II}}`, the
   value of :math:`\mathcal{U}` is bounded from above. Increasing the
   wind parameter :math:`\Omega` can allow arbitrarily small values of
   :math:`\mathcal{U}`, but not arbitrarily large ones. If the user
   requests a physically impossible combination of parameters,
   cloudy_slug will note the problem and react as specified by the
   options given to the cloudy_slug script.
#. Even for parameters that are not physically impossible, the results
   may not be sensible, and may cause cloudy to crash in extreme
   cases. For example, if one sets :math:`\Omega = 0` and
   :math:`\mathcal{U} = 10^{-4}`, then for an ionizing lumnosity of
   :math:`Q(\mathrm{H}^0) = 10^{50}` photons/s (typical for a cluster
   of :math:`\sim 10^4M_\odot`), the corresponding density is
   :math:`n_{\mathrm{II}} \approx 10^{-5}\mbox{ cm}^{-3}`! As this
   density the gas will be fully ionized by cosmic rays and the
   extragalactic background, and it makes no sense to think of it as
   an HII region. Caution is required.
#. The parameter combinations :math:`(r_0,\mathcal{U})` and
   :math:`(r_1,\mathcal{U}_0)` are not allowed
   because they do not define a unique solution for the other
   parameters (the resulting equations have multiple physically-valid
   solutions).
#. The relations given above are only valid if radiation pressure is
   not dynamically significant. If it is, then there are no known
   analytic relations between the various quantities. The cloudy_slug
   code will still run, and will use the relations above, but the
   actual HII region properties may be markedly different from those
   requested. In cases where radiation pressure is important, it is
   generally advisable to save the HII region physical conditions
   output by cloudy to compute quatities from them directly. The
   cloudy_slug script will issue a warning if radiation pressure is
   expected to be significant for the HII region being computed. As a rule
   of thumb, radiation pressure is significant if

.. math:: \zeta \equiv \frac{r_{\mathrm{ch}}}{r_1} > 1

where :math:`r_{\mathrm{ch}}` is the
characteristic radius defined by `Krumholz & Matzner (2009, ApJ,
703, 1352) <http://adsabs.harvard.edu/abs/2009ApJ...703.1352K>`_ as

.. math:: r_{\mathrm{ch}} & = 
	  \frac{\alpha_B}{12 \pi \phi}
	  \left(\frac{\epsilon_0}{2 f_e k_B T_{\mathrm{II}}}\right)^2
	  f_{\mathrm{trap}}^2 \frac{\psi^2 Q(\mathrm{H}^0)}{c^2}
	  \\

Here  :math:`\phi = 0.73` is the fraction of ionizing photons absorbed
by hydrogen atoms rather than dust, :math:`\epsilon_0 =
13.6\;\mathrm{eV}` is the hydrogen ionization potential,
:math:`T_{\mathrm{II}} = 10^4\;\mathrm{K}` is the temperature inside
the HII region, :math:`f_{\mathrm{trap}} = 2` is the trapping factor
that accounts for stellar wind and trapped infrared radiation
pressure, and :math:`\psi = 3.2` is the mean photon energy in Rydberg for
a fully sampled IMF at zero age. 


.. _sssec-cloudy-dynamical-cluster-mode:

Dynamical Mode
^^^^^^^^^^^^^^

In cluster mode, cloudy_slug allows an additional option to derive the
physical properties of the HII region. They can be computed from a
dynamical model of HII region expansion, taken from `Krumholz &
Matzner (2009, ApJ, 703, 1352)
<http://adsabs.harvard.edu/abs/2009ApJ...703.1352K>`_.  In this model,
the radius of an HII region can be computed as a function of the
ionizing luminosity :math:`Q(\mathrm{H}^0)`, ambient hydrogen number
density :math:`n_{\mathrm{H}}`, and star cluster age :math:`t` as

.. math::

   r_1 & = r_{\mathrm{ch}}
   \left(x_{\mathrm{rad}}^{7/2} +
   x_{\mathrm{gas}}^{7/2}\right)^{2/7} \\

   x_{\mathrm{rad}} &= (2\tau^2)^{1/4} \\

   x_{\mathrm{gas}} &= (49\tau^2/36)^{2/7} \\

   \tau &= t/t_{\mathrm{ch}} \\

   t_{\mathrm{ch}} & = \left(\frac{4\pi \mu m_{\mathrm{H}}
   n_{\mathrm{H}} c r_{\mathrm{ch}}^4}{3 f_{\mathrm{trap}}
   Q(\mathrm{H}^0) \psi \epsilon_0}\right)^{1/2} \\

Definitions of various quantities appearing in these equations are
given above. The quantity :math:`\mu = 1.4` is the mean
mass per hydrogen nucleus for gas of the standard cosmic
composition.

We refer to this method of computing HII region properties as
dynamical mode. In this mode, a user can specify the properties of the
nebula in terms of an ambient density :math:`n_{\mathrm{H}}` and
a wind parameter :math:`\Omega`. All other quantities are
derived from these two and from the ionizing luminosity
:math:`Q(\mathrm{H}^0)` and age :math:`t` of each cluster. Dynamical
mode can only be used in combination with cluster mode, not integrated
mode, because composite stellar populations do not have well-defined
ages.


.. _ssec-cloudy-template:

The cloudy_slug Input Template
------------------------------

The cloudy_slug interface operates by reading SLUG output spectra and
using them as inputs to a cloudy calculation. However, cloudy
obviously requires many input parameters beyond simply the spectrum of
the input radiation field. These parameters are normally provided by
an input file whose format is as described in the `cloudy documentation
<http://nublado.org>`_. The cloudy_slug interface works by reading a
*template* input file that specifies all these parameter, and which
will be used as a basis for the final cloudy input files that will
contain the SLUG spectra.

In general the template input file looks just like an ordinary cloudy
input file, subject to the following restrictions:

#. The input file *must not* contain any commands that specify the
   luminosity, intensity, or the spectral shape. These will be
   inserted automatically by the cloudy_slug script.
#. The input file *may* contain a radius command specifying the inner
   radius of the HII region. If it does not, then the user must
   specify the radius in another way, by setting 2 of the 6 inputs
   described in :ref:`sssec-cloudy-nebular properties` (for
   simulations not done in dynamic mode) or by setting an ambient
   density and wind parameters in
   :ref:`sssec-cloudy-dynamical-cluster-mode`. If the user does
   set these quantities, any radius command in the template file will be
   ignored, and a warning message will be issued if one is
   found. Finally, note that cloudy_slug will only compute derived
   parameters correctly from a radius in the template file if the
   radius is specified in cloudy's default format, by giving a log of
   the radius in cm; the keywords "linear" and "parsecs" are not
   currently supported.
#. The input file *may* contain a hydrogen density command specifying
   the starting hydrogen density. The rules for this are the same as
   for the radius command.
#. Any outputs to be written (specified using the ``save`` or
   ``punch`` keywords) must give file names containing the string
   ``$OUTPUT_FILENAME``. This string will be replaced by the
   cloudy_slug script to generate a unique file name for each cloudy
   run, and to read back these outputs for post-processing.
#. The cloudy_slug output will contain output spectra only if the
   cloudy input file contains a ``save last continuum`` command. See
   :ref:`ssec-cloudy-output`.
#. The cloudy_slug output will contain output line luminosities only
   if the cloudy input file contains a ``save last line list emergent
   absolute column`` command. See :ref:`ssec-cloudy-output`.
#. The cloudy_slug output will contain output physical conditions and
   dimensionless values only if the cloudy input file contains a
   ``save last hydrogen conditions`` command. See
   :ref:`ssec-cloudy-output`.
#. If any other outputs are produced by the input file, they will
   neither be processed nor moved, deleted, or otherwise changed by
   the cloudy_slug script.
#. Running cloudy in grid mode is not currently supported.

An example cloudy input file with reasonable parameter choices is
provided as ``cloudy_slug/cloudy_in.template`` in the main directory
of the SLUG repository.

In addition to the input file, the default template makes use of a
cloudy line list file to specify which line luminosities should be
output (see the `cloudy documentation <http://nublado.org>`_ for
details). The template points to the file
``cloudy_slug/LineList_HII.data`` (which is identical to cloudy's
default line list for HII regions), but any other valid cloudy line
list file would work as well.

.. _ssec-cloudy-slug-options:

The cloudy_slug Interface Script
--------------------------------

The ``cloudy_slug.py`` script provides the interface between SLUG and
cloudy. Usage for this script is as follows::

  cloudy_slug.py [-h] [-a AGEMAX] [--cloudypath CLOUDYPATH]
                 [--cloudytemplate CLOUDYTEMPLATE] [-cm]
                 [-cf COVERINGFAC] [-d] [-hd HDEN] [-ip IONPARAM]
                 [-ip0 IONPARAM0] [-ipm IONPARAMMAX]
		 [--ionparammin IONPARAMMIN] [-nl NICELEVEL]
                 [-n NPROC] [-ps PARAMSAFETY] [-qm QH0MIN] [-r0 R0]
                 [-r1 R1] [-s] [--slugformat SLUGFORMAT]
                 [--slugpath SLUGPATH] [-t TMPDIR] [-v] [-wp WINDPARAM]
                 [-wr]
                 slug_model_name [start_spec] [end_spec]

The positional arguments are as follows:

* ``slug_model_name``: this is the name of the SLUG output to be used
  as a basis for the cloudy calculation. This should be the same as
  the ``model_name`` parameter used in the SLUG simulation, with the
  optional addition of a path specification in front.
* ``start_spec``: default behavior is to run cloudy on all the
  integrated spectra or cluster spectra (see
  :ref:`sssec-cloudy-integrated-cluster`). If this
  argument is set, cloudy will only be run in spectra starting with
  the specified trial number or cluster number; numbers are
  0-offset, so the first trial/cluster is 0, the next is 1, etc.
* ``end_spec``: same as ``start_spec``, but specifying the last
  cluster to be processed. Per standard python convention, the spectra
  processed will go up to but not include ``end_spec``.

The following optional arguments control paths and file locations:

* ``--cloudypath CLOUDYPATH``: path to the cloudy executable; default
  is ``$CLOUDY_DIR/cloudy.exe``
* ``--cloudytemplate CLOUDYTEMPLATE``: cloudy input file template (see
  :ref:`ssec-cloudy-template`); default is
  ``$SLUG_DIR/cloudy_slug/cloudy.in_template``
* ``--slugformat SLUGFORMAT``: the format of slug output data to use;
  valid values are ``ascii``, ``bin``, ``binary``, and ``fits``. By
  default ``cloudy_slug`` checks for any output whose name and path
  match the model name and search path, regardless of format.
* ``--slugpath SLUGPATH``: path to the SLUG output data. If not set,
  cloudy_slug searches for an appropriately-named set of output files
  first in the current working directory, and next in
  ``$SLUG_DIR/output``
* ``-t TMPDIR, --tmpdir TMPDIR``: location of the temporary directory
  where temporary files should be stored; defaults to
  `./cloudy_tmp_MODEL_NAME`.

The following arguments control how HII regions are processed:
  
* ``-a AGEMAX, --agemax AGEMAX``: maximum cluster age in Myr for
  cloudy computation. Cloudy will not be run on clusters older than
  this value, and the predicted nebular emission for such clusters
  will be recorded as zero. Default value is 10 Myr. This argument only
  has an effect if running in cluster mode (see
  :ref:`sssec-cloudy-integrated-cluster`); otherwise it is ignored.
* ``-cf COVERINGFRAC, --coveringfrac COVERINGFRAC``: this sets the
  covering fraction of the HII region, i.e., the fraction of ionizing
  photons that are assumed to produce nebular emission; the output
  luminosity is decreased by a factor of the covering fraction
* ``-cm, --clustermode``: if this argument is set, then cloudy_slug
  will process cluster spectra; the default behavior is to process
  integrated spectra
* ``--ionparammax IONPARAMMAX``: maximum value for the inner radius
  ionization parameter :math:`\mathcal{U}_0`. If the value falls
  outside this range, the behavior is controlled by the setting of the
  ``paramsafety`` option (see below).
* ``--ionparammmin IONPARAMMIN``: same as ``ionparammax``, but sets a
  minimum instead of a maximum. 
* ``-qm QH0MIN, --qH0min QH0MIN``: minimum ionizing luminosity for
  which to run cloudy (default = 0). As with ``--agemax``, for
  clusters / times where cloudy is not run, that case will still
  appear in the output, but the nebular spectra and nebular line
  luminosities will all be set to zero.

The following parameters specify the physical properties of HII
regions, as explained in :ref:`sssec-cloudy-nebular
properties`. Parameters can be set to either fixed values, or to the
names of PDF files. Any numerical value given is interpreted as a
fixed constant, while non-numerical values are interpreted as the
names of PDF files that specify a PDF from which the corresponding
parameter is to be drawn. See :ref:`sec-pdfs` for details on PDF file
formats.

* ``-hd HDEN, --hden HDEN``: hydrogen density in HII region,
  :math:`n_{\mathrm{II}}`
* ``-ip IONPARAM, --ionparam IONPARAM``: volume-averaged ionization
  parameter, :math:`\mathcal{U}`
* ``-ip0 IONPARAM, --ionparam0 IONPARAM0``: ionization parameter at
  HII region inner edge, :math:`\mathcal{U}_0`
* ``-r0 R0``: inner radius of the HII region
* ``-r1 R1``: outer radius of the HII region
* ``-wp WINDPARAM, --windparam WINDPARAM``: wind parameter,
  :math:`\Omega`

The following arguments control general code behavior:

* ``-h, --help``: prints a help message and then exits
* ``-nl NICELEVEL, --nicelevel NICELEVEL``: if this is set, then the
  cloudy processes launched by the script will be run at this nice
  level. If it is not set, they will not be nice'd. Note that this
  option will only work correctly on platforms that support nice.
* ``-n NPROC, --nproc NPROC``: number of simultaneous cloudy processes
  to run; default is the number of cores available on the system
* ``--ps PARAMSAFETY, --paramsafety PARAMSAFETY``: specifies how to
  handle situations where the combination of input HII region
  parameters is not physically allowed, or falls outside the bounds
  set by ``ionparammin`` and ``ionparammax``. Available options are:
  
  *  ``warn``: one of the input parameters is adjusted to bring it
     to a physically-allowed value, and a warning is issued; the run
     continues. This is the default behaviour.
  *  ``skip``: runs with unphysical parameter choices are skipped; the
     parameters that were chosen are recorded in the output file, but
     the output spectrum, all line luminosities, and all other
     parameters are set to 0, and cloudy is not run.
  *  ``halt``: if a forbidden parameter combination is found, the
     entire cloudy_slug run is halted.
  *  ``redraw``: if a forbidden parameter combination is found, and
     one or more parameters are being drawn from a PDF, a new set of
     parameters will be drawn from the PDF. Redrawing will continue up
     to 100 times until a physically-allowed parameter combination is
     found. If no valid parameter combination is found after 100
     attempts, revert to ``skip``.
* ``-s, --save``: by default, cloudy_slug will extract line and
  spectral data from the cloudy outputs and store them as described in
  :ref:`ssec-cloudy-output`, then delete the cloudy output files. If
  this option is set, the cloudy output files will NOT be deleted, and
  will be left in place, in sub-directory of the working directory
  called ``cloudy_tmp_MODEL_NAME`` where ``MODEL_NAME`` is the SLUG model
  name. WARNING: cloudy's outputs are written in ASCII and are quite
  voluminous, so choose this option only if you are running
  cloudy on a small number of SLUG spectra and/or you are prepared to
  store hundreds of GB or more. The data that ``cloudy_slug`` extract
  are much, much smaller, and (if you do not use ASCII format) are
  stored in a much more compact form.
* ``-v, --verbose``: if this option is set, cloudy_slug produces
  verbose output as it runs
* ``-wr, --writeparams``: if set, this option causes ``cloudy_slug``
  to write out a file beginning with ``cloudy_slug.param`` for each
  cloudy run. This file is written in the same directory used by the
  save command, and it contains an ASCII printout of the various
  parameters. This option is only applied if ``--save`` is also set.

.. _ssec-cloudy-output:

Full Description of cloudy_slug Output
--------------------------------------

The cloudy_slug script will automatically process the cloudy output
and produce a series of new output files, which will be written to the
same directory where the input SLUG files are located, and using the
same output mode (ASCII text, raw binary, or FITS -- see
:ref:`sec-output`). If cloudy_slug is called to process integrated
spectra, the four output files will be
``MODEL_NAME_integrated_cloudyparams.ext``,
``MODEL_NAME_integrated_cloudylines.ext``, 
``MODEL_NAME_integrated_cloudyphot.ext``, and 
``MODEL_NAME_integrated_cloudyspec.ext``, where the extension ``.ext``
is one of ``.txt``, ``.bin``, or ``.fits``, depending on the
``output_mode``. If cloudy_slug is run on cluster spectra, the four
output files will be
``MODEL_NAME_cluster_cloudyparams.ext``,
``MODEL_NAME_cluster_cloudylines.ext``, 
``MODEL_NAME_cluster_cloudyphot.ext``, and 
``MODEL_NAME_cluster_cloudyspec.ext``. All of these output files will
be read and processed automatically if the outputs are read using
``read_integrated`` or ``read_cluster`` in the :ref:`sec-slugpy`
library.

The format of these files is described below.

The ``integrated_cloudyparams`` File
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This file contains the input parameters for the cloudy runs, and
quantities derived from them. All parameters are as defined in
:ref:`sssec-cloudy-nebular properties`. The output file consists of a
series of entries containin the following fields:

* ``Trial``: which trial these data are from
* ``Time``: evolution time at which the output is produced
* ``Hden``: number density of hydrogen nuclei at the inner edge of the
  HII region, in H/cm^3
* ``R0``: radius of the inner edge of the HII region, in cm
* ``R1``: radius of the outer edge of the HII region, in cm
* ``QH0``: ionizing luminosity input to cloudy, in photons/s
* ``CovFac``: covering factor used
* ``U``: volume-averaged ionization parameter :math:`\mathcal{U}`
* ``U0``: inner edge ionization parameter :math:`\mathcal{U}_0`
* ``Omega``: wind parameter :math:`\Omega`
* ``zeta``: radiation pressure parameter :math:`\zeta`

It will also contain the following fields if the cloudy template file
includes a ``save last hydrogen conditions`` command:

* ``Hden_out``: mean H number density for the HII region structure
  computed by cloudy, in H/cm^3; the average is weighted by the
  ionized volume, i.e., it is weighted by :math:`x\, dV`, where
  :math:`x` is the hydrogen ionization fraction.
* ``R1_out``: HII region outer radius returned by cloudy
* ``Omega_out``: wind parameter :math:`\Omega`, computed using
  ``R1_out`` instead of ``R1``
* ``zeta_out``: radiation pressure parameter :math:`\zeta`, computing
  using ``R1_out`` instead of ``R1``

If the SLUG data input to cloudy_slug were written in ``ascii`` mode,
these data are output as a text file containing a series of columns,
with different trials separated by lines of dashes.

If the SLUG data input to cloudy_slug were written in ``fits`` mode,
the data are written in a FITS file containing a binary table
extension. The table contains one column whose name corresponds to the
list of fields above.

If the SLUG data input to cloudy_slug were writtin in ``binary`` mode,
these data are written in a raw binary file that is formatted as
follows. First, there are four bytes specifying if the optional fields
are included:

* ``Hden_out_set`` (byte): 0 if the data do not include ``Hden_out``, 1 if
  they do include it
* ``R1_out_set`` (byte): 0 if the data do not include ``R1_out``, 1 if
  they do include it
* ``Omega_out_set`` (byte): 0 if the data do not include ``Omega_out``, 1 if
  they do include it
* ``zeta_out_set`` (byte): 0 if the data do not include ``zeta_out``, 1 if
  they do include it

This is followed by a series of records containing the following fields:

* ``Trial`` (numpy ``uint64``)
* ``Time`` (numpy ``float64``)
* ``Hden`` (numpy ``float64``)
* ``R0`` (numpy ``float64``)
* ``R1`` (numpy ``float64``)
* ``QH0`` (numpy ``float64``)
* ``covFac`` (numpy ``float64``)
* ``U`` (numpy ``float64``)
* ``U0`` (numpy ``float64``)
* ``Omega`` (numpy ``float64``)
* ``zeta`` (numpy ``float64``)
* ``Hden_out`` (numpy ``float64``; optional, only if the relevant byte
  in the header is set to )
* ``R1_out`` (numpy ``float64``; optional, only if the relevant byte
  in the header is set to )
* ``Omega_out`` (numpy ``float64``; optional, only if the relevant byte
  in the header is set to )
* ``zeta_out`` (numpy ``float64``; optional, only if the relevant byte
  in the header is set to 1)

There is one such record for each output time, with different trials
ordered sequentially, so that all the times for one trial are output
before the first time for the next trial.

The ``integrated_cloudylines`` File
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This file contains data on the nebular line emission produced by the
interaction of the stellar radiation field with the ISM. It consists
of a series of entries containing the following fields:

* ``Time``: evolution time at which the output is produced
* ``LineLabel``: four letter code labeling each line. These codes
  are the codes used by cloudy (see the `cloudy documentation
  <http://nublado.org>`_)
* ``Wavelength``: wavelength of the line, in Angstrom. Note that
  default cloudy behavior is to round wavelengths to the nearest
  Angstrom.
* ``Luminosity``: line luminosity, in erg/s

If the SLUG data input to cloudy_slug were written in ``ascii`` mode,
these data are output as a text file containing a series of columns,
with different trials separated by lines of dashes.

If the SLUG data input to cloudy_slug were written in ``fits`` mode,
the data are written in a FITS file containing two binary table
extensions. The first extension contains two fields, ``Line_label`` and
``Wavelength``, giving the four-letter cloudy line codes and central
wavelengths. The second extension contains three columns, giving the
trial number, time, and line luminosity for each line at each time in
each trial.

If the SLUG data input to cloudy_slug were written in ``binary`` mode,
the data are written in a raw binary file. The file starts with a
header consisting of

* ``NLine`` (python ``int``, equivalent to C ``long``): number of lines
* ``LineLabel`` (``NLine`` entries stored as ``ASCII text``): line
  labels listed in ASCII, one label per line

This is followed by a series of entries of the form

* ``Trial`` (numpy ``uint64``)
* ``Time`` (numpy ``float64``)
* ``LineLum`` (``NLine`` entries of type numpy ``float64``)

There is one such record for each output time, with different trials
ordered sequentially, so that all the times for one trial are output
before the first time for the next trial.

.. _sssec-int-cloudyspec-file:

The ``integrated_cloudyspec`` File
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This file contains data on the spectrum produced by interaction
between the stellar radiation field and the nebula. Each entry in the
output file contains the folling fields:

* ``Trial``: which trial these data are from
* ``Time``: evolution time at which the output is produced
* ``Wavelength``: the wavelength at which the spectrum is evaluated,
  in Angstrom
* ``Incident``: specific luminosity in erg/s/Angstrom at the specified
  wavelength. In cloudy's terminology, this is the *incident*
  spectrum, i.e., the stellar radiation field entering the nebula. It
  should be the same as the spectrum contained in the SLUG
  ``integrated_spec`` file for the corresponding time and trial,
  except interpolated onto the wavelength grid used by cloudy.
* ``Transmitted``:  specific luminosity in erg/s/Angstrom at the specified
  wavelength. In cloudy's terminology, this is the *transmitted*
  spectrum, i.e., the stellar spectrum exiting the HII region, not
  including any emission produced within the nebula. This is what
  would be detected by an observing aperture that included only the
  stars, and none of the nebula.
* ``Emitted``:  specific luminosity in erg/s/Angstrom at the specified
  wavelength. In cloudy's terminology, this is the *emitted*
  spectrum, i.e., the spectrum emitted by the diffuse gas in the HII
  region, excluding any light from the stars themselves. This is what
  would be seen by an observer whose aperture covered the nebula, but
  masked the stars.
* ``Transmitted_plus_emitted``: this is just the sum of
  ``Transmitted`` and ``Emitted``. It represents what would be
  observed in an aperture including both the stars and the HII
  region.

If the SLUG data input to cloudy_slug were written in ``ascii`` mode,
these data are output as a text file containing a series of columns,
with different trials separated by lines of dashes.

If the SLUG data input to cloudy_slug were written in ``fits`` mode,
these data are written in a FITS file containing two binary table
extensions. The first extension contains one field, ``Wavelength``,
which gives the wavelengths of the spectra in Angstrom. The second
extension contains six fields: ``Trial``, ``Time``, 
``Incident_spectrum``, ``Transmitted_spectrum``, ``Emitted_spectrum``,
and ``Transmitted_plus_emitted_spectrum``. The first two of these give
the trial number and time, and the remaining four give the incident,
transmitted, emitted, and transmitted plus emitted spectra for the
corresponding time and trial.

If the SLUG data input to cloudy_slug were written in ``binary`` mode,
these data are written in a raw binary file that is formatted as
follows. The file begins with a header consisting of

* ``NWavelength`` (numpy ``int64``): number of wavelengths
* ``Wavelength`` (``NWavelength`` entries of numpy ``float64``)

and then contains a series of records of the form

* ``Trial`` (numpy ``uint64``)
* ``Time`` (numpy ``float64``)
* ``Incident`` (``NWavelength`` entries of numpy ``float64``)
* ``Transmitted`` (``NWavelength`` entries of numpy ``float64``)
* ``Emitted`` (``NWavelength`` entries of numpy ``float64``)
* ``Transmitted_plus_emitted`` (``NWavelength`` entries of numpy
  ``float64``)

There is one such record for each output time, with different trials
ordered sequentially, so that all the times for one trial are output
before the first time for the next trial.

The ``integrated_cloudyphot`` File
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This file contains photometric data computed for the spectra produced
by the interaction between the stellar radiation field and the HII
region. The file consists of a series of entries containing the
following fields:

* ``Trial``: which trial these data are from
* ``Time``: evolution time at which the output is computed
* ``PhotFilter1_trans``: photometric value for the *Transmitted*
  radiation field through filter 1, where filter 1 here is the same as
  filter 1 in :ref:`ssec-int-phot-file`; units are also the same as
  in that file.
* ``PhotFilter1_emit``: photometric value for the *Emitted*
  radiation field through filter 1
* ``PhotFilter1_trans_emit``: photometric value for the
  *Transmitted_plus_emitted* radiation field through filter 1
* ``PhotFilter2_trans``
* ``PhotFilter2_emit``
* ``PhotFilter2_trans_emit``
* ``...``

For distinctions between the *Transmitted*, *Emitted*, and
*Transmitted_plus_emitted* radiation fields, see
:ref:`sssec-int-cloudyspec-file`, or the `cloudy documentaiton
<http://nublado.org>`_. Note that we do not record photometry for the
incident spectrum, since that would be, up to the accuracy of the
numerical integration, identical to the photometry already recorded in
the :ref:`ssec-int-phot-file`.

If the SLUG data input to cloudy_slug were written in ``ascii`` mode,
these data are output as a text file containing a series of columns,
with different trials separated by lines of dashes.

If the SLUG data input to cloudy_slug were written in ``fits`` mode,
these data are written in a FITS file containing one binary table
extension, consisting of a series of columns. The columns are
``Trial``, ``Time``, ``Filter1_Transmitted``, ``Filter1_Emitted``,
``Filter1_Transmitted_plus_emitted``, ``...``. The first two columns
give the trial number and the time, and the remainder give the
photometric values for the transmitted, emitted, and transmitted plus
emitted spectra in each filter.

If the SLUG data input to cloudy_slug were written in ``binary`` mode,
these data are written to a raw binary file that is formatted as
follows. The file starts with an ASCII header consisting of the
following, each on a separate line:

* ``NFilter`` (stored as ``ASCII text``): number of filters used
* ``FilterName`` ``FilterUnit`` (``NFilter`` entries stored as ``ASCII
  text``): the name and units for each filter are listed in ASCII, one
  filter-unit pair per line

This is followed by a series of entries of the form:

* ``PhotFilter_Transmitted`` (``NFilter`` entries of numpy
  ``float64``), giving the transmitted photometry in each filter
* ``PhotFilter_Emitted`` (``NFilter`` entries of numpy
  ``float64``), giving the emitted photometry in each filter
* ``PhotFilter_Transmitted_plus_emitted`` (``NFilter`` entries of numpy
  ``float64``), giving the transmitted plus emitted photometry in each
  filter

There is one such record for each output time, with different trials
ordered sequentially, so that all the times for one trial are output
before the first time for the next trial.

The ``cluster_cloudyparams`` File
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This file contains the input parameters for the cloudy runs, and
quantities derived from them. It consists of a series of entries
containin the following fields:

* ``UniqueID``: a unique identifier number for each cluster that is
  preserved across times and output files
* ``Time``: evolution time at which the output is produced
* ``Hden``: number density of hydrogen nuclei at the inner edge of the
  HII region whose structure cloudy computes, in H/cm^3
* ``R0``: radius of the inner edge of the HII region, in cm
* ``R1``: radius of the outer edge of the HII region, cm
* ``QH0``: ionizing luminosity input to cloudy, in photons/s
* ``CovFac``: covering factor used
* ``U``: volume-averaged ionization parameter :math:`\mathcal{U}`
* ``U0``: inner edge ionization parameter :math:`\mathcal{U}_0`
* ``Omega``: wind parameter :math:`\Omega`
* ``zeta``: radiation pressure parameter :math:`\zeta`

It will also contain the following fields if the cloudy template file
includes a ``save last hydrogen conditions`` command:

* ``Hden_out``: mean H number density for the HII region structure
  computed by cloudy, in H/cm^3; the average is weighted by the
  ionized volume, i.e., it is weighted by :math:`x\, dV`, where
  :math:`x` is the hydrogen ionization fraction.
* ``R1_out``: HII region outer radius returned by cloudy
* ``Omega_out``: wind parameter :math:`\Omega`, computed using
  ``R1_out`` instead of ``R1``
* ``zeta_out``: radiation pressure parameter :math:`\zeta`, computing
  using ``R1_out`` instead of ``R1``

If the SLUG data input to cloudy_slug were written in ``ascii`` mode,
these data are output as a text file containing a series of columns,
with different trials separated by lines of dashes.

If the SLUG data input to cloudy_slug were written in ``fits`` mode,
the data are written in a FITS file containing a binary table
extension. The table contains one column whose name corresponds to the
list of fields above.

If the SLUG data input to cloudy_slug were writtin in ``binary`` mode,
these data are written in a raw binary file that is formatted as
follows. First, there are four bytes specifying if the optional fields
are included:

* ``Hden_out_set`` (byte): 0 if the data do not include ``Hden_out``, 1 if
  they do include it
* ``R1_out_set`` (byte): 0 if the data do not include ``R1_out``, 1 if
  they do include it
* ``Omega_out_set`` (byte): 0 if the data do not include ``Omega_out``, 1 if
  they do include it
* ``zeta_out_set`` (byte): 0 if the data do not include ``zeta_out``, 1 if
  they do include it


Next there are a series of records, one for each output time, with
different trials ordered sequentially, so that all the times for one
trial are output before the first time for the next trial. Each record
consists of a header containing

* ``Time`` (``double``)
* ``NCluster`` (``std::vector<double>::size_type``, usually ``unsigned long long``): number of non-disrupted clusters present at this time

This is followed by ``NCluster`` entries of the following form:

* ``UniqueID`` (numpy ``uint64``)
* ``Time`` (numpy ``float64``)
* ``Hden`` (numpy ``float64``)
* ``R0`` (numpy ``float64``)
* ``R1`` (numpy ``float64``)
* ``QH0`` (numpy ``float64``)
* ``covFac`` (numpy ``float64``)
* ``U`` (numpy ``float64``)
* ``U0`` (numpy ``float64``)
* ``Omega`` (numpy ``float64``)
* ``zeta`` (numpy ``float64``)
* ``Hden_out`` (numpy ``float64``; optional, only if the relevant byte
  in the header is set to 1)
* ``R1_out`` (numpy ``float64``; optional, only if the relevant byte
  in the header is set to 1)
* ``Omega_out`` (numpy ``float64``; optional, only if the relevant byte
  in the header is set to 1)
* ``zeta_out`` (numpy ``float64``; optional, only if the relevant byte
  in the header is set to 1)

The ``cluster_cloudylines`` File
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This file contains data on the nebular line emission produced by the
interaction of the stellar radiation field with the ISM around each
cluster. It consists of a series of entries containing the following
fields:

* ``UniqueID``: a unique identifier number for each cluster that is
  preserved across times and output files
* ``Time``: evolution time at which the output is produced
* ``LineLabel``: four letter code labeling each line. These codes
  are the codes used by cloudy (see the `cloudy documentation
  <http://nublado.org>`_)
* `` Wavelength``: wavelength of the line, in Angstrom. Note that
  default cloudy behavior is to round wavelengths to the nearest
  Angstrom.
* ``Luminosity``: line luminosity, in erg/s

If the SLUG data input to cloudy_slug were written in ``ascii`` mode,
these data are output as a text file containing a series of columns,
with different trials separated by lines of dashes.

If the SLUG data input to cloudy_slug were written in ``fits`` mode,
the data are written in a FITS file containing two binary table
extensions. The first extension contains two fields, ``Line_label`` and
``Wavelength``, giving the four-letter cloudy line codes and central
wavelengths. The second extension contains four columns, giving the
unique ID, trial number, time, and line luminosity for each line at
each time in each trial.

If the SLUG data input to cloudy_slug were written in ``binary`` mode,
the data are written in a raw binary file. The file starts with a
header consisting of

* ``NLine`` (python ``int``, equivalent to C ``long``): number of lines
* ``LineLabel`` (``NLine`` entries stored as ``ASCII text``): line
  labels listed in ASCII, one label per line

This is followed by a series of records, one for each output time,
with different trials ordered sequentially, so that all the times for
one trial are output before the first time for the next trial. Each
record consists of a header containing

* ``Time`` (``double``)
* ``NCluster`` (``std::vector<double>::size_type``, usually ``unsigned long long``): number of non-disrupted clusters present at this time

This is followed by ``NCluster`` entries of the following form:

* ``UniqueID`` (numpy ``uint64``)
* ``LineLum`` (``NLine`` entries of numpy ``float64``)

The ``cluster_cloudyspec`` File
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This file contains data on the spectra produced by the interaction of
the stellar radiation field with the ISM around each cluster. It
consists of a series of entries containing the following fields:

* ``UniqueID``: a unique identifier number for each cluster that is
  preserved across times and output files
* ``Time``: evolution time at which the output is produced
* ``Wavelength``: observed frame wavelength at which the spectrum is evaluated
* ``Incident``: specific luminosity in erg/s/Angstrom at the specified
  wavelength for the *incident* radiation field
* ``Transmitted``: specific luminosity in erg/s/Angstrom at the specified
  wavelength for the *transmitted* radiation field
* ``Emitted``: specific luminosity in erg/s/Angstrom at the specified
  wavelength for the *emitted* radiation field
* ``Transmitted_plus_emitted``: specific luminosity in erg/s/Angstrom
  at the specified wavelength for the *transmitted plus emitted*
  radiation field

For explanations of the distinction between the incident, transmitted,
emitted, and transmitted plus emitted radiation fields, see
:ref:`sssec-int-cloudyspec-file`.

If the SLUG data input to cloudy_slug were written in ``ascii`` mode,
these data are output as a text file containing a series of columns,
with different trials separated by lines of dashes.

If the SLUG data input to cloudy_slug were written in ``fits`` mode,
these data are written in a FITS file containing two binary table
extensions. The first table contains a column ``Wavelength`` listing
the wavelengths at which the spectra are given. The second table
consists of seven columns: ``Trial``, ``UniqueID``, ``Time``,
``Incident_spectrum``, ``Transmitted_spectrum``, ``Emitted_spectrum``,
and ``Transmitted_plus_emitted_spectrum``. The first three of these
give the trial number, unique ID of the cluster, and the time. The
remaining four give the incident, transmitted, emitted, and
transmitted plus emitted spectra for the corresponding cluster.

If the SLUG data input to cloudy_slug were written in ``binary`` mode,
these data are written to a raw binary file formatted as follows. The
file starts with

* ``NWavelength`` (numpy ``int64``): the number of wavelength entries in the spectra
* ``Wavelength`` (``NWavelength`` entries of type ``double``)

and then contains a series of records, one for each output time, with
different trials ordered sequentially, so that all the times for one
trial are output before the first time for the next trial. Each record
consists of a header containing

* ``Time`` (``double``)
* ``NCluster`` (python ``int``): number of non-disrupted clusters present at this time

This is followed by ``NCluster`` entries of the following form:

* ``UniqueID`` (``unsigned long``)
* ``Incident`` (``NWavelength`` entries of numpy ``float64``)
* ``Transmitted`` (``NWavelength`` entries of numpy ``float64``)
* ``Emitted`` (``NWavelength`` entries of numpy ``float64``)
* ``Transmitted_plus_emitted`` (``NWavelength`` entries of numpy
  ``float64``)

The ``cluster_cloudyphot`` File
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This file contains data on the photometry of the spectra produced by
the interaction of the stellar radiation field with the ISM around
each cluster. It consists of a series of entries containing the
following fields:

* ``UniqueID``: a unique identifier number for each cluster that is
  preserved across times and output files
* ``Time``: evolution time at which the output is produced
* ``PhotFilter1_trans``: photometric value for the *Transmitted*
  radiation field through filter 1, where filter 1 here is the same as
  filter 1 in :ref:`ssec-int-phot-file`; units are also the same as
  in that file.
* ``PhotFilter1_emit``: photometric value for the *Emitted*
  radiation field through filter 1
* ``PhotFilter1_trans_emit``: photometric value for the
  *Transmitted_plus_emitted* radiation field through filter 1
* ``PhotFilter2_trans``
* ``PhotFilter2_emit``
* ``PhotFilter2_trans_emit``
* ``...``

For distinctions between the *Transmitted*, *Emitted*, and
*Transmitted_plus_emitted* radiation fields, see
:ref:`sssec-int-cloudyspec-file`, or the `cloudy documentation
<http://nublado.org>`_. Note that we do not record photometry for the
incident spectrum, since that would be, up to the accuracy of the
numerical integration, identical to the photometry already recorded in
the :ref:`ssec-cluster-phot-file`.

If the SLUG data input to cloudy_slug were written in ``ascii`` mode,
these data are output as a text file containing a series of columns,
with different trials separated by lines of dashes.

If the SLUG data input to cloudy_slug were written in ``fits`` mode,
these data are written in a FITS file containing one binary table
extension. The columns in this FITS file are ``Trial``, ``UniqueID``,
``Time``, ``Filter1_Transmitted``, ``Filter1_Emitted``,
``Filter1_Transmitted_plus_emitted``, ``...``. The first three columns
give the trial number, cluster unique ID, and the time, and the
remainder give the photometric values for the transmitted, emitted,
and transmitted plus emitted spectra in each filter.

If the SLUG data input to cloudy_slug were written in ``binary`` mode,
these data are written in a raw binary file that is formatted as
follows. The file starts with an ASCII text header consisting of the
following, each on a separate line:

* ``NFilter`` (stored as ``ASCII text``): number of filters used
* ``FilterName`` ``FilterUnit`` (``NFilter`` entries stored as ``ASCII
  text``): the name and units for each filter are listed in ASCII, one
  filter-unit pair per line

This is followed by a series of entries of that each begin with a
header

* ``Time`` (``double``)
* ``NCluster`` (``std::vector<double>::size_type``, usually ``unsigned long long``): number of non-disrupted clusters present at this time

This is followed by ``NCluster`` entries of the following form:

* ``UniqueID`` (``unsigned long``)
* ``PhotFilter_Transmitted`` (``NFilter`` entries of numpy
  ``float64``), giving the transmitted photometry in each filter
* ``PhotFilter_Emitted`` (``NFilter`` entries of numpy
  ``float64``), giving the emitted photometry in each filter
* ``PhotFilter_Transmitted_plus_emitted`` (``NFilter`` entries of numpy
  ``float64``), giving the transmitted plus emitted photometry in each
  filter

Full Documentation of slugpy.cloudy
-----------------------------------

.. automodule:: slugpy.cloudy
   :members:
