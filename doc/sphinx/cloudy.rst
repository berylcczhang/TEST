.. highlight:: rest

.. _sec-cloudy-slug:

cloudy_slug: An Automated Interface to cloudy
=============================================

SLUG stochastically generates stellar spectra, but it does not compute
the nebular lines produced when those photons interact with the
interstellar medium. To perform such calculations, SLUG includes an
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

5. Invoke the cloudy_slug interface script via::

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
   :ref:`sec-slugpy`). See :ref:`ssec-cloudy-output` for a description
   of the outputs.

.. _ssec-cloudy-cluster:

The cloudy_slug Physical Model: Integrated Mode Versus Cluster Mode
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
required by cloudy. The cloudy_slug interface can compute nebular
emission under either of these scenarios; we refer to the former as
integrated mode, and to the latter as cluster mode. Note that, in
either mode, the spectrum that is used to compute the nebular emission
will be the *unextincted, non-redshifted* spectrum computed by SLUG.

.. _sssec-cloudy-integrated-mode:

Integrated Mode
^^^^^^^^^^^^^^^

In integrated mode, cloudy_slug will read all the spectra contained in
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

.. _sssec-cloudy-cluster-mode:

Cluster Mode
^^^^^^^^^^^^

In cluster mode, cloudy_slug will read all the individual cluster
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
by hydrogen atoms rather than dust, :math:`\epsilon_0 =
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

Once the outer radius is known, cloudy_slug sets the starting radius
for the cloudy calculation to :math:`10^{-3} r_{\mathrm{II}}`, and
sets the starting density to the value expected for photoionization
equilibrum in a uniform HII region,

.. math:: n_{\mathrm{II}} = \left(\frac{3
	  Q(\mathrm{H}^0)}{4.4 \pi\alpha_B
	  r_{\mathrm{II}}^3}\right)^{1/2}

Note that this approximation will be highly inaccurate if
:math:`r_{\mathrm{II}} \ll r_{\mathrm{ch}}`, but no better analytic
approximation is available, and this phase should be very short-lived
for most clusters.


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
#. The input file *must not* contain a radius command. This too will
   be computed automatically by the cloudy_slug script.
#. The input file *must* contain an entry ``hden N`` where ``N`` is
   the log base 10 of the hydrogen density. This will be interpreted
   differently depending on whether cloudy_slug is being run in
   cluster mode or integrated mode -- see :ref:`ssec-cloudy-cluster`.
#. Any outputs to be written (specified using the ``save`` or
   ``punch`` keywords) must give file names containing the string
   ``OUTPUT_FILENAME``. This string will be replaced by the
   cloudy_slug script to generate a unique file name for each cloudy
   run, and to read back these outputs for post-processing.
#. The cloudy_slug output will contain output spectra only if the
   cloudy input file contains a ``save last continuum`` command. See
   :ref:`ssec-cloudy-output`.
#. The cloudy_slug output will contain output line luminosities only
   if the cloudy input file contains a ``save last line list emergent
   absolute column`` command. See :ref:`ssec-cloudy-output`.
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
                  [--cloudytemplate CLOUDYTEMPLATE] [-cm] [-nl NICELEVEL]
                  [-n NPROC] [-s] [--slugpath SLUGPATH] [-v]
                  slug_model_name [start_spec] [end_spec]

The positional arguments are as follows:

* ``slug_model_name``: this is the name of the SLUG output to be used
  as a basis for the cloudy calculation. This should be the same as
  the ``model_name`` parameter used in the SLUG simulation, with the
  optional addition of a path specification in front.
* ``start_spec``: default behavior is to run cloudy on all the
  integrated spectra (in :ref:`sssec-cloudy-integrated-mode`) or
  cluster spectra (in :ref:`sssec-cloudy-cluster-mode`). If this
  argument is set, cloudy will only be run in spectra starting with
  the specified trial number (in :ref:`sssec-cloudy-integrated-mode`)
  or cluster number (in :ref:`sssec-cloudy-cluster-mode`); numbers are
  0-offset, to the first trial/cluster is 0, the next is 1, etc.
* ``end_spec``: default behavior is to run cloudy on all the
  integrated spectra (in :ref:`sssec-cloudy-integrated-mode`) or
  cluster spectra (in :ref:`sssec-cloudy-cluster-mode`). If this
  argument is set, cloudy will only be run on spectra up to the
  specified trial number (in :ref:`sssec-cloudy-integrated-mode`) or
  cluster number (in :ref:`sssec-cloudy-cluster-mode`); numbers are
  0-offset, to the first trial/cluster is 0, the next is 1, etc.

The optional arguments are as follows:

* ``-h, --help``: prints a help message and then exits
* ``-a AGEMAX, --agemax AGEMAX``: maximum cluster age in Myr for
  cloudy computation. Cloudy will not be run on clusters older than
  this value, and the predicted nebular emission for such clusters
  will be recorded as zero. Default value is 4 Myr. This argument only
  has an effect if running in :ref:`sssec-cloudy-cluster-mode`;
  otherwise it is ignored.
* ``--cloudypath CLOUDYPATH``: path to the cloudy executable; default
  is ``$CLOUDY_DIR/cloudy.exe``
* ``--cloudytemplate CLOUDYTEMPLATE``: cloudy input file template (see
  :ref:`ssec-cloudy-template`); default is
  ``$SLUG_DIR/cloudy_slug/cloudy.in_template``
* ``-cm, --clustermode``: if this argument is set, then cloudy_slug
  will run in :ref:`sssec-cloudy-cluster-mode`; default behavior is to
  run in :ref:`sssec-cloudy-integrated-mode`
* ``-nl NICELEVEL, --nicelevel NICELEVEL``: if this is set, then the
  cloudy processes launched by the script will be run at this nice
  level. If it is not set, they will not be nice'd. Note that this
  option will only work correctly on platforms that support nice.
* ``-n NPROC, --nproc NPROC``: number of simultaneous cloudy processes
  to run; default is the number of cores available on the system
* ``-s, --save``: by default, cloudy_slug will extract line and
  spectral data from the cloudy outputs and store them as described in
  :ref:`ssec-cloudy-output`, then delete the cloudy output files. If
  this option is set, the cloudy output files will NOT be deleted, and
  will be left in place. WARNING: cloudy's outputs are written in
  ASCII and are quite voluminous, so only choose this option if you
  are only running cloudy on a small number of SLUG spectra and/or you
  are prepared to store hundreds of GB more more.
* ``--slugpath SLUGPATH``: path to the SLUG output data. If not set,
  cloudy_slug searches for an appropriately-named set of output files
  first in the current working directory, and next in
  ``$SLUG_DIR/output``
* ``-v, --verbose``: if this option is set, cloudy_slug produces
  verbose output as it runs

.. _ssec-cloudy-output:

Full Description of cloudy_slug Output
--------------------------------------

The cloudy_slug script will automatically process the cloudy output
and produce a series of new output files, which will be written to the
same directory where the input SLUG files are located, and using the
same output mode (ASCII text, raw binary, or FITS -- see
:ref:`sec-output`). If cloudy_slug is run in
:ref:`sssec-cloudy-integrated-mode`, the three output files will be
``MODEL_NAME_integrated_cloudylines.ext``, 
``MODEL_NAME_integrated_cloudyphot.ext``, and 
``MODEL_NAME_integrated_cloudyspec.ext``, where the extension ``.ext``
is one of ``.txt``, ``.bin``, or ``.fits``, depending on the
``output_mode``. If cloudy_slug is run in
:ref:`sssec-cloudy-cluster-mode`, the three output files will be
``MODEL_NAME_cluster_cloudylines.ext``, 
``MODEL_NAME_cluster_cloudyphot.ext``, and 
``MODEL_NAME_cluster_cloudyspec.ext``. All of these output files will
be read and processed automatically if the outputs are read using
``read_integrated`` or ``read_cluster`` in the :ref:`sec-slugpy`
library.

The format of those files is described below.

The ``integrated_cloudylines`` File
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This file contains data on the nebular line emission produced by the
interaction of the stellar radiation field with the ISM. It consists
of a series of entries containing the following fields:

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

* ``Time`` (``double``)
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

* ``Time`` (numpy ``float64``)
* ``Incident`` (``NWavelength`` entries of numpy ``float64``)
* ``Transmitted`` (``NWavelength`` entries of numpy ``float64``)
* ``Emitted`` (``NWavelength`` entries of numpy ``float64``)
* ``Transmitted_plus_emitted`` (``NWavelength`` entries of numpy
  ``float64``)

There is one such record for each output time, with different trials ordered sequentially, so that all the times for one trial are output before the first time for the next trial.

The ``integrated_cloudyphot`` File
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This file contains photometric data computed for the spectra produced
by the interaction between the stellar radiation field and the HII
region. The file consists of a series of entries containing the
following fields:

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

This is followed by a series of records, one for each output time, with different trials ordered sequentially, so that all the times for one trial are output before the first time for the next trial. Each record consists of a header containing

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
:ref:`sssec-int-cloudyspec-file`, or the `cloudy documentaiton
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
