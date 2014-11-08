.. highlight:: rest

.. _sec-parameters:

Parameter Specification
=======================

File Format
-----------

An example parameter file is included as ``param/example.param`` in the source tree. Parameter files for SLUG are generically formatted as a series of entries of the form::

   keyword    value

Any line starting with ``#`` is considered to be a comment and is ignored, and anything on a line after a ``#`` is similarly treated as a comment and ignored. Some general rules on keywords are:

* Keywords may appear in any order.
* Some keywords have default values, indicated in parenthesis in the list below. These keywords are optional and need not appear in the parameter file. All others are required. 
* Keywords and values are case-insensitive. 
* Unless explicitly stated otherwise, units for mass are always :math:`M_\odot`, units for time are always yr.
* Any time a file or directory is specified, if it is given as a relative rather than absolute path, it is assumed to be relative to the environment variable ``$SLUG_DIR``. If this environment variable is not set, it is assumed to be relative to the current working directory.

The keywords recognized by SLUG can be categorized as described in the remainder of this section.

.. _ssec-basic-keywords:

Basic Keywords
--------------

These specify basic data for the run.

* ``model_name`` (default: ``SLUG_DEF``): name of the model. This will become the base filename for the output files.
* ``out_dir`` (default: ``output``): name of the directory into which output should be written.
* ``verbosity`` (default: ``1``): level of verbosity when running, with 0 indicating no output, 1 indicating some output, and 2 indicating a great deal of output.

Simulation Control Keywords
---------------------------

These control the operation of the simulation.

* ``sim_type`` (default: ``galaxy``): set to ``galaxy`` to run a galaxy simulation (a composite stellar population), or to ``cluster`` to run a cluster simulation (a simple stellar population)
* ``n_trials`` (default: ``1``): number of trials to run
* ``log_time`` (default: ``0``): set to 1 for logarithmic time step, 0 for linear time steps
* ``time_step``: size of the time step. If ``log_time`` is set to 0, this is in yr. If ``log_time`` is set to 1, this is in dex (i.e., a value of 0.2 indicates that every 5 time steps correspond to a factor of 10 increase in time).
* ``start_time``: first output time. This may be omitted if ``log_time`` is set to 0, in which case it defaults to a value equal to ``time_step``.
* ``end_time``: last output time, in yr. Note that not all the tracks include entries going out to times >1 Gyr, and the results will become inaccurate if the final time is larger than the tracks allow.
* ``sfr``: star formation rate. Only used if ``sim_type`` is ``galaxy``; for ``cluster``, it will be ignored, and can be omitted. If, instead of specifying a numerical value for this parameter, you specify the string ``sfh``, the code will interpret this as a flag that a star formation history should be read from the file specified by the ``sfh`` keyword.
* ``sfh``: name of star formation history file. This file is a PDF file, formatted as described in :ref:`sec-pdfs`. This is ignored, and can be omitted, if ``sim_type`` is ``cluster``, or if ``sfr`` is not set to ``sfh``.
* ``cluster_mass``: mass of the star cluster for simulations with ``sim_type`` set to ``cluster``. This can be omitted, and will be ignored, if ``sim_type`` is ``galaxy``. This parameter can be set to either a positive number or to the string ``cmf``. If it is set to a numerical value, that value will be used as the cluster mass, in :math:`M_\odot` for each trial. If it is set to ``cmf``, then a new cluster mass will be drawn from the CMF for each trial.
* ``redshift`` (default: ``0``): place the system at the specified redshift. The computed spectra and photometry will then be computed in the observed rather than the rest frame of the system.


Output Control Keywords
-----------------------

These control what quantities are computed and written to disk. Full a full description of the output files and how they are formatted, see :ref:`sec-output`.

* ``out_cluster`` (default: ``1``): write out the physical properties of star clusters? Set to 1 for yes, 0 for no.
* ``out_cluster_phot`` (default: ``1``): write out the photometry of star clusters? Set to 1 for yes, 0 for no.
* ``out_cluster_spec`` (default: ``1``): write out the spectra of star clusters? Set to 1 for yes, 0 for no.
* ``out_integrated`` (default: ``1``): write out the integrated physical properties of the whole galaxy? Set to 1 for yes, 0 for no. This keyword is ignored if ``sim_type`` is ``cluster``.
* ``out_integrated_phot`` (default: ``1``): write out the integrated photometry of the entire galaxy? Set to 1 for yes, 0 for no. This keyword is ignored if ``sim_type`` is ``cluster``.
* ``out_integrated_spec`` (default: ``1``): write out the integrated spectra of the entire galaxy? Set to 1 for yes, 0 for no. This keyword is ignored if ``sim_type`` is ``cluster``.
* ``output_mode`` (default: ``ascii``): set to ``ascii``, ``binary``, or ``fits``. Selecting ``ascii`` causes the output to be written in ASCII text, which is human-readable, but produces much larger files. Selecting ``binary`` causes the output to be written in raw binary. Selecting ``fits`` causes the output to be written FITS format. This will be somewhat larger than raw binary output, but the resulting files will be portable between machines, which the raw binary files are not guaranteed to be. All three output modes can be read by the python library, though with varying speed -- ASCII output is slowest, FITS is intermediate, and binary is fastest.

.. _ssec-phys-keywords:

Physical Model Keywords
-----------------------

These specify the physical models to be used for stellar evolution, atmospheres, the IMF, extinction, etc.

* ``imf`` (default: ``lib/imf/chabrier.imf``): name of the IMF descriptor file; this is a PDF file, formatted as described in :ref:`sec-pdfs`. Note that SLUG ships with the following IMF files pre-defined (in the directory ``lib/imf``)
   * ``chabrier.imf`` (single-star IMF from `Chabrier, 2005, in "The Initial Mass Function 50 Years Later", eds. E. Corbelli, F. Palla, & H. Zinnecker, Springer: Dordrecht, p. 41 <http://adsabs.harvard.edu/abs/2005ASSL..327...41C>`_)
   * ``chabrier03.imf`` (single-star IMF from `Chabrier, 2003, PASP, 115, 763-795 <http://adsabs.harvard.edu/abs/2003PASP..115..763C>`_)
   * ``kroupa.imf`` (IMF from `Kroupa, 2002, Science, 295, 82-91 <http://adsabs.harvard.edu/abs/2002Sci...295...82K>`_)
   * ``kroupa_sb99.imf`` (simplified version of the Kroupa, 2002 IMF used by default by `starburst99 <http://www.stsci.edu/science/starburst99/docs/default.htm>`_)
   * ``salpeter.imf`` (single-component power law IMF from `Salpeter, 1955, ApJ, 121, 161 <http://adsabs.harvard.edu/abs/1955ApJ...121..161S>`_)
* ``cmf`` (default: ``lib/cmf/slug_default.cmf``): name of the CMF descriptor file; this is a PDF file, formatted as described in :ref:`sec-pdfs`. The default selection is a power law :math:`dN/dM \propto M^{-2}` from :math:`M = 10^2 - 10^7\;M_\odot`. This is ignored, and may be omitted, if ``sim_type`` is set to ``cluster`` and ``cluster_mass`` is set to a numerical value.
* ``clf`` (default: ``lib/clf/slug_default.clf``): name of the CLF descriptor file; this is a PDF file, formatted as described in :ref:`sec-pdfs`. The default gives a power law distribution of lifetimes :math:`t` with :math:`dN/dt\propto t^{-1.9}` from 1 Myr to 1 Gyr. Note that this corresponds to a cluster age distribution of slope -0.9. The SLUG source also ships with an alternative CLF file, ``lib/clf/nodisrupt.clf``, which disables cluster disruption entirely (by setting the lifetime distribution to a :math:`\delta` function at :math:`10^{300}` yr).
* ``tracks`` (default: ``lib/tracks/Z0140v00.txt``): stellar evolution tracks to use. The following tracks ship with SLUG (all in the directory ``lib/tracks``):
   * ``ZXXXXvYY.txt``: Geneva (2013) tracks; metallicities are Solar (``XXXX = 0140``) and 1/7 Solar (``XXXX = 0020``), and rotation rates are 0 (``YY = 00``) and 40% of breakup (``YY = 40``).
   * ``modcXXX.dat``: Geneva tracks with standard mass loss, for metallicities of :math:`2\times` Solar (``040``), Solar (``020``), :math:`0.4\times` Solar (``008``), :math:`0.2\times` Solar (``004``), and :math:`0.05\times` Solar (``001``).
   * ``modeXXX.dat``: same as ``modcXXX.dat``, but with higher mass loss rates.
   * ``modpXXX.dat``: Padova tracks with thermally pulsing AGB stars; metallicities use the same scale as ``modcXXX.dat`` files (i.e., ``020`` is Solar).
   * ``modsXXX.dat``: same as ``modpXXX.dat``, but without thermally pulsing AGB stars
* ``atmospheres`` (default: ``lib/atmospheres``): directory where the stellar atmosphere library is located. Note that file names are hard-coded, so if you want to use different atmosphere models with a different format, you will have to write new source code to do so.
* ``specsyn_mode`` (default: ``sb99``): spectral synthesis mode. Allowed values are:
   * ``planck``: treat all stars as black bodies
   * ``Kurucz``: use Kurucz atmospheres, as compiled by `Lejeune et al. (1997, A&AS, 125, 229) <http://adsabs.harvard.edu/abs/1997A%26AS..125..229L>`_, for all stars
   * ``Kurucz+Hillier``: use Kurucz atmospheres for all stars except Wolf-Rayet stars; WR stars use Hillier model atmospheres (`Hillier & Miller, 1998, ApJ, 496, 407 <http://adsabs.harvard.edu/abs/1998ApJ...496..407H>`_)
   * ``Kurucz+Pauldrach``: use Kurucz atmospheres for all stars except OB stars; OB stars use Pauldrach model atmospheres (`Pauldrach et al., 2001, A&A, 375, 161 <http://adsabs.harvard.edu/abs/2001A%26A...375..161P>`_)
   * ``SB99``: emulate the behavior of ``starburst99``: use Pauldrach for OB stars, Hillier for WR stars, and Kurucz for all other stars
* ``clust_frac`` (default: ``1.0``): fraction of stars formed in clusters
* ``min_stoch_mass`` (default: ``0.0``): minimum stellar mass to be treated stochastically. All stars with masses below this value are assumed to be sampled continuously from the IMF.
* ``metallicity``: metallicity of the stellar population, relative to solar. This may be omitted if ``tracks`` is set to one of the default sets of tracks that ships with SLUG, as the metallicities for these tracks are hardwired in. This keyword is provided to allow users to supply their own tracks.
* ``WR_mass``: minimum starting mass that stars must have in order to pass through a Wolf-Rayet phase. This can be omitted if ``tracks`` is set to one of the default sets of tracks that ships with SLUG, as the WR cutoff masses for these tracks are hardwired in. This keyword is provided to allow users to supply their own tracks.
* ``A_V`` (default: no extinction): extinction distribution. This parameter has three possible behaviors. If the parameter ``A_V`` is omitted entirely, then the code will not compute extinction-corrected spectra or photometry at all; only unextincted values will be reported. If this parameter is specified as a real number, it will be interepreted as specifying a uniform extinction value :math:`A_V`, in mag, and this extinction will be applied to all predicted light output. Finally, if this parameter is a string that cannot be converted to a real number, it will be interpreted as the name of a PDF file, formatted as described in :ref:`sec-pdfs`, specifying the probability distribution of :math:`A_V` values, in mag.
* ``extinction_curve`` (default: ``lib/extinct/SB_ATT_SLUG.dat``) file specifying the extinction curve; the file format is two columns of numbers in ASCII, the first giving the wavelength in Angstrom and the second giving the exintction :math:`\kappa_\nu` at that wavelength / frequency in :math:`\mathrm{cm}^2`. Note that the absolute normalization of the exitnction curve is unimportant; only the wavelength-dependence matters (see :ref:`ssec-spec-phot`). SLUG ships with the following extinction curves (all in ``lib/extinct``):
   * ``LMC_EXT_SLUG.dat`` : LMC extinction curve; optical-UV from `Fitzpatrick, E. L., 1999, PASP, 111, 63 <http://adsabs.harvard.edu/abs/1999PASP..111...63F>`_, IR from `Landini, M., et al., 1984, A&A, 134, 284 <http://adsabs.harvard.edu/abs/1984A%26A...134..284L>`_; parts combined by D. Calzetti
   * ``MW_EXT_SLUG.dat`` : MW extinction curve; optical-UV from `Fitzpatrick, E. L., 1999 PASP, 111, 63 <http://adsabs.harvard.edu/abs/1999PASP..111...63F>`_, IR from `Landini, M., et al., 1984, A&A, 134, 284 <http://adsabs.harvard.edu/abs/1984A%26A...134..284L>`_; parts combined by D. Calzetti
   * ``SB_ATT_SLUG.dat`` : "starburst" extinction curve from `Calzetti, D., et al., 2000, ApJ, 533, 682 <http://adsabs.harvard.edu/abs/2000ApJ...533..682C>`_
   * ``SMC_EXT_SLUG.dat`` : SMC extinction curve from `Bouchet, P., et al., 1985, A&A, 149, 330 <http://adsabs.harvard.edu/abs/1985A%26A...149..330B>`_

.. _ssec-phot-keywords:

Photometric Filter Keywords
---------------------------

These describe the photometry to be computed. Note that none of these keywords have any effect unless ``out_integrated_phot`` or ``out_cluster_phot`` is set to 1.

* ``phot_bands``: photometric bands for which photometry is to be computed. The values listed here can be comma- or whitespace-separated. For a list of available photometric filters, see the file ``lib/filters/FILTER_LIST``. In addition to these filters, SLUG always allows four special "bands":
   * ``QH0``: the :math:`\mathrm{H}^0` ionizing luminosity, in photons/sec
   * ``QHe0``: the :math:`\mathrm{He}^0` ionizing luminosity, in photons/sec
   * ``QHe1``: the :math:`\mathrm{He}^+` ionizing luminosity, in photons/sec
   * ``Lbol``: the bolometric luminosity, in :math:`L_\odot`
* ``filters`` (default: ``lib/filters``): directory containing photometric filter data
* ``phot_mode`` (default: ``L_nu``): photometric system to be used when writing photometric outputs. Full definitions of the quantities computed for each of the choices listed below are given in :ref:`ssec-spec-phot`. Note that these values are ignored for the four special bands ``QH0``, ``QHe0``, ``QHe1``, and ``Lbol``. These four bands are always written out in the units specified above. Allowed values are:
   * ``L_nu``: report frequency-averaged luminosity in the band, in units of erg/s/Hz
   * ``L_lambda``: report wavelength-averaged luminosity in the band, in units of erg/s/Angstrom
   * ``AB``: report AB magnitude
   * ``STMAG``: report ST magnitude
   * ``VEGA``: report Vega magnitude


