.. highlight:: rest

.. _sec-parameters:

Parameter Specification
=======================

Automated Parameter File Generation
-----------------------------------

The remainder of this section contains information on how parameter files are formatted, and exactly how parameter choices specify code behavior. However, as a convenience slug comes with a python script that provides a simple menu-driven interface to write parameter files automatically. The script can be started by doing::

   python tools/python/write_param.py

Once started, the script provides a series of menus that allow the user to set all the keywords specified below. The script can then write a validly-formatted parameter file based on the options chosen.


File Format
-----------

An example parameter file is included as ``param/example.param`` in the source tree. Parameter files for slug are generically formatted as a series of entries of the form::

   keyword    value

Any line starting with ``#`` is considered to be a comment and is ignored, and anything on a line after a ``#`` is similarly treated as a comment and ignored. Some general rules on keywords are:

* Keywords may appear in any order.
* Many keywords have default values, indicated in parentheses in the
  list below. These keywords are optional and need not appear in the
  parameter file; missing keywords are set to their default values.
  Keywords that do not have a default are required.
* Keywords and values are case-insensitive. 
* Unless explicitly stated otherwise, units for mass are always
  :math:`M_\odot`, units for time are always yr, and units for
  extinction are always magnitudes.
* Any time a file or directory is specified, the path is resolved by
  the following rules:
  
  1. Absolute paths are always treated exactly as written.
  2. For relative paths, the code first checks for a file or directory
     of that name in the directory where slug is running. If no file
     of that name is found, and the environment variable ``$SLUG_DIR``
     is set, teh code will search for the file at
     ``$SLUG_DIR/user_specified_path``. This is to allow users to
     use default options that ship with the code simply by setting
     their ``$SLUG_DIR`` environment variable to wherever the code is
     installed.
  3. Notwithstanding rule 2, the output directory (``out_dir``) is
     always assumed to be relative to the location where slug is
     running; it will never be prepended with ``$SLUG_DIR``.

The keywords recognized by slug can be categorized as described in the remainder of this section.

.. _ssec-basic-keywords:

Basic Keywords
--------------

These specify basic data for the run.

* ``model_name`` (default: ``SLUG_DEF``): name of the model. This will become the base filename for the output files.
* ``out_dir`` (default: current working direcory): name of the directory into which output should be written. If not specified, output is written into the directory from which the slug executable is called.
* ``verbosity`` (default: ``1``): level of verbosity when running, with 0 indicating no output, 1 indicating some output, and 2 indicating a great deal of output.

Simulation Control Keywords
---------------------------

These control the operation of the simulation.

* ``sim_type`` (default: ``galaxy``): set to ``galaxy`` to run a galaxy simulation (a composite stellar population), or to ``cluster`` to run a cluster simulation (a simple stellar population)
* ``n_trials`` (default: ``1``): number of trials to run
* ``checkpoint_interval`` (default: checkpointing off): output a checkpoint every ``checkpoint_interval`` trials
* ``log_time`` (default: ``0``): set to 1 for logarithmic time step, 0 for linear time steps
* ``time_step``: size of the time step. If ``log_time`` is set to 0, this is in yr. If ``log_time`` is set to 1, this is in dex (i.e., a value of 0.2 indicates that every 5 time steps correspond to a factor of 10 increase in time). Alternately, if ``time_step`` is set to any value that cannot be converted to a real number, then this is interpreted as giving the name of a PDF file, which must be formatted as described in :ref:`sec-pdfs`. In this case one output time will be selected randomly for each trial from the specified PDF. This option is useful, for example, for generating a library of simulations that are randomly sampled in stellar population age. For the PDF option, the options ``log_time``, ``start_time`` and ``end_time`` will all be ignored, as the relevant parameters will be taken from the specified PDF file. This keyword may be omitted, and will be ignored, if ``output_times`` is set.
* ``start_time``: first output time. This may be omitted if ``log_time`` is set to 0, in which case it defaults to a value equal to ``time_step``. It may also be omitted if ``output_times`` is set.
* ``end_time``: last output time, in yr. This may be omitted of ``output_times`` is set. Note that not all the tracks include entries going out to times >1 Gyr, and the results will become inaccurate if the final time is larger than the tracks allow.
* ``output_times``: an optional parameter giving an exact list of times in yr at which to write output. This must be specified as a comma-separated list, i.e., time0, time1, time2, ..., where all times are positive and are in strictly increasing order. There is no limit on the number of output times that may be included. This parameter may be omitted if ``end_time`` and ``time_step`` are set. If this parameter is set, it overrides ``start_time``, ``end_time``, and ``time_step``.
* ``sfr``: star formation rate. Only used if ``sim_type`` is ``galaxy``; for ``cluster``, it will be ignored, and can be omitted. This parameter can be set in three ways. If the value given is a number, it will be interpreted as specifying a constant star formation rate. If it is the string ``sfh``, the code will interpret this as a flag that a star formation history should be read from the file specified by the ``sfh`` keyword. If the parameter value is any other string that cannot be converted to a numerical value, it will be interpreted as the name of a PDF file (see :ref:`sec-pdfs`); a (constant) value of the star formation rate for each trial will be drawn from this PDF.
* ``sfh``: name of star formation history file. This file is a PDF file, formatted as described in :ref:`sec-pdfs`. This is ignored, and can be omitted, if ``sim_type`` is ``cluster``, or if ``sfr`` is not set to ``sfh``.
* ``cluster_mass``: mass of the star cluster for simulations with ``sim_type`` set to ``cluster``. This can be omitted, and will be ignored, if ``sim_type`` is ``galaxy``. This parameter can be set to either a positive number or to the string ``cmf``. If it is set to a numerical value, that value will be used as the cluster mass, in :math:`M_\odot` for each trial. If it is set to ``cmf``, then a new cluster mass will be drawn from the CMF for each trial.
* ``redshift`` (default: ``0``): place the system at the specified redshift. The computed spectra and photometry will then be computed in the observed rather than the rest frame of the system.


Output Control Keywords
-----------------------

These control what quantities are computed and written to disk. Full a full description of the output files and how they are formatted, see :ref:`sec-output`.

* ``out_cluster`` (default: ``1``): write out the physical properties of star clusters? Set to 1 for yes, 0 for no.
* ``out_cluster_phot`` (default: ``1``): write out the photometry of star clusters? Set to 1 for yes, 0 for no.
* ``out_cluster_spec`` (default: ``1``): write out the spectra of star clusters? Set to 1 for yes, 0 for no.
* ``out_cluster_yield`` (default: ``1``): write out the yield of star clusters? Set to 1 for yes, 0 for no.
* ``out_integrated`` (default: ``1``): write out the integrated physical properties of the whole galaxy? Set to 1 for yes, 0 for no. This keyword is ignored if ``sim_type`` is ``cluster``.
* ``out_integrated_phot`` (default: ``1``): write out the integrated photometry of the entire galaxy? Set to 1 for yes, 0 for no. This keyword is ignored if ``sim_type`` is ``cluster``.
* ``out_integrated_spec`` (default: ``1``): write out the integrated spectra of the entire galaxy? Set to 1 for yes, 0 for no. This keyword is ignored if ``sim_type`` is ``cluster``.
* ``out_integrated_yield`` (default: ``1``): write out the integrated yield of the entire galaxy? Set to 1 for yes, 0 for no. This keyword is ignored if ``sim_type`` is ``cluster``.
* ``output_mode`` (default: ``ascii``): set to ``ascii``, ``binary``, or ``fits``. Selecting ``ascii`` causes the output to be written in ASCII text, which is human-readable, but produces much larger files. Selecting ``binary`` causes the output to be written in raw binary. Selecting ``fits`` causes the output to be written FITS format. This will be somewhat larger than raw binary output, but the resulting files will be portable between machines, which the raw binary files are not guaranteed to be. All three output modes can be read by the python library, though with varying speed -- ASCII output is slowest, FITS is intermediate, and binary is fastest.

.. _ssec-stellar-keywords:

Stellar Model Keywords
----------------------

These specify the physical models to be used for stellar evolution, atmospheres, the IMF, extinction, etc.

* ``imf`` (default: ``lib/imf/chabrier.imf``): name of the IMF descriptor file; this is a PDF file, formatted as described in :ref:`sec-pdfs`. Note that slug ships with the following IMF files pre-defined (in the directory ``lib/imf``)
   * ``chabrier.imf`` (single-star IMF from `Chabrier, 2005, in "The Initial Mass Function 50 Years Later", eds. E. Corbelli, F. Palla, & H. Zinnecker, Springer: Dordrecht, p. 41 <http://adsabs.harvard.edu/abs/2005ASSL..327...41C>`_)
   * ``chabrier03.imf`` (single-star IMF from `Chabrier, 2003, PASP, 115, 763-795 <http://adsabs.harvard.edu/abs/2003PASP..115..763C>`_)
   * ``kroupa.imf`` (IMF from `Kroupa, 2002, Science, 295, 82-91 <http://adsabs.harvard.edu/abs/2002Sci...295...82K>`_)
   * ``kroupa_sb99.imf`` (simplified version of the Kroupa, 2002 IMF used by default by `starburst99 <http://www.stsci.edu/science/starburst99/docs/default.htm>`_)
   * ``salpeter.imf`` (single-component power law IMF from `Salpeter, 1955, ApJ, 121, 161 <http://adsabs.harvard.edu/abs/1955ApJ...121..161S>`_)
* ``cmf`` (default: ``lib/cmf/slug_default.cmf``): name of the CMF descriptor file; this is a PDF file, formatted as described in :ref:`sec-pdfs`. The default selection is a power law :math:`dN/dM \propto M^{-2}` from :math:`M = 10^2 - 10^7\;M_\odot`. This is ignored, and may be omitted, if ``sim_type`` is set to ``cluster`` and ``cluster_mass`` is set to a numerical value.
* ``clf`` (default: ``lib/clf/slug_default.clf``): name of the CLF descriptor file; this is a PDF file, formatted as described in :ref:`sec-pdfs`. The default gives a power law distribution of lifetimes :math:`t` with :math:`dN/dt\propto t^{-1.9}` from 1 Myr to 1 Gyr. Note that this corresponds to a cluster age distribution of slope -0.9. The slug source also ships with an alternative CLF file, ``lib/clf/nodisrupt.clf``, which disables cluster disruption entirely (by setting the lifetime distribution to a :math:`\delta` function at :math:`10^{300}` yr).
* ``tracks`` (default: ``geneva_2013_vvcrit_00``): stellar evolution tracks to use. This can be specified either by giving the name of a particular set of tracks (i.e., a set of tracks computed using the same code or group, but at a range of metallicities), or by giving the name of a particular file (a particular set of tracks at a particular metallicity). When specifying a track set, the user can also specify the metallicity (see below), and the tracks will be read for (and if necessary interpolated to) the specified metallicity automatically. The following track sets and indvidual track files are available:
   * ``geneva_2013_vvcrit_00`` and ``geneva_2013_vvcrit_00``: Geneva (2013) track set, rotating at 0% and 40% of breakup, respectively. These tracks are available at metallicities of Solar and 1/7 Solar.
   * ``geneva_mdot_std`` and ``geneva_mdot_enhanced``: pre-2013 Geneva track set, no rotation, with standard and 2 times standard mass loss rates, respectively. These models are available in metallicities of (relative to Solar) :math:`Z = 0.05, 0.2, 0.4, 1.0, 2.0`.
   * ``padova_tpagb_yes`` and ``padova_tpagb_no``: Padova tracks, no rotation, with and withouth thermally-pulsing AGB stars. These models are available in metallicities of (relative to Solar) :math:`Z = 0.02, 0.2, 0.4, 1.0, 2.5`.
   * ``mist_2016_vvcrit_00`` and ``mist_2016_vvcrit_40``: MIST v1.0 models, rotating and 0% and 40% of breakup, respectively; these models are available at Solar-normalised metallcities of :math:`\log Z = -4` to 0.5, in steps of 0.5 dex from :math:`\log Z = -4` to :math:`\log Z = -2`, and 0.25 dex from :math:`\log Z = -2` to :math:`\log Z = 0.5`.
   * ``lib/tracks/sb99/ZXXXXvYY.txt``: individual files for Geneva (2013) tracks; metallicities are Solar (``XXXX = 0140``) and 1/7 Solar (``XXXX = 0020``), and rotation rates are 0 (``YY = 00``) and 40% of breakup (``YY = 40``).
   * ``lib/tracks/sb99/modcXXX.dat``: individual files Geneva tracks with standard mass loss, for metallicities of :math:`2\times` Solar (``040``), Solar (``020``), :math:`0.4\times` Solar (``008``), :math:`0.2\times` Solar (``004``), and :math:`0.05\times` Solar (``001``).
   * ``lib/tracks/sb99/modeXXX.dat``: same as ``modcXXX.dat``, but twice standard mass loss rates.
   * ``lib/tracks/sb99/modpXXX.dat``: individual files for Padova tracks with thermally pulsing AGB stars; metallicities use the same scale as ``modcXXX.dat`` files (i.e., ``020`` is Solar).
   * ``lib/tracks/sb99/modsXXX.dat``: same as ``modpXXX.dat``, but without thermally pulsing AGB stars
   * ``lib/tracks/mist/vvcrit000/MIST_v1.0_feh_XXXXX_afe_p0.0_vvcrit0.0_EEPS.fits.gz``: individual files for MIST (2016, v1.0) non-rotating tracks. The ``XXXXX`` specifies the metallicity; the first letter is ``p`` or ``m`` for plus or minus, and the following letters give the numerical value of the log metallicity in Solar-scaled units (e.g., ``p0.00`` is Solar, ``m1.00`` is 1/10 solar, ``m2.00`` is 1/100th Solar, etc.).
   * ``lib/tracks/mist/vvcrit040//MIST_v1.0_feh_XXXXX_afe_p0.0_vvcrit0.4_EEPS.fits.gz``: same as ``/MIST_v1.0_feh_XXXXX_afe_p0.0_vvcrit0.0_EEPS.fits.gz``, but rotating at 40% of breakup
* ``atmospheres`` (default: ``lib/atmospheres``): directory where the stellar atmosphere library is located. Note that file names are hard-coded, so if you want to use different atmosphere models with a different format, you will have to write new source code to do so.
* ``specsyn_mode`` (default: ``sb99``): spectral synthesis mode. Allowed values are:
   * ``planck``: treat all stars as black bodies
   * ``Kurucz``: use Kurucz atmospheres, as compiled by `Lejeune et al. (1997, A&AS, 125, 229) <http://adsabs.harvard.edu/abs/1997A%26AS..125..229L>`_, for all stars
   * ``Kurucz+Hillier``: use Kurucz atmospheres for all stars except Wolf-Rayet stars; WR stars use Hillier model atmospheres (`Hillier & Miller, 1998, ApJ, 496, 407 <http://adsabs.harvard.edu/abs/1998ApJ...496..407H>`_)
   * ``Kurucz+Pauldrach``: use Kurucz atmospheres for all stars except OB stars; OB stars use Pauldrach model atmospheres (`Pauldrach et al., 2001, A&A, 375, 161 <http://adsabs.harvard.edu/abs/2001A%26A...375..161P>`_)
   * ``SB99``: emulate the behavior of ``starburst99``: use Pauldrach for OB stars, Hillier for WR stars, and Kurucz for all other stars
* ``clust_frac`` (default: ``1.0``): fraction of stars formed in clusters
* ``min_stoch_mass`` (default: ``0.0``): minimum stellar mass to be treated stochastically. All stars with masses below this value are assumed to be sampled continuously from the IMF.
* ``metallicity`` (default: ``1.0``): metallicity of the stellar population, relative to Solar. If the tracks are specified by giving a track set, this value must be within the metallicity range covered by the chosen track set. If the tracks are set by specifying a particular track file, this keyword will be ignored in favor of the metallicity used for that track file, and a warning will be issued if it is set.

.. _ssec-extinction-keywords:

Extinction Keywords
-------------------

* ``A_V`` (default: no extinction): extinction distribution. This parameter has three possible behaviors. If the parameter ``A_V`` is omitted entirely, then the code will not compute extinction-corrected spectra or photometry at all; only unextincted values will be reported. If this parameter is specified as a real number, it will be interepreted as specifying a uniform extinction value :math:`A_V`, in mag, and this extinction will be applied to all predicted light output. Finally, if this parameter is a string that cannot be converted to a real number, it will be interpreted as the name of a PDF file, formatted as described in :ref:`sec-pdfs`, specifying the probability distribution of :math:`A_V` values, in mag.
* ``extinction_curve`` (default: ``lib/extinct/SB_ATT_SLUG.dat``) file specifying the extinction curve; the file format is two columns of numbers in ASCII, the first giving the wavelength in Angstrom and the second giving the exintction :math:`\kappa_\nu` at that wavelength / frequency in :math:`\mathrm{cm}^2`. Note that the absolute normalization of the exitnction curve is unimportant; only the wavelength-dependence matters (see :ref:`ssec-spec-phot`). Slug ships with the following extinction curves (all in ``lib/extinct``):
   * ``LMC_EXT_SLUG.dat`` : LMC extinction curve; optical-UV from `Fitzpatrick, E. L., 1999, PASP, 111, 63 <http://adsabs.harvard.edu/abs/1999PASP..111...63F>`_, IR from `Landini, M., et al., 1984, A&A, 134, 284 <http://adsabs.harvard.edu/abs/1984A%26A...134..284L>`_; parts combined by D. Calzetti
   * ``MW_EXT_SLUG.dat`` : MW extinction curve; optical-UV from `Fitzpatrick, E. L., 1999 PASP, 111, 63 <http://adsabs.harvard.edu/abs/1999PASP..111...63F>`_, IR from `Landini, M., et al., 1984, A&A, 134, 284 <http://adsabs.harvard.edu/abs/1984A%26A...134..284L>`_; parts combined by D. Calzetti
   * ``SB_ATT_SLUG.dat`` : "starburst" extinction curve from `Calzetti, D., et al., 2000, ApJ, 533, 682 <http://adsabs.harvard.edu/abs/2000ApJ...533..682C>`_
   * ``SMC_EXT_SLUG.dat`` : SMC extinction curve from `Bouchet, P., et al., 1985, A&A, 149, 330 <http://adsabs.harvard.edu/abs/1985A%26A...149..330B>`_
   * ``MW_draine_RV3.1.dat`` : MW extinction curve for reddening :math:`R_V = 3.1`, taken from the model of `Draine, 2003, ARA&A, 41, 241 <http://adsabs.harvard.edu/abs/2003ARA%26A..41..241D>`_, and retrieved from B. Draine's `personal web page <https://www.astro.princeton.edu/~draine/dust/dustmix.html>`_
* ``nebular_extinction_factor`` (default: 1.0): nebular extinction excess factor. This parameter specifies the ratio of the extinction applied to the nebular light to that applied to the starlight, i.e., it gives :math:`f_{\mathrm{neb,ex}} = A_{V,\mathrm{neb}} / A_{V,*}`, as defined in :ref:`ssec-extinction`. As with ``A_V``, this parameter can be set either to a real number, in which case this ratio is treated as constant and equal to the input number, or to the name of a PDF file that specified the distribution of this ratio, formatted as described in :ref:`sec-pdfs`. If this keyword is omitted entirely, the nebular and stellar extinctions are set equal to one another.

.. _ssec-nebular-keywords:

Nebular Keywords
----------------
 
* ``compute_nebular`` (default: ``1``): compute the spectrum that results after starlight is processed through the nebula surrounding each cluster or star? Set to 1 for yes, 0 for no.
* ``atomic_data`` (default: ``lib/atomic/``): directory where the atomic data used for nebular emission calculations is located
* ``nebular_no_metals`` (default: 0): if set to 1, metal lines are not used when computing nebular emission
* ``nebular_den`` (default: ``1e2``): hydrogen number density in :math:`\mathrm{cm}^{-3}` to use in nebular emission computations
* ``nebular_temp`` (default: ``-1``): gas kinetic temperature in K to use in nebular emission computations; if set to non-positive value, the temperature will be determined via the lookup table of cloudy runs for fully sampled IMFs
* ``nebular_logU`` (default: ``-3``): log of dimensionless volume-weighted ionization parameter to assume when computing metal line emission and HII region temperatures from the tabulated cloudy data. At present the allowed values are -3, -2.5, and -2.
* ``nebular_phi`` (default: ``0.73``): fraction of ionizing photons absorbed by H atoms rather than being absorbed by dust grains or rescaping; the default value of ``0.73``, taken from `McKee & Williams (1997, ApJ, 476, 144) <http://adsabs.harvard.edu/abs/1997ApJ...476..144M>`_ means that 73% of ionizing photons are absorbed by H


.. _ssec-phot-keywords:

Photometric Filter Keywords
---------------------------

These describe the photometry to be computed. Note that none of these keywords have any effect unless ``out_integrated_phot`` or ``out_cluster_phot`` is set to 1.

* ``phot_bands``: photometric bands for which photometry is to be computed. The values listed here can be comma- or whitespace-separated. For a list of available photometric filters, see the file ``lib/filters/FILTER_LIST``. In addition to these filters, slug always allows four special "bands":
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

.. _ssec-yield-keywords:

Yield Keywords
--------------

These keywords control the calculation of chemical yields. See
:ref:`ssec-yields` for explanations of the physical models
corresponding to these choices.

* ``yield_dir`` (default: ``lib/yields``): directory where the
  stellar yield tables are located. Note that the file name and
  format is hardcoded, so if you want to use a different format,
  you will have to write source code to do so.
* ``yield_mode`` (default: ``sukhbold16+karakas16+doherty14``):
  sources for yields information. Valid options are:
  
  * ``sukhbold16+karakas16+doherty14``: core collapse supernova yields
    from Sukhbold et al. (2016) plus AGB star yields from Karakas &
    Lugaro (2016) plus Doherty et al. (2014)
  * ``karakas16+doherty14``: AGB star yields as in the first option,
    no core collapse supernova yields
  * ``sukhbold16``: core collapse superonva yields as in the first
    option, no AGB star yields
     
* ``no_decay_isotopes`` (default: ``0``): if set to a non-zero value,
  this option disables radioactive decay of unstable isotopes
* ``isotopes_included`` (default: ``intersection``): controls how to
  handle isotopes that are present in some yield tables but not
  others. Valid options are:
  
  * ``intersection``: only include isotopes present in all yield
    tables
  * ``union``: include all isotopes found in any of the yield tables


