.. highlight:: rest

.. _sec-output:

Output Files and Format
=======================

SLUG can produce 7 output files, though the actual number produced depends on the setting for the ``out_*`` keywords in the parameter file. (Additional output files can be produced by :ref:`sec-cloudy-slug`, and are documented in that section rather than here.)

The only file that is always produced is the summary file, which is named ``MODEL_NAME_summary.txt``, where ``MODEL_NAME`` is the value given by the ``model_name`` keyword in the parameter file. This file contains some basic summary information for the run, and is always formatted as ASCII text regardless of the output format requested.

The other eight output files all have names of the form ``MODEL_NAME_xxx.ext``, where the extension ``.ext`` is one of ``.txt``, ``.bin``, or ``.fits`` depending on the ``output_mode`` specified in the parameter file, and ``xxx`` is ``integrated_prop``, ``integrated_spec``, ``integrated_phot``, ``integrated_yield``, ``cluster_prop``, ``cluster_spec``, ``cluster_phot``, or ``cluster_yield``. The production of these output files is controlled by the parameters ``out_integrated``, ``out_integrated_spec``, ``out_integrated_phot``, ``out_integrated_yield``, ``out_cluster``, ``out_cluster_spec``, ``out_cluster_phot``, and ``out_cluster_yield`` in the parameter file. 

The easiest way to read these output files is with :ref:`sec-slugpy`, which can parse them and store the information in python structures. However, for users who wish to write their own parsers or examine the data directly, the format is documented below. The following conventions are used throughout, unless noted otherwise:

* Masses are in :math:`M_\odot`
* Times in year
* Wavelengths are in Angstrom
* Specific luminosities are in erg/s/Angstrom
* For ``binary`` outputs, variable types refer to C++ types


The ``integrated_prop`` File
----------------------------

This file contains data on the bulk physical properties of the galaxy as a whole. It consists of a series of entries containing the following fields:

* ``Trial``: which trial these data are from
* ``Time``: evolution time at which the output is produced
* ``TargetMass``: target mass of stars in the galaxy up that time, if the IMF and SFH were perfectly sampled
* ``ActualMass``: actual mass of stars produced in the galaxy up to that time; generally not exactly equal to ``TargetMass`` due to finite sampling of the IMF and SFH
* ``LiveMass``:  current mass of all stars in the galaxy, accounting for the effects of stellar evolution (mass loss, supernovae); does not include the mass of stellar remnants (black holes, neutron stars, white dwarfs)
* ``StellarMass``: same as ``LiveMass``, but including the mass of stellar remnants; at present, remnant mass calculation is hard-coded to use the initial-final mass relation of `Kruijssen (2009, A&A, 507, 1409) <http://adsabs.harvard.edu/abs/2009A%26A...507.1409K>`_
* ``ClusterMass``: current mass of all stars in the galaxy that are presently in clusters; same as ``LiveMass``, but only including the stars in clusters
* ``NumClusters``: number of non-disrupted clusters present in the galaxy at this time
* ``NumDisClust``: number of disrupted clusters present in the galaxy at this time
* ``NumFldStars``: number of field stars present in the galaxy at this time; this count only includes those stars being treated stochastically (see the parameter ``min_stoch_mass`` in :ref:`ssec-stellar-keywords`)
* ``VPx``: values drawn for variable parameter ``x`` (0,1,2 etc...); present only if SLUG was run with a variable mode IMF

If ``output_mode`` is ``ascii``, these data are output in a series of columns, with different trials separated by lines of dashes. If ``output_mode`` is ``fits``, the data are stored as a FITS binary table extension, with one column for each of the variables above, plus an additional column giving the trial number for that entry. Both the ASCII- and FITS-formatted output should be fairly self-documenting.

For ``binary`` output, the file consists of a series of records containing the following variables

* ``Trial`` (``unsigned int``)
* ``Time`` (``double``)
* ``TargetMass`` (``double``)
* ``ActualMass`` (``double``)
* ``LiveMass`` (``double``)
* ``StellarMass`` (``double``)
* ``ClusterMass`` (``double``)
* ``NumClusters`` (``std::vector<double>::size_type``, usually ``unsigned long long``)
* ``NumDisClust`` (``std::vector<double>::size_type``, usually ``unsigned long long``)
* ``NumFldStars`` (``std::vector<double>::size_type``, usually ``unsigned long long``)
* ``VPx`` (``double``); present only if ``nvps`` is greater than 0, with one entry present for each variable parameter ``x``

The first record in the binary file indicates the number of variable parameters used.

* ``nvps`` (``integer``): the number of variable parameters for the IMF.

There is one record of this form for each output time, with different trials ordered sequentially, so that all the times for one trial are output before the first time for the next trial.

.. _ssec-int-spec-file:

The ``integrated_spec`` File
----------------------------

This file contains data on the spectra of the entire galaxy, and consists of a series of entries containing the following fields:

* ``Trial``: which trial these data are from
* ``Time``: evolution time at which the output is produced
* ``Wavelength``: observed frame wavelength at which the stellar spectrum is evaluated
* ``L_lambda``: specific luminosity at the specified wavelength, before extinction or nebular effects are applied
* ``Wavelength_neb``: observed frame wavelength at which the stellar plus nebular spectrum is evaluated (present only if SLUG was run with nebular emission enabled)
* ``L_lambda_neb``: specific luminosity at the specified wavelength, after the light has been processed through the nebula (only present if SLUG was run with nebular emission enabled)
* ``Wavelength_ex``: observed frame wavelength at which the extincted stellar spectrum is evaluated (present only if SLUG was run with extinction enabled)
* ``L_lambda_ex``: specific luminosity at the specified wavelength after extinction is applied, but without the effects of the nebula (only present if SLUG was run with extinction enabled)
* ``Wavelength_neb_ex``: observed frame wavelength at which the extincted stellar plus nebular spectrum is evaluated (present only if SLUG was run with nebular processing and  extinction enabled)
* ``L_lambda_neb_ex``: specific luminosity at the specified wavelength, after the light is first processed by the nebular and then subjected to dust extinction (only present if SLUG was run with both extinction and nebular emission enabled)

If ``output_mode`` is ``ascii``, these data are output in a series of columns, with different trials separated by lines of dashes. In ``ascii`` mode, only a single ``Wavelength`` column is present, and for those wavelengths that are not included in one of the grids, some entries may be blank.

If ``output_mode`` is ``fits``, the output FITS file has two binary table extensions. The first table contains a field ``Wavelength`` listing the wavelengths at which the stellar spectra are given; if nebular emission was enabled in the SLUG calculation, there is also a field ``Wavelength_neb`` giving the nebular wavelength grid, and if extinction was enabled the table has a field ``Wavelength_ex`` listing the wavelengths at which the extincted spectrum is computed. If both nebular emission and extinction were included, the field ``Wavelength_neb_ex`` exists as well, giving the wavelength grid for that spectrum. The second table has three fields, ``Trial``, ``Time``, and ``L_lambda`` giving the trial number, time, and stellar spectrum. It may also contain fields ``L_lambda_neb``, ``L_lambda_ex``, and ``L_lambda_neb_ex`` giving the stellar plus nebular spectrum, extincted stellar spectrum, and extincted stellar plus nebular spectrum. Both the ASCII- and FITS-formatted output should be fairly self-documenting.

For binary output, the file is formatted as follows. The file starts with

* ``Nebular`` (``byte``): a single byte, with a value of 0 indicating that nebular processing was not enabled for this run, and a value of 1 indicating that it was enabled
* ``Extinct`` (``byte``): a single byte, with a value of 0 indicating that extinction was not enabled for this run, and a value of 1 indicating that it was enabled
* ``NWavelength`` (``std::vector<double>::size_type``, usually ``unsigned long long``): the number of wavelength entries in the stellar spectra
* ``Wavelength`` (``NWavelength`` entries of type ``double``)
* ``NWavelength_neb`` (``std::vector<double>::size_type``, usually ``unsigned long long``): the number of wavelength entries in the stellar plus nebular spectra; only present if ``Nebular`` is 1
* ``Wavelength_neb`` (``NWavelength_neb`` entries of type ``double``)
* ``NWavelength_ex`` (``std::vector<double>::size_type``, usually ``unsigned long long``): the number of wavelength entries in the extincted spectra; only present if ``Extinct`` is 1
* ``Wavelength_ex`` (``NWavelength_ex`` entries of type ``double``); only present if ``Extinct`` is 1
* ``NWavelength_neb_ex`` (``std::vector<double>::size_type``, usually ``unsigned long long``): the number of wavelength entries in the extincted nebular plus stellar spectra; only present if ``Nebular`` and ``Extinct`` are both 1
* ``Wavelength_ex`` (``NWavelength_neb_ex`` entries of type ``double``); only present if ``Nebular`` and ``Extinct`` are both 1

and then contains a series of records in the format

* ``Trial`` (``unsigned int``)
* ``Time`` (``double``)
* ``L_lambda`` (``NWavelength`` entries of type ``double``)
* ``L_lambda_neb`` (``NWavelength_neb`` entries of type ``double``); only present if ``Nebular`` is 1
* ``L_lambda_ex`` (``NWavelength_ex`` entries of type ``double``); only present if ``Extinct`` is 1
* ``L_lambda_neb_ex`` (``NWavelength_neb_ex`` entries of type ``double``); only present if ``Nebular`` and ``Extinct`` are both 1

There is one such record for each output time, with different trials ordered sequentially, so that all the times for one trial are output before the first time for the next trial.

.. _ssec-int-phot-file:

The ``integrated_phot`` File
----------------------------

This file contains data on the photometric properties of the entire galaxy, and consists of a series of entries containing the following fields:

* ``Trial``: which trial these data are from
* ``Time``: evolution time at which the output is produced
* ``PhotFilter1``: photometric value through filter 1, where filters follow the order in which they are specified by the ``phot_bands`` keyword; units depend on the value of ``phot_mode`` (see :ref:`ssec-phot-keywords`)
* ``PhotFilter2``
* ``PhotFilter3``
* ``...``
* ``PhotFilter1_neb``: photometric value through filter 1 for the spectrum after nebular processing, in the same units as ``PhotFilter1``; only present if SLUG was run with nebular processing enabled
* ``PhotFilter2_neb``
* ``PhotFilter3_neb``
* ``...``
* ``PhotFilter1_ex``: photometric value through filter 1 for the extincted spectrum, in the same units as ``PhotFilter1``; only present if SLUG was run with extinction enabled
* ``PhotFilter2_ex``
* ``PhotFilter3_ex``
* ``...``
* ``PhotFilter1_neb_ex``: photometric value through filter 1 for the spectrum after nebular processing and extinction, in the same units as ``PhotFilter1``; only present if SLUG was run with both nebular processing and extinction enabled
* ``PhotFilter2_neb_ex``
* ``PhotFilter3_neb_ex``
* ``...``

If ``output_mode`` is ``ascii``, these data are output in a series of
columns, with different trials separated by lines of dashes. The
columns for photometry of the extincted spectrum are present only if
extinction was enabled when SLUG was run. Entries for some filters may
be left blank. If so, this indicates that the photon response function
provided for that filter extends beyond the wavelength range covered
by the provided extinction curve. Since the extincted spectrum cannot
be computed over the full range of the filter in this case, photometry
for that filter cannot be computed either.

If ``output_mode`` is ``fits``, the data are stored as a series of
columns in a binary table extension to the FITS file; the filter names
and units are included in the header information for the columns. If
SLUG was run with nebular emission enabled, for each filter ``FILTERNAME``
there is a corresponding column ``FILTERNAME_neb`` giving the photometric
value for the nebular-processed spectrum. Similarly, the columns
``FILTERNAME_ex`` and ``FILTERNAME_neb_ex`` give the photometric values
for the extincted stellar and stellar + nebular spectra, respectively.
Some of the extincted values may be ``NaN``; this
indicates that the photon response function provided for that filter
extends beyond the wavelength range covered by the provided extinction
curve. In addition to the time and photometric filter values, the FITS
file contains a column specifying the trial number for that
entry. Both the ASCII- and FITS-formatted output should be fairly
self-documenting.
 
For binary output, the file is formatted as follows. The file starts with

* ``NFilter`` (stored as ``ASCII text``): number of filters used
* ``FilterName`` ``FilterUnit`` (``NFilter`` entries stored as ``ASCII
  text``): the name and units for each filter are listed in ASCII, one
  filter-unit pair per line
* ``Nebular`` (``byte``): a single byte, with a value of 0 indicating
  that nebular processing was not enabled for this run, and a value of 1
  indicating that it was enabled
* ``Extinct`` (``byte``): a single byte, with a value of 0 indicating
  that extinction was not enabled for this run, and a value of 1
  indicating that it was enabled

This is followed by a series of entries of the form

* ``Trial`` (``unsigned int``)
* ``Time`` (``double``)
* ``PhotFilter`` (``NFilter`` entries of type ``double``)
* ``PhotFilter_neb`` (``NFilter`` entries of type ``double``); only present if ``Nebular`` is 1.
* ``PhotFilter_ex`` (``NFilter`` entries of type ``double``); only present if ``Extinct`` is 1. Note that some values may be ``NaN`` if photometry could not be computed for that filter (see above).
* ``PhotFilter_neb_ex`` (``NFilter`` entries of type ``double``); only present if ``Nebular`` and ``Extinct`` are both 1. Note that some values may be ``NaN`` if photometry could not be computed for that filter (see above).

There is one such record for each output time, with different trials ordered sequentially, so that all the times for one trial are output before the first time for the next trial.

.. _ssec-int-yield-file:

The ``integrated_yield`` File
-----------------------------

This file contains data on the integrated chemical yield of the entire
galaxy, and consists of a series of entries containing the following
fields:

* ``Trial``: which trial these data are from
* ``Time``: evolution time at which the output is produced
* ``Name``: name (i.e., atomic symbol) of an isotope being produced
* ``Z``: atomic number of an isotope being produced
* ``A``: mass number of an isotope being produced
* ``Yield``: mass of a particular isotope produced up to the specified time; for unstable isotopes, this includes the effects of radioactive decay, so yield can decrease with time under some circumstances

If ``output_mode`` is ``ascii``, these data are output in a series of
columns, with different trials separated by lines of dashes. If
``output_mode`` is ``fits``, the data are stored as a series of
columns in a binary table extension to the FITS file; the ``Name``,
``Z``, and ``A`` fields are placed in the first binary table
extension, and are the same for every output. The ``Time`` and
``Yields`` fields are in the second binary table extension. In
addition to these two fields, the second binary table contains a
column specifying the trial number for each entry.

For binary output, the file is formatted as follows. It starts with

* ``NIso`` (``std::vector<double>::size_type``, usually ``unsigned
  long long``): number of isotopes in the output

This is followed by ``NIso`` entries of the form

* ``Name`` (``char[4]``): isotope name (i.e., element symbol)
* ``Z`` (``unsigned int``)
* ``A`` (``unsigned int``)

The remainder of the file contains records of the from

* ``Trial`` (``unsigned int``)
* ``Time`` (``double``)
* ``Yield`` (``double[NIso]``)

There is one such record for each output time, with different trials
ordered sequentially, so that all the times for one trial are output
before the first time for the next trial.


The ``cluster_prop`` File
-------------------------

This file contains data on the bulk physical properties of the non-disrupted star clusters in the galaxy, with one entry per cluster per time at which that cluster exists. Each entry contains the following fields

* ``UniqueID``: a unique identifier number for each cluster that is preserved across times and output files
* ``Time``: evolution time at which the output is produced
* ``FormTime``: time at which that cluster formed
* ``Lifetime``: amount of time from birth to when the cluster will disrupt
* ``TargetMass``: target mass of stars in the cluster, if the IMF were perfectly sampled
* ``BirthMass``: actual mass of stars present in the cluster at formation
* ``LiveMass``: current mass of all stars in the cluster, accounting for the effects of stellar evolution (mass loss, supernovae); does not include the mass of stellar remnants (black holes, neutron stars, white dwarfs)
* ``StellarMass``: same as ``LiveMass``, but including the mass of stellar remnants; at present, remnant mass calculation is hard-coded to use the initial-final mass relation of `Kruijssen (2009, A&A, 507, 1409) <http://adsabs.harvard.edu/abs/2009A%26A...507.1409K>`_
* ``NumStar``: number of living stars in the cluster at this time; this count only includes those stars being treated stochastically (see the parameter ``min_stoch_mass`` in :ref:`ssec-stellar-keywords`)
* ``MaxStarMass``: mass of most massive star still living in the cluster; this only includes those stars being treated stochastically (see the parameter ``min_stoch_mass`` in :ref:`ssec-stellar-keywords`)
* ``A_V``: visual extinction for that cluster, in mag; present only if SLUG was run with extinction enabled
* ``VPx``: values drawn for variable parameter ``x`` (0,1,2 etc...); present only if SLUG was run with a variable mode IMF

If ``output_mode`` is ``ascii``, these data are output in a series of columns, with different trials separated by lines of dashes. If ``output_mode`` is ``fits``, the data are stored as a FITS binary table extension, with one column for each of the variables above, plus an additional column giving the trial number for that entry. Both the ASCII- and FITS-formatted output should be fairly self-documenting.

For ``binary`` output, the first entry in the file is a header containing

* ``Extinct`` (``byte``): a single byte, with a value of 0 indicating that extinction was not enabled for this run, and a value of 1 indicating that it was enabled
* ``nvps`` (``integer``): the number of variable parameters for the IMF.


Thereafter, the file consists of a series of records, one for each output time, with different trials ordered sequentially, so that all the times for one trial are output before the first time for the next trial. Each record consists of a header containing

* ``Time`` (``double``)
* ``NCluster`` (``std::vector<double>::size_type``, usually ``unsigned long long``): number of non-disrupted clusters present at this time

This is followed by ``NCluster`` entries of the following form:

* ``UniqueID`` (``unsigned long``)
* ``FormationTime`` (``double``)
* ``Lifetime`` (``double``)
* ``TargetMass`` (``double``)
* ``BirthMass`` (``double``)
* ``LiveMass`` (``double``)
* ``StellarMass`` (``double``)
* ``NumStar`` (``std::vector<double>::size_type``, usually ``unsigned long long``)
* ``MaxStarMass`` (``double``)
* ``A_V`` (``double``); present only if ``Extinct`` is 1
* ``VPx`` (``double``); present only if ``nvps`` is greater than 0, with one entry present for each variable parameter ``x``

The ``cluster_spec`` File
-------------------------

This file contains the spectra of the individual clusters, and each entry contains the following fields:

* ``UniqueID``: a unique identifier number for each cluster that is preserved across times and output files
* ``Time``: evolution time at which the output is produced
* ``Wavelength``: observed frame wavelength at which the stellar spectrum is evaluated
* ``L_lambda``: specific luminosity at the specified wavelength, before extinction or nebular effects are applied
* ``Wavelength_neb``: observed frame wavelength at which the stellar plus nebular spectrum is evaluated (present only if SLUG was run with nebular emission enabled)
* ``L_lambda_neb``: specific luminosity at the specified wavelength, after the light has been processed through the nebula (only present if SLUG was run with nebular emission enabled)
* ``Wavelength_ex``: observed frame wavelength at which the extincted stellar spectrum is evaluated (present only if SLUG was run with extinction enabled)
* ``L_lambda_ex``: specific luminosity at the specified wavelength after extinction is applied, but without the effects of the nebula (only present if SLUG was run with extinction enabled)
* ``Wavelength_neb_ex``: observed frame wavelength at which the extincted stellar plus nebular spectrum is evaluated (present only if SLUG was run with nebular processing and  extinction enabled)
* ``L_lambda_neb_ex``: specific luminosity at the specified wavelength, after the light is first processed by the nebular and then subjected to dust extinction (only present if SLUG was run with both extinction and nebular emission enabled)

If ``output_mode`` is ``ascii``, these data are output in a series of columns, with different trials separated by lines of dashes. The columns ``L_lambda_neb``, ``L_lambda_ex``, and ``L_lambda_neb_ex`` are present only if SLUG was run with the appropriate options enabled. Some entries in these fields may be empty; see :ref:`ssec-int-spec-file`.

If ``output_mode`` is ``fits``, the output FITS file has two binary table extensions. The first table contains a field listing the wavelengths at which the spectra are given, in the same format as for :ref:`ssec-int-spec-file`. The second table has always contains the fields ``UniqueId``, ``Time``, ``Trial``, and ``L_lambda`` giving the cluster unique ID, time, trial number, and stellar spectrum. Depending on whether nebular processing and/or extinction were enabled when SLUG was run, it may also contain the fields ``L_lambda_neb``, ``L_lambda_ex``, and ``L_lambda_neb_ex`` giving the nebular-processed, extincted, and nebular-processed plus extincted spectra. Both the ASCII- and FITS-formatted output should be fairly self-documenting.

Output in ``binary`` mode is formatted as follows.  The file starts with

* ``Nebular`` (``byte``): a single byte, with a value of 0 indicating that nebular processing was not enabled for this run, and a value of 1 indicating that it was enabled
* ``Extinct`` (``byte``): a single byte, with a value of 0 indicating that extinction was not enabled for this run, and a value of 1 indicating that it was enabled
* ``NWavelength`` (``std::vector<double>::size_type``, usually ``unsigned long long``): the number of wavelength entries in the stellar spectra
* ``Wavelength`` (``NWavelength`` entries of type ``double``)
* ``NWavelength_neb`` (``std::vector<double>::size_type``, usually ``unsigned long long``): the number of wavelength entries in the stellar plus nebular spectra; only present if ``Nebular`` is 1
* ``Wavelength_neb`` (``NWavelength_neb`` entries of type ``double``)
* ``NWavelength_ex`` (``std::vector<double>::size_type``, usually ``unsigned long long``): the number of wavelength entries in the extincted spectra; only present if ``Extinct`` is 1
* ``Wavelength_ex`` (``NWavelength_ex`` entries of type ``double``); only present if ``Extinct`` is 1
* ``NWavelength_neb_ex`` (``std::vector<double>::size_type``, usually ``unsigned long long``): the number of wavelength entries in the extincted nebular plus stellar spectra; only present if ``Nebular`` and ``Extinct`` are both 1
* ``Wavelength_ex`` (``NWavelength_neb_ex`` entries of type ``double``); only present if ``Nebular`` and ``Extinct`` are both 1

and then contains a series of records, one for each output time, with different trials ordered sequentially, so that all the times for one trial are output before the first time for the next trial. Each record consists of a header containing

* ``Time`` (``double``)
* ``NCluster`` (``std::vector<double>::size_type``, usually ``unsigned long long``): number of non-disrupted clusters present at this time

This is followed by ``NCluster`` entries of the following form:

* ``UniqueID`` (``unsigned long``)
* ``L_lambda`` (``NWavelength`` entries of type ``double``)
* ``L_lambda_neb`` (``NWavelength_neb`` entries of type ``double``); only present if ``Nebular`` is 1
* ``L_lambda_ex`` (``NWavelength_ex`` entries of type ``double``); only present if ``Extinct`` is 1
* ``L_lambda_neb_ex`` (``NWavelength_neb_ex`` entries of type ``double``); only present if ``Nebular`` and ``Extinct`` are both 1


.. _ssec-cluster-phot-file:

The ``cluster_phot`` File
-------------------------

This file contains the photometric values for the individual clusters. Each entry contains the following fields:

* ``UniqueID``: a unique identifier number for each cluster that is preserved across times and output files
* ``Time``: evolution time at which the output is produced
* ``PhotFilter1``: photometric value through filter 1, where filters follow the order in which they are specified by the ``phot_bands`` keyword; units depend on the value of ``phot_mode`` (see :ref:`ssec-phot-keywords`)
* ``PhotFilter2``
* ``PhotFilter3``
* ``...``
* ``PhotFilter1_neb``: photometric value through filter 1 for the spectrum after nebular processing, in the same units as ``PhotFilter1``; only present if SLUG was run with nebular processing enabled
* ``PhotFilter2_neb``
* ``PhotFilter3_neb``
* ``...``
* ``PhotFilter1_ex``: photometric value through filter 1 for the extincted spectrum, in the same units as ``PhotFilter1``; only present if SLUG was run with extinction enabled
* ``PhotFilter2_ex``
* ``PhotFilter3_ex``
* ``...``
* ``PhotFilter1_neb_ex``: photometric value through filter 1 for the spectrum after nebular processing and extinction, in the same units as ``PhotFilter1``; only present if SLUG was run with both nebular processing and extinction enabled
* ``PhotFilter2_neb_ex``
* ``PhotFilter3_neb_ex``
* ``...``

If ``output_mode`` is ``ascii``, these data are output in a series of columns, with different trials separated by lines of dashes. Some of the extincted photometry columns may be blank; see :ref:`ssec-int-phot-file`.

If ``output_mode`` is ``fits``, the data are stored as a series of
columns in a binary table extension to the FITS file; the filter names
and units are included in the header information for the columns. If
SLUG was run with nebular emission enabled, for each filter ``FILTERNAME``
there is a corresponding column ``FILTERNAME_neb`` giving the photometric
value for the nebular-processed spectrum. Similarly, the columns
``FILTERNAME_ex`` and ``FILTERNAME_neb_ex`` give the photometric values
for the extincted stellar and stellar + nebular spectra, respectively.
Some of the extincted values may be ``NaN``; this
indicates that the photon response function provided for that filter
extends beyond the wavelength range covered by the provided extinction
curve. In addition to the time and photometric filter values, the FITS
file contains a column specifying the trial number for that
entry. Both the ASCII- and FITS-formatted output should be fairly
self-documenting.

In ``binary`` output mode, the binary data file starts with

* ``NFilter`` (stored as ``ASCII text``): number of filters used
* ``FilterName`` ``FilterUnit`` (``NFilter`` entries stored as ``ASCII text``): the name and units for each filter are listed in ASCII, one filter-unit pair per line
* ``Nebular`` (``byte``): a single byte, with a value of 0 indicating that nebular processing was not enabled for this run, and a value of 1 indicating that it was enabled
* ``Extinct`` (``byte``): a single byte, with a value of 0 indicating that extinction was not enabled for this run, and a value of 1 indicating that it was enabled

and then contains a series of records, one for each output time , with different trials ordered sequentially, so that all the times for one trial are output before the first time for the next trial. Each record consists of a header containing

* ``Time`` (``double``)
* ``NCluster`` (``std::vector<double>::size_type``, usually ``unsigned long long``): number of non-disrupted clusters present at this time

This is followed by ``NCluster`` entries of the following form:

* ``UniqueID`` (``unsigned long``)
* ``PhotFilter`` (``NFilter`` entries of type ``double``)
* ``PhotFilter_neb`` (``NFilter`` entries of type ``double``); only present if ``Nebular`` is 1.
* ``PhotFilter_ex`` (``NFilter`` entries of type ``double``); only present if ``Extinct`` is 1. Note that some values may be ``NaN`` if photometry could not be computed for that filter (see above).
* ``PhotFilter_neb_ex`` (``NFilter`` entries of type ``double``); only present if ``Nebular`` and ``Extinct`` are both 1. Note that some values may be ``NaN`` if photometry could not be computed for that filter (see above).

.. _ssec-cluster-yield-file:

The ``cluster_yield`` File
--------------------------

This file contains data on the chemical yield of individual star
clusters, and consists of a series of entries containing the following
fields:

* ``UniqueID``: a unique identifier number for each cluster that is
  preserved across times and output files
* ``Time``: evolution time at which the output is produced
* ``Name``: name (i.e., atomic symbol) of an isotope being produced
* ``Z``: atomic number of an isotope being produced
* ``A``: mass number of an isotope being produced
* ``Yield``: mass of a particular isotope produced up to the specified
  time; for unstable isotopes, this includes the effects of
  radioactive decay, so yield can decrease with time under some
  circumstances

If ``output_mode`` is ``ascii``, these data are output in a series of
columns, with different trials separated by lines of dashes. If
``output_mode`` is ``fits``, the data are stored as a series of
columns in a binary table extension to the FITS file; the ``Name``,
``Z``, and ``A`` fields are placed in the first binary table
extension, and are the same for every output. The ``Time`` and
``Yields`` fields are in the second binary table extension. In
addition to these two fields, the second binary table contains a
column specifying the trial number for each entry.

For binary output, the file is formatted as follows. It starts with

* ``NIso`` (``std::vector<double>::size_type``, usually ``unsigned
  long long``): number of isotopes in the output

This is followed by ``NIso`` entries of the form

* ``Name`` (``char[4]``): isotope name (i.e., element symbol)
* ``Z`` (``unsigned int``)
* ``A`` (``unsigned int``)

Thereafter, the file consists of a series of records, one for each output time, with different trials ordered sequentially, so that all the times for one trial are output before the first time for the next trial. Each record consists of a header containing

* ``Time`` (``double``)
* ``NCluster`` (``std::vector<double>::size_type``, usually ``unsigned long long``): number of non-disrupted clusters present at this time

This is followed by ``NCluster`` entries of the following form:

* ``UniqueID`` (``unsigned long``)
* ``Yield`` (``double[NIso]``)

There is one such record for each output time, with different trials
ordered sequentially, so that all the times for one trial are output
before the first time for the next trial.

.. _ssec-checkpoint-files:

Checkpoint Files
----------------

Checkpoint files are identical to regular output files, except that
they begin with a statement of the number of trials they contain. For
ASCII files, this is indicated by a first line of the form ``N_Trials =
N``. For binary files, the file begins with an unsigned integer that
gives the number of trials. For FITS files, the file has a keyword
``N_Trials`` in the first binary table HDU that gives the number of
trials.

The slugpy library can read checkpoint files as well as regular output
files (see :ref:`sec-slugpy`).
