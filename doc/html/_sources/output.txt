.. highlight:: rest

.. _sec-output:

Output Files and Format
=======================

SLUG can produce 7 output files, though the actual number produced depends on the setting for the ``out_*`` keywords in the parameter file. (Additional output files can be produced by :ref:`sec-cloudy-slug`, and are documented in that section rather than here.)

The only file that is always produced is the summary file, which is named ``MODEL_NAME_summary.txt``, where ``MODEL_NAME`` is the value given by the ``model_name`` keyword in the parameter file. This file contains some basic summary information for the run, and is always formatted as ASCII text regardless of the output format requested.

The other six output files all have names of the form ``MODEL_NAME_xxx.ext``, where the extension ``.ext`` is one of ``.txt``, ``.bin``, or ``.fits`` depending on the ``output_mode`` specified in the parameter file, and ``xxx`` is ``integrated_prop``, ``integrated_spec``, ``integrated_phot``, ``cluster_prop``, ``cluster_spec``, or ``cluster_phot``. The production of these output files is controlled by the parameters ``out_integrated``, ``out_integrated_spec``, ``out_integrated_phot``, ``out_cluster``, ``out_cluster_spec``, and ``out_cluster_phot`` in the parameter file. The files are formatted as described below. 

The following conventions are used throughout, unless noted otherwise:

* Masses are in :math:`M_\odot`
* Times in year
* Wavelengths are in Angstrom
* Specific luminosities are in erg/s/Angstrom
* For ``binary`` outputs, variable types refer to C++ types


The ``integrated_prop`` File
----------------------------

This file contains data on the bulk physical properties of the galaxy as a whole. It consists of a series of entries containing the following fields:

* ``Time``: evolution time at which the output is produced
* ``TargetMass``: target mass of stars in the galaxy up that time, if the IMF and SFH were perfectly sampled
* ``ActualMass``: actual mass of stars produced in the galaxy up to that time; generally not exactly equal to ``TargetMass`` due to finite sampling of the IMF and SFH
* ``LiveMass``: actual mass of stars produced in the galaxy up to that time, and which have not yet reached the end of their lives (as marked by the final entry in the stellar evolution tracks)
* ``ClusterMass``: actual mass of stars produced in the galaxy up to that time that are still members of non-disrupted clusters
* ``NumClusters``: number of non-disrupted clusters present in the galaxy at this time
* ``NumDisClust``: number of disrupted clusters present in the galaxy at this time
* ``NumFldStars``: number of field stars present in the galaxy at this time; this count only includes those stars being treated stochastically (see the parameter ``min_stoch_mass`` in :ref:`ssec-phys-keywords`)


If ``output_mode`` is ``ascii``, these data are output in a series of columns, with different trials separated by lines of dashes. If ``output_mode`` is ``fits``, the data are stored as a FITS binary table extension, with one column for each of the variables above, plus an additional column giving the trial number for that entry. Both the ASCII- and FITS-formatted output should be fairly self-documenting.

For ``binary`` output, the file consists of a series of records containing the following variables

* ``Time`` (``double``)
* ``TargetMass`` (``double``)
* ``ActualMass`` (``double``)
* ``LiveMass`` (``double``)
* ``ClusterMass`` (``double``)
* ``NumClusters`` (``std::vector<double>::size_type``, usually ``unsigned long long``)
* ``NumDisClust`` (``std::vector<double>::size_type``, usually ``unsigned long long``)
* ``NumFldStars`` (``std::vector<double>::size_type``, usually ``unsigned long long``)

There is one record of this form for each output time, with different trials ordered sequentially, so that all the times for one trial are output before the first time for the next trial.

.. _ssec-int-spec-file:

The ``integrated_spec`` File
----------------------------

This file contains data on the spectra of the entire galaxy, and consists of a series of entries containing the following fields:

* ``Time``: evolution time at which the output is produced
* ``Wavelength``: observed frame wavelength at which the spectrum is evaluated
* ``L_lambda``: specific luminosity at the specified wavelength, without extinction
* ``L_lambda_ex``: specific luminosity at the specified wavelength after extinction is applied (only present if SLUG was run with extinction enabled)

If ``output_mode`` is ``ascii``, these data are output in a series of columns, with different trials separated by lines of dashes. Note that some entries for ``L_lambda_ex`` may be blank. This indicates wavelengths for which a stellar atmosphere model was available, but outside the wavelength range included in the user-specified extinction curve file. Extincted luminosities are unavailable at these wavelengths due to the lack of an extinction curve.

 If ``output_mode`` is ``fits``, the output FITS file has two binary table extensions. The first table contains a field listing the wavelengths at which the spectra are given; if extinction was enabled in the SLUG calculation, it also contains a field ``Wavelength_ex`` listing the wavelengths at which the extincted spectrum is computed. The second table has three fields if SLUG was run without extinction, or four fields if it was run with extinction. The first three fields, always present, give the trial number, the time, and the spectrum ``L_lambda`` at that time. The final field, ``L_lambda_ex``, gives the extincted spectrum at that time and trial. Both the ASCII- and FITS-formatted output should be fairly self-documenting.

For binary output, the file is formatted as follows. The file starts with

* ``Extinct`` (``byte``): a single byte, with a value of 0 indicating that extinction was not enabled for this run, and a value of 1 indicating that it was enabled
* ``NWavelength`` (``std::vector<double>::size_type``, usually ``unsigned long long``): the number of wavelength entries in the spectra
* ``Wavelength`` (``NWavelength`` entries of type ``double``)
* ``NWavelength_ex`` (``std::vector<double>::size_type``, usually ``unsigned long long``): the number of wavelength entries in the extinted spectra; only present if ``Extinct`` is 1
* ``Wavelength_ex`` (``NWavelength_ex`` entries of type ``double``); only present if ``Extinct`` is 1

and then contains a series of records in the format

* ``Time`` (``double``)
* ``L_lambda`` (``NWavelength`` entries of type ``double``)
* ``L_lambda_ex`` (``NWavelength_ex`` entries of type ``double``); only present if ``Extinct`` is 1

There is one such record for each output time, with different trials ordered sequentially, so that all the times for one trial are output before the first time for the next trial.

.. _ssec-int-phot-file:

The ``integrated_phot`` File
----------------------------

This file contains data on the photometric properties of the entire galaxy, and consists of a series of entries containing the following fields:

* ``Time``: evolution time at which the output is produced
* ``PhotFilter1``: photometric value through filter 1, where filters follow the order in which they are specified by the ``phot_bands`` keyword; units depend on the value of ``phot_mode`` (see :ref:`ssec-phot-keywords`)
* ``PhotFilter2``
* ``PhotFilter3``
* ``...``
* ``PhotFilter1_ex``: photometric value through filter 1 for the extincted spectrum, in the same units as ``PhotFilter1``; only present if SLUG was run with extinction enabled
* ``PhotFilter2_ex``
* ``PhotFilter3_ex``
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
SLUG was run with extinction enabled, in for each filter
``FILTERNAME`` there is a corresponding column ``FILTERNAME_ex``
containing the photometric value for that filter applied to the
extincted spectrum. Some of these values may be ``NaN``; this
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
* ``Extinct`` (``byte``): a single byte, with a value of 0 indicating
  that extinction was not enabled for this run, and a value of 1
  indicating that it was enabled

This is followed by a series of entries of the form

* ``Time`` (``double``)
* ``PhotFilter`` (``NFilter`` entries of type ``double``)
* ``PhotFilter_ex`` (``NFilter`` entries of type ``double``); only present if ``Extinct`` is 1. Note that some values may be ``NaN`` if photometry could not be computed for that filter (see above).

There is one such record for each output time, with different trials ordered sequentially, so that all the times for one trial are output before the first time for the next trial.

The ``cluster_prop`` File
-------------------------

This file contains data on the bulk physical properties of the non-disrupted star clusters in the galaxy, with one entry per cluster per time at which that cluster exists. Each entry contains the following fields

* ``UniqueID``: a unique identifier number for each cluster that is preserved across times and output files
* ``Time``: evolution time at which the output is produced
* ``FormTime``: time at which that cluster formed
* ``Lifetime``: amount of time from birth to when the cluster will disrupt
* ``TargetMass``: target mass of stars in the cluster, if the IMF were perfectly sampled
* ``BirthMass``: actual mass of stars present in the cluster at formation
* ``LiveMass``: actual mass of stars produced in the cluster at this output time that have not yet reached the end of their lives (as marked by the final entry in the stellar evolution tracks)
* ``NumStar``: number of living stars in the cluster at this time; this count only includes those stars being treated stochastically (see the parameter ``min_stoch_mass`` in :ref:`ssec-phys-keywords`)
* ``MaxStarMass``: mass of most massive star still living in the cluster; this only includes those stars being treated stochastically (see the parameter ``min_stoch_mass`` in :ref:`ssec-phys-keywords`)
* ``A_V``: visual extinction for that cluster, in mag; present only if SLUG was run with extinction enabled

If ``output_mode`` is ``ascii``, these data are output in a series of columns, with different trials separated by lines of dashes. If ``output_mode`` is ``fits``, the data are stored as a FITS binary table extension, with one column for each of the variables above, plus an additional column giving the trial number for that entry. Both the ASCII- and FITS-formatted output should be fairly self-documenting.

For ``binary`` output, the first entry in the file is a header containing

* ``Extinct`` (``byte``): a single byte, with a value of 0 indicating that extinction was not enabled for this run, and a value of 1 indicating that it was enabled

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
* ``NumStar`` (``std::vector<double>::size_type``, usually ``unsigned long long``)
* ``MaxStarMass`` (``double``)
* ``A_V`` (``double``); present only if ``Extinct`` is 1


The ``cluster_spec`` File
-------------------------

This file contains the spectra of the individual clusters, and each entry contains the following fields:

* ``UniqueID``: a unique identifier number for each cluster that is preserved across times and output files
* ``Time``: evolution time at which the output is produced
* ``Wavelength``: observed frame wavelength at which the spectrum is evaluated
* ``L_lambda``: specific luminosity at the specified wavelength
* ``L_lambda_ex``: specific luminosity at the specified wavelength after extinction is applied (only present if SLUG was run with extinction enabled)

If ``output_mode`` is ``ascii``, these data are output in a series of columns, with different trials separated by lines of dashes. The column ``L_lambda_ex`` is present only if SLUG was run with extinction enabled. Some entries for ``L_lambda_ex`` may be empty; see :ref:`ssec-int-spec-file`.

If ``output_mode`` is ``fits``, the output FITS file has two binary table extensions. The first table contains a field listing the wavelengths at which the spectra are given; if extinction was enabled in the SLUG calculation, it also contains a field ``Wavelength_ex`` listing the wavelengths at which the extincted spectrum is computed. The second table has four fields if SLUG was run without extinction, or five fields if it was run with extinction. The first four fields, always present, give the trial number, unique ID of the cluster, the time, and the spectrum ``L_lambda`` at that time. The final field, ``L_lambda_ex``, gives the extincted spectrum at that time and trial. Both the ASCII- and FITS-formatted output should be fairly self-documenting.

Output in ``binary`` mode is formatted as follows.  The file starts with

* ``Extinct`` (``byte``): a single byte, with a value of 0 indicating that extinction was not enabled for this run, and a value of 1 indicating that it was enabled
* ``NWavelength`` (``std::vector<double>::size_type``, usually ``unsigned long long``): the number of wavelength entries in the spectra
* ``Wavelength`` (``NWavelength`` entries of type ``double``)
* ``NWavelength_ex`` (``std::vector<double>::size_type``, usually ``unsigned long long``): the number of wavelength entries in the extinted spectra; only present if ``Extinct`` is 1
* ``Wavelength_ex`` (``NWavelength_ex`` entries of type ``double``); only present if ``Extinct`` is 1

and then contains a series of records, one for each output time, with different trials ordered sequentially, so that all the times for one trial are output before the first time for the next trial. Each record consists of a header containing

* ``Time`` (``double``)
* ``NCluster`` (``std::vector<double>::size_type``, usually ``unsigned long long``): number of non-disrupted clusters present at this time

This is followed by ``NCluster`` entries of the following form:

* ``UniqueID`` (``unsigned long``)
* ``L_lambda`` (``NWavelength`` entries of type ``double``)
* ``L_lambda_ex`` (``NWavelength_ex`` entries of type ``double``); only present if ``Extinct`` is 1


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
* ``PhotFilter1_ex``: photometric value through filter 1 for the extincted spectrum, in the same units as ``PhotFilter1``; only present if SLUG was run with extinction enabled
* ``PhotFilter2_ex``
* ``PhotFilter3_ex``
* ``...``

If ``output_mode`` is ``ascii``, these data are output in a series of columns, with different trials separated by lines of dashes. Some of the extincted photometry columns may be blank; see :ref:`ssec-int-phot-file`.

If ``output_mode`` is ``fits``, the data are stored as a series of columns in a binary table extension to the FITS file; the filter names and units are included in the header information for the columns. If SLUG was run with extinction enabled, in for each filter ``FILTERNAME`` there is a corresponding column ``FILTERNAME_ex`` containing the photometric value for that filter applied to the extincted spectrum. Some of these values may be ``NaN``; see :ref:`ssec-int-phot-file`. In addition to the time, unique ID, and photometric filter values, the FITS file contains a column specifying the trial number for that entry. Both the ASCII- and FITS-formatted output should be fairly self-documenting.

In ``binary`` output mode, the binary data file starts with

* ``NFilter`` (stored as ``ASCII text``): number of filters used
* ``FilterName`` ``FilterUnit`` (``NFilter`` entries stored as ``ASCII text``): the name and units for each filter are listed in ASCII, one filter-unit pair per line
* ``Extinct`` (``byte``): a single byte, with a value of 0 indicating that extinction was not enabled for this run, and a value of 1 indicating that it was enabled

and then contains a series of records, one for each output time , with different trials ordered sequentially, so that all the times for one trial are output before the first time for the next trial. Each record consists of a header containing

* ``Time`` (``double``)
* ``NCluster`` (``std::vector<double>::size_type``, usually ``unsigned long long``): number of non-disrupted clusters present at this time

This is followed by ``NCluster`` entries of the following form:

* ``UniqueID`` (``unsigned long``)
* ``PhotFilter`` (``NFilter`` entries of type ``double``)
* ``PhotFilter_ex`` (``NFilter`` entries of type ``double``); only present if ``Extinct`` is 1. Note that some values may be ``NaN`` if photometry could not be computed for that filter (see above).



