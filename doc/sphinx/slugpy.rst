.. highlight:: rest

.. _sec-slugpy:

slugpy -- The Python Helper Library
===================================

Installing slugpy
-----------------

SLUG comes with the python module slugpy, which contains an extensive set of routines for reading, writing, and manipulating SLUG outputs. You can install slugpy one of two ways.

1. Using ``Make``. If you compile the main slug code by doing ``Make`` in the main slug directory, the c slugpy extensions will be build automatically. Once that is done, you will be able to use slugpy just by importing it, provided that the slugpy directory is in your python import path.

2. Using ``setup.py``. The slug distribution comes with a ``setup.py`` script that follows the standard python package convensions. Just do::

     python setup.py build

to build the c extensions in place, which will let you import slugpy from the directory where it is located. Alternately, do::

  python setup.py install

to install as a site package. Installing as a site package often requires root permissions, in which case you should do::

  sudo python setup.py install

or::

  python setup.py install --user

instead.

A note on compiling slugpy with setup.py: the slugpy c extensions require the `GNU Scientific Library <http://www.gnu.org/software/gsl/>`_; to build or install slugpy, the appropriate headers must be in your default include path, and the appropriate libraries in your default link path. If they are not, you can tell setup where they are located by creating a file named ``setup.cfg`` in the slug2 directory, which contains the lines::

  [build_ext]
  include_dirs=/PATH/TO/GSL/HEADER
  library_dirs=/PATH/TO/GSL/LIBRARIES

Then you should be able to build and install slugpy with ``setup.py``.

Basic Usage
-----------

The most common task is to read a set of SLUG outputs into memory so that they can be processed. To read the data from a SLUG run using slugpy, one can simply do the following::

   from slugpy import *
   idata = read_integrated('SLUG_MODEL_NAME')
   cdata = read_cluster('SLUG_MODEL_NAME')

The ``read_integrated`` function reads all the integrated-light data (i.e., the data stored in the ``_integrated_*`` files -- see :ref:`sec-output`) for a SLUG output whose name is given as the argument. This is the base name specified by the ``model_name`` keyword (see :ref:`ssec-basic-keywords`), without any extensions; the slugpy library will automatically determine which outputs are available and in what format, and read the appropriate files. It returns a ``namedtuple`` containing all the output data available for that simulation. Note that some of these fields will only be present if the cloudy-slug interface (see :ref:`sec-cloudy-slug`) was used to process the SLUG output through cloudy to predict nebular emission, and some will be present only if extinction was enabled when SLUG was run. The fields returned are as follows:

* time: output times
* target_mass: target stellar mass at each time
* actual_mass: actual stellar mass at each time
* live_mass: mass of currently-alive stars
* cluster_mass: mass of living stars in non-disrupted clusters
* num_clusters: number of non-disrupted clusters
* num_dis_clusters: number of disrupted clusters
* num_fld_stars: number of still-living stars that formed in the field
* wl: wavelengths of output stellar spectra (in Angstrom)
* spec: integrated spectrum of all stars, expressed as a specific luminosity (erg/s/Angstrom)
* filter_names: list of photometric filter names
* filter_units: list of units for photometric outputs
* filter_wl_eff: effective wavelength for each photometric filter
* filter_wl: list of wavelengths for each filter at which the response function is given (in Angstrom)
* filter_response: photon response function for each filter at each wavelength (dimensionless)
* filter_beta: index :math:`\beta` used to set the normalization for each filter -- see :ref:`ssec-spec-phot`
* filter_wl_c: pivot wavelength used to set the normalization for each filter for which :math:`\beta \neq 0` -- see :ref:`ssec-spec-phot`
* phot: photometry of the stars in each filter
* isotope_name: element symbols for the isotopes whose yields are reported
* isotope_Z: atomic numbers for the isotopes whose yields are reported
* isotope_A: atomic numbers for the isotopes whose yields are reported
* yld: yield of each isotope at each time

The following fields are present only if SLUG was run with nebular processing enabled:

* wl_neb: same as wl, but for the spectrum that emerges after the starlight has passed through the nebulae around the emitting clusters and field stars. The nebular grid is finer than the stellar grid, because it contains wavelength extra entries around prominent lines so that the lines are resolved on the grid
* spec_neb: same as spec, but for the nebular-processed spectrum
* phot_neb: same as phot, but for the nebular-processed spectrum

The following fields are present only if SLUG was run with extinction enabled:

* wl_ex: wavelengths of output stellar spectra after extinction has been applied(in Angstrom). Note that wl_ex will generally cover a smaller wavelength range than wl, because the extinction curve used may not cover the full wavelength range of the stellar spectra. Extincted spectra are computed only over the range covered by the extinction curve.
* spec_ex: same as spec, but for the extincted spectrum
* phot_ex: same as phot, but for the extincted spectrum. Note that some values may be ``NaN``. This indicates that photometry of the extincted spectrum could not be computed for that filter, because the filter response curve extends to wavelengths outside the range covered by the extinction curve.

The following fields are present only if SLUG was run with both nebular processing and extinction enabled:

* wl_neb_ex: same as wl_neb, but for the extincted, nebular-processed spectrum. Will be limited to the same wavelength range as wl_ex.
* spec_neb_ex: same as spec_neb, but with extinction applied
* phot_neb_ex: same as phot_neb, but wtih extinction applied. Note that some values may be ``NaN``. This indicates that photometry of the extincted spectrum could not be computed for that filter, because the filter response curve extends to wavelengths outside the range covered by the extinction curve.

The following fields are present only for runs that have been processed through the cloudy_slug interface (see :ref:`sec-cloudy-slug`):

* cloudy_wl: wavelengths of the output nebular spectra (in Angstrom)
* cloudy_inc: incident stellar radiation field, expressed as a specific luminosity (erg/s/Angstrom) -- should be the same as spec, but binned onto cloudy's wavelength grid; provided mainly as a bug-checking diagnostic
* cloudy_trans: the transmitted stellar radiation field computed by cloudy, expressed as a specific luminosity (erg/s/Angstrom) -- this is the radiation field of the stars after it has passed through the HII region, and is what one would see in an observational aperture centered on the stars with negligible contribution from the nebula
* cloudy_emit: the emitted nebular radiation field computed by cloudy, expressed as a specific luminosity (erg/s/Angstrom) -- this is the radiation emitted by the nebula excluding the stars, and is what one would see in an observational aperture that included the nebula but masked out the stars
* cloudy_trans_emit: the sum of the transmitted stellar and emitted nebular radiation, expressed as a specific luminosity (erg/s/Angstrom) -- this is what one would see in an observational aperture covering the both the stars and the nebula
* cloudy_linelabel: list of emitting species for the line luminosities computed by cloudy, following cloudy's 4-letter notation
* cloudy_linewl: wavelengths of all the lines computed by cloudy (in Angstrom)
* cloudy_linelum: luminosities of the lines computed by cloudy (in erg/s)
* cloudy_filter_names, cloudy_filter_units, cloudy_filter_wl_eff, cloudy_filter_wl, cloudy_filter_response, cloudy_filter_beta, cloudy_filter_wl_c: exactly the same as the corresponding fields without the cloudy prefix, but for the photometric filters applied to the cloudy output
* cloudy_phot_trans, cloudy_phot_emit, and cloudy_phot_trans_emit: photometry of the transmitted, emitted, and transmitted+emitted radiation field provided by cloudy_trans, cloudy_emit, and cloudy_trans_emit

For the above fields, quantities that are different for each trial and each time are stored as numpy arrays with a shape (N_times, N_trials) for scalar quantities (e.g., actual_mass), or a shape (N, N_times, N_trials) for quantities that are vectors of length N (e.g., the spectrum).

The ``read_cluster`` function is analogous, except that instead of reading the whole-galaxy data, it reads data on the individual star clusters, as stored in the ``_cluster_*`` output files. It returns the following fields:

* id: a unique identifier number for each cluster; this is guaranteed to be unique across both times and trials, so that if two clusters in the list have the same id number, that means that the data given are for the same cluster at two different times in its evolution
* trial: the trial number in which that cluster appeared
* time: the time at which the data for that cluster are computed
* form_time: the time at which that cluster formed
* lifetime: the between when the cluster formed and when it will disrupt
* target_mass: the target stellar mass of the cluster
* actual_mass: the actual stellar mass of the cluter
* live_mass: the mass of all still-living stars in the cluster
* num_star: the number of stars in the cluster
* max_star_mass: the mass of the single most massive still-living star in the cluster
* A_V: the visual extinction for this cluster, in mag; present only if SLUG was run with extinction enabled
* All the remaining fields are identical to those listed above for integrated quantities, starting with wl

For all these fields, scalar quantities that are different for each cluster (e.g., actual_mass) will be stored as arrays of shape (N_cluster); vector quantities that are different for each cluster (e.g., spec) will be stored as arrays of shape (N_cluster, N).

The following fields are present only if SLUG was run with a Variable Mode IMF:

* VPx: The value drawn for variable parameter x (0,1,2...) in each trial. The parameters are numbered in the order they are defined in the IMF definition file.

These fields are present in both the cluster and integrated outputs if a simulation has been run using the variable mode IMF.

Full Documentation of slugpy
----------------------------

.. automodule:: slugpy
   :members:
