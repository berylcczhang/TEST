"""
This defines a class that can be used to estimate the PDF of star
cluster properties (mass, age, extinction) from a set of input
photometry in various bands.
"""

import numpy as np
import copy
import os
import os.path as osp
import warnings
from copy import deepcopy
import errno
try:
    # Python 3
    from urllib.request import urlopen
except ImportError:
    # Python 2
    from urllib2 import urlopen

# Import the data reading and Bayesian inference stuff we need
from ..bayesphot import bp
from ..read_cluster_prop import read_cluster_prop
from ..read_cluster_phot import read_cluster_phot

##################################################################
# Sample density of the default library                          #
##################################################################
def _default_sample_density(physprop):
    logm = physprop[:,0]
    logt = physprop[:,1]
    sden = np.ones(len(logm))
    sden[logm > 4] = sden[logm > 4] * 1.0/10.**(logm[logm > 4]-4)
    sden[logt > 8] = sden[logt > 8] * 1.0/10.**(logt[logt > 8]-8)
    return sden

##################################################################
# Define the cluster_slug class                                  #
##################################################################

class cluster_slug(object):
    """
    A class that can be used to estimate the PDF of star cluster
    properties (mass, age, extinction) from a set of input photometry
    in various bands.

    Properties
       priors : array, shape (N) | callable | None
          prior probability on each data point; interpretation
          depends on the type passed; array, shape (N): values are
          interpreted as the prior probability of each data point;
          callable: the callable must take as an argument an array
          of shape (N, nphys), and return an array of shape (N)
          giving the prior probability at each data point; None:
          all data points have equal prior probability
       abstol : float
          absolute error tolerance for kernel density estimation
       reltol : float
          relative error tolerance for kernel density estimation
       thread_safe : bool
          if True, the computation routines will run in thread-safe
          mode, allowing use with multiprocessing; this incurs a small
          performance penalty

    Methods
       filters() : 
          returns list of filters available in the library
       filtersets() :
          return a list of the currently-loaded filter sets
       filter_units() :
          returns units for available filters
       add_filters() : 
          adds a set of filters for use in parameter estimation
       logL() : 
          compute log likelihood at a particular set of physical and
          photometric parameters
       mpdf() : 
          computer marginal posterior probability distribution for a
          set of photometric measurements
       mcmc() :
          due MCMC estimation of the posterior PDF on a set of
          photometric measurments
       bestmatch() : 
          find the simulations in the library that are the closest
          matches to the input photometry
       make_approx_phot() :
          given a set of physical properties, return a set of points
          that can be used for fast approximation of the corresponding
          photometric properties
       make_approx_phys() :
          given a set of photometric properties, return a set of points
          that can be used for fast approximation of the corresponding
          physical properties
       draw_sample():
          draw a random sample of cluster physical and photometric
          properties
       draw_phot():
          given a set of physical properties, return a
          randomly-selected set of photometric properties
    """

    ##################################################################
    # Initializer method
    ##################################################################
    def __init__(self, libname=None, filters=None, photsystem=None,
                 lib=None, bw_phys=0.1, bw_phot=None, ktype='gaussian', 
                 priors=None, sample_density=None, pobs=None, reltol=1.0e-2,
                 abstol=1.0e-8, leafsize=16, use_nebular=True,
                 use_extinction=True, thread_safe=True,
                 pruning=False, caching='none', vp_list=[]):
        """
        Initialize a cluster_slug object.

        Parameters
           libname : string
              name of the SLUG model to load; if left as None, the default
              is $SLUG_DIR/cluster_slug/modp020_chabrier_MW
           lib : object
              a library read by the read_cluster function; if specified
              this overrides the libname option; the library must
              contain both physical properties and photometry, and
              must include filter data; if this is not None, then the 
              photsystem keyword is ignored
           filters : iterable of stringlike
              list of filter names to be used for inferenence
           photsystem : None or string
              If photsystem is None, the library will be left in
              whatever photometric system was used to write
              it. Alternately, if it is a string, the data will be
              converted to the specified photometric system. Allowable
              values are 'L_nu', 'L_lambda', 'AB', 'STMAG', and
              'Vega', corresponding to the options defined in the SLUG
              code. Once this is set, any subsequent photometric data
              input are assumed to be in the same photometric system.
           bw_phys : 'auto' | float | array, shape (2) | array, shape (3)
              bandwidth for the physical quantities in the kernel
              density estimation; if set to 'auto', the bandwidth will
              be estimated automatically; if set to a scalar quantity,
              this will be used for all physical quantities; if set to
              an array, the array must have 2 elements if
              use_extinction is False, or 3 if it is True
           bw_phot : None | 'auto' | float | array
              bandwidth for the photometric quantities; if set to
              None, defaults to 0.25 mag / 0.1 dex; if set to 'auto',
              bandwidth is estimated automatically; if set to a float,
              this bandwidth is used for all photometric dimensions;
              if set to an array, the array must have the same number
              of dimensions as len(filters)
           ktype : string
              type of kernel to be used in densty estimation; allowed
              values are 'gaussian' (default), 'epanechnikov', and
              'tophat'; only Gaussian can be used with error bars
           priors : array, shape (N) | callable | None
              prior probability on each data point; interpretation
              depends on the type passed; array, shape (N): values are
              interpreted as the prior probability of each data point;
              callable: the callable must take as an argument an array
              of shape (N, nphys), and return an array of shape (N)
              giving the prior probability at each data point; None:
              all data points have equal prior probability
           pobs : array, shape (N) | callable | None
              probability of being observed on each data point; interpretation
              depends on the type passed; array, shape (N): values are
              interpreted as the prior probability of each data point;
              callable: the callable must take as an argument an array
              of shape (N, nphys), and return an array of shape (N)
              giving the prior probability at each data point; None:
              all data points have equal prior probability
           sample_density : array, shape (N) | callable | 'auto' | None
              the density of the data samples at each data point; this
              need not match the prior density; interpretation depends
              on the type passed; array, shape (N): values are
              interpreted as the density of data sampling at each
              sample point; callable: the callable must take as an
              argument an array of shape (N, nphys), and return an
              array of shape (N) giving the sampling density at each
              point; 'auto': the sample density will be computed
              directly from the data set; note that this can be quite
              slow for large data sets, so it is preferable to specify
              this analytically if it is known; None: data are assumed
              to be uniformly sampled, or to be sampled as the default
              library is if libname is also None
           reltol : float
              relative error tolerance; errors on all returned
              probabilities p will satisfy either
              abs(p_est - p_true) <= reltol * p_est   OR
              abs(p_est - p_true) <= abstol,
              where p_est is the returned estimate and p_true is the
              true value
           abstol : float
              absolute error tolerance; see above
           leafsize : int
              number of data points in each leaf of the KD tree
           use_nebular : bool
              if True, photometry including nebular emission will be
              used if available; if not, nebular emission will be
              omitted
           use_extinction : bool
              if True, photometry including extinction will be used;
              if not, it will be omitted, and in this case no results
              making use of the A_V dimension will be available
           thread_safe : bool
              if True, cluster_slug will make extra copies of internals
              as needed to ensure thread safety when the computation
              routines (logL, mpdf, mcmc, bestmatch, make_approx_phot,
              make_approx_phys, mpdf_approx) are used with
              multiprocessing; this incurs a minor performance
              penalty, and can be disabled by setting to False if the
              code will not be run with the multiprocessing module
           pruning : bool
              if True, the underlying bayesphot objects will be pruned
              of clusters for which pobs or priors are zero, speeding
              up evaluations; the current implementation is limited in
              that pruning for priors is only done when the
              cluster_slug object is instantiated, and pruning for
              pobs is only done when each filter set is added, so the
              list of pruned clusters cannot be modified later
           caching : 'aggressive' | 'lazy' | 'none'
              strategy for caching subsets of the data with some
              dimensions marginalised out; behavior is as follows:

                 'agressive'
                    on construction, store sorted data for fast
                    calculation of 1D PDFs of variables by themselves,
                    and 1D PDFs of all physical variables marginalised
                    over all other physical variables; this
                    significantly increases the memory footprint and
                    construction time, but greatly speeds up
                    subsequent evaluation of any of these quantities,
                    and is generally the best choice for prudction
                    work in parallel
                 'lazy'
                    sorted data sets for fast computation are created
                    and cached as needed; this will make the first
                    computation of any marginal PDF slower, but speed
                    up all subsequent ones, without imposing any
                    extra time at initial construction; this mode is
                    generally best for interactive work, but is not
                    thread_safe = True; memory cost depends on how
                    many different marginal PDF combinations are
                    calculated, but is always less than aggressive
                 'none'
                    no caching is performed automatically; the user
                    may still manually cache data by calling
                    the make_cache method
           vp_list : list
              A list with an element for each of the variable parameters
              in the data. An element is set to True if we wish to use
              that parameter here, or False if we do not.

        Returns
           Nothing

        Raises
           IOError, if the library cannot be found
        """

        # If using the default library, assign the library name
        if libname is None and lib is None:
            self.__libname = osp.join('cluster_slug', 'modp020_chabrier_MW')
            if 'SLUG_DIR' in os.environ:
                self.__libname = osp.join(os.environ['SLUG_DIR'], 
                                          self.__libname)
        else:
            self.__libname = libname

        # If using a library passed in, store a pointer to is
        if lib is not None:
            self.__lib = lib
            self.__libname = None
        else:
            self.__lib = None

        # Load the cluster physical properties
        try:
            if self.__lib is None:
                prop = read_cluster_prop(self.__libname)
            else:
                prop = self.__lib
        except IOError:

            # If we're here, we failed to load the library. If we were
            # given a library file name explicitly, just raise an
            # error.
            if libname is not None:
                raise IOError(errno.ENOENT,
                              "unable to open library {}".
                              format(self.__libname))

            # If we've made it to here, we were asked to open the
            # default library but failed to do so. Check if the
            # failure could be because we don't have astropy and thus
            # can't open fits files. If that's the cause, print out a
            # helpful error message.
            try:
                import astropy.io.fits as fits
            except ImportError:
                raise IOError(errno.EIO,
                              "failed to read default cluster_slug " +
                              "library cluster_slug/clusterslug_mw " +
                              "due to missing " +
                              "astropy.io.fits; install astropy " +
                              "or specify a library in a non-FITS " +
                              "format")

            # If we're here, we couldn't open the default library, and
            # it's not because we don't have FITS capability. The file
            # must not exist, or must be damaged. Check if we're in
            # interactive mode. If not, just raise an error and
            # suggest the user to go get the library file.
            errstr = "Unable to open default cluster_slug " + \
                     "library file cluster_slug/clusterslug_mw."
            import __main__ as main
            if hasattr(main, '__file__'):
                # We're not interactive; just raise an error
                raise IOError(errno.EIO, errstr + " " +
                              "Try downloading it from " +
                              "https://sites.google.com/site/runslug/data")

            # If we're here, we don't have the library file, but we
            # are in interactive mode. Thus offer the user an option
            # to go download the file now.
            usr_response \
                = raw_input(errstr + " Would you like to download it "
                            "now (warning: 20 GB)? [y/n] ").\
                lower().strip()
            if not usr_response in ['yes', 'y', 'ye']:
                # User didn't say yes, so raise error
                raise IOError(errno.EIO, "Unable to proceeed")

            # If we're here, download the files
            print("Fetching modp020_chabrier_MW_cluster_prop " +
                  "(this may take a while)...")
            url = urllib2.urlopen(
                'https://www.dropbox.com/s/tg2ogad713nfux0/modp020_chabrier_MW_cluster_prop.fits?dl=0')
            rawdata = url.read()
            url.close()
            fp = open(osp.join(osp.dirname(self.__libname),
                               'modp020_chabrier_MW_cluster_phot.fits'), 'wb')
            fp.write(rawdata)
            fp.close()
            print("Fetching modp020_chabrier_MW_cluster_phot.fits " +
                  "(this make take a while)...")
            url = urllib2.urlopen(
                'https://www.dropbox.com/s/inmeeuldkxsazzs/'
                'modp020_chabrier_MW_cluster_phot.fits?dl=0')
            rawdata = url.read()
            url.close()
            fp = open(osp.join(osp.dirname(self.__libname),
                               'modp020_chabrier_MW_cluster_phot.fits'), 'wb')
            fp.write(rawdata)
            fp.close()

            # Now try reading the data
            try:
                prop = read_integrated_prop(self.__libname)
            except IOError:
                raise IOError(errno.EIO,
                              "still unable to open default library")
             
        # Find out how many variable parameters we have and how many
        # we want to use
        nvp_total = np.size(vp_list)
        nvp = np.size(np.where(vp_list))
        if nvp_total > 0:
            print("Total Number of Variable Parameters in data: "
                  +str(nvp_total))
            print("Number of Variable Parameters to use: "+str(nvp))

        # Set number of physical properties
        if use_extinction:
            self.__nphys = 3 + nvp
        else:
            self.__nphys = 2 + nvp
            
        # Store the physical properties
        self.__ds_phys = np.zeros((len(prop.id), self.__nphys))
        self.__ds_phys[:,0] = np.log10(prop.actual_mass)
        self.__ds_phys[:,1] = np.log10(prop.time - prop.form_time)
        pnum = 1
        if use_extinction:
            self.__ds_phys[:,2] = prop.A_V
            pnum = 2
            
        # Grab the variable parameters we want to use
        for vpi in range (0,nvp_total,1):
            if vp_list[vpi] is True:
                self.__ds_phys[:,pnum+(vpi+1)] = getattr(prop, "VP"+repr(vpi))

        # Record available filters
        if self.__lib is None:
            filter_info = read_cluster_phot(self.__libname,
                                            filters_only=True, 
                                            nofilterdata=True)
        else:
            filter_info = self.__lib
        self.__allfilters = filter_info.filter_names
        self.__allunits = filter_info.filter_units

        # Record other stuff that we'll use later
        if self.__lib is None:
            self.__photsystem = photsystem
        else:
            self.__photsystem = None
        self.__use_nebular = use_nebular
        self.__use_extinction = use_extinction
        self.__ktype = ktype
        self.__priors = priors
        if (sample_density is not None) or \
           (libname is not None):
            self.__sample_density = sample_density
        else:
            self.__sample_density = _default_sample_density
        self.__reltol = reltol
        self.__abstol = abstol
        self.__bw_phot_default = bw_phot
        self.__thread_safe = thread_safe
        self.__pruning = pruning

        # If we are pruning, use the priors object to figure out which
        # clusters should be pruned
        if self.__pruning:
            if self.__priors is not None:
                if hasattr(self.__priors, '__call__'):
                    pr = self.__priors(self.__ds_phys)
                    self.__prior_nonzero = pr > 0.0
                else:
                    self.__prior_nonzero = priors > 0.0
            else:
                self.__prior_nonzero = np.ones(len(self.__ds_phys),
                                               dtype=np.bool)

        # Set the physical bandwidth
        self.__bw_phys = copy.deepcopy(bw_phys)

        # Initialize list of photometric data we've read to empty dict
        self.__photdata = {}
        self.__photbw = {}

        # Initialize an empty list of filter sets
        self.__filtersets = []

        # Save caching mode
        self.__caching = caching

        # If we have been given a filter list, create the data set to
        # go with it
        if filters is not None:
            self.add_filters(filters, bandwidth=bw_phot, pobs=pobs)


    ##################################################################
    # Method to load the data off disk for a particular filter, or to
    # set up arrays pointing to data if using a library stored in
    # memory
    ##################################################################
    def load_data(self, filter_name, bandwidth=None, 
                  force_reload=False):
        """
        Loads photometric data for the specified filter into memory

        Parameters:
           filter_name : string
              name of filter to load
           bandwidth : float
              default bandwidth for this filter
           force_reload : bool
              if True, reinitialize the data even if has already been
              stored

        Returns:
           None

        Raises:
           ValueError, if filter_name is not one of the available
           filters
        """

        # Do nothing if data has already been read, unless we've been
        # told to force a re-read
        if filter_name in self.__photdata.keys() and not force_reload:
            return

        # Make sure we have this filter; if not, raise error
        if filter_name not in self.__allfilters:
            raise ValueError("no data available for filter {}".
                             format(filter_name))

        # Suppress obnoxious numpy warning messages here
        errstate = np.geterr()
        np.seterr(divide='ignore', invalid='ignore', over='ignore',
                  under='ignore')

        # Special case: for ionizing fluxes, always load the
        # non-nebular, non-extincted value, and don't do any
        # photometric system conversion
        if filter_name == 'QH0' or filter_name == 'QHe0' or \
           filter_name == 'QHe1':
            if self.__lib is None:
                phot = read_cluster_phot(self.__libname, 
                                           read_filters=filter_name,
                                           read_nebular=False,
                                           read_extinct=False,
                                           phot_only=True)
                phdata = np.squeeze(phot.phot)
            else:
                phot = self.__lib
                phdata = np.squeeze(
                    np.copy(self.__lib.phot
                            [:,self.__allfilters.index(filter_name)]))
        else:

            # Load data; first try reading the requested combination
            # of nebular and extincted values
            if self.__lib is None:
                phot = read_cluster_phot(self.__libname, 
                                         read_filters=filter_name,
                                         read_nebular=self.__use_nebular,
                                         read_extinct=self.__use_extinction,
                                         phot_only=True,
                                         photsystem=self.__photsystem)
            else:
                phot = self.__lib

            # Make sure we have the data we want; if not, lots of ugly
            # special cases as fallbacks
            if (self.__use_nebular and self.__use_extinction):
                # Wanted nebular + extinction
                if 'phot_neb_ex' in phot._fields:
                    # Got it
                    if self.__lib is None:
                        phdata = np.squeeze(phot.phot_neb_ex)
                    else:
                        phdata = np.copy(phot.phot_neb_ex
                                         [:,self.__allfilters.index(filter_name)])
                    warn_nebular = False
                    warn_extinct = False
                    nebular = True
                    extinct = True
                else:
                    # Didn't get it; try without nebular
                    if self.__lib is None:
                        phot = read_cluster_phot(
                            self.__libname, 
                            read_filters=filter_name,
                            read_nebular=False,
                            read_extinct=True,
                            phot_only=True,
                            photsystem=self.__photsystem)
                    if 'phot_ex' in phot._fields:
                        if self.__lib is None:
                            phdata = np.squeeze(phot.phot_ex)
                        else:
                            phdata = np.copy(phot.phot_ex
                                             [:,self.__allfilters.index(filter_name)])
                        warn_nebular = True
                        warn_extinct = False
                        nebular = False
                        extinct = True
                    else:
                        # Couldn't get extinction only; try nebular only
                        if self.__lib is None:
                            phot = read_cluster_phot(
                                self.__libname, 
                                read_filters=filter_name,
                                read_nebular=True,
                                read_extinct=False,
                                phot_only=True,
                                photsystem=self.__photsystem)
                        if 'phot_neb' in phot._fields:
                            if self.__lib is None:
                                phdata = np.squeeze(phot.phot_neb)
                            else:
                                phdata = np.copy(phot.phot_neb
                                                 [:,self.__allfilters.index(filter_name)])
                            warn_nebular = False
                            warn_extinct = True
                            nebular = True
                            extinct = False
                        else:
                            # Ultimate fallback: no nebular or extinction
                            if self.__lib is None:
                                phot = read_cluster_phot(
                                    self.__libname, 
                                    read_filters=filter_name,
                                    read_nebular=False,
                                    read_extinct=False,
                                    phot_only=True,
                                    photsystem=self.__photsystem)
                                phdata = np.squeeze(phot.phot)
                            else:
                                phdata = np.copy(phot.phot
                                                 [:,self.__allfilters.index(filter_name)])
                            warn_nebular = True
                            warn_extinct = True
                            nebular = False
                            extinct = False

            elif self.__use_nebular:
                # Want nebular, no extinction
                if 'phot_neb' in phot._fields:
                    # Got it
                    if self.__lib is None:
                        phdata = np.squeeze(phot.phot_neb)
                    else:
                        phdata = np.copy(phot.phot_neb
                                         [:,self.__allfilters.index(filter_name)])
                    warn_nebular = False
                    warn_extinct = False
                    nebular = True
                    extinct = False
                else:
                    # Didn't get it; go to no nebular or extinction
                    if self.__lib is None:
                        phot = read_cluster_phot(
                            self.__libname, 
                            read_filters=filter_name,
                            read_nebular=False,
                            read_extinct=False,
                            phot_only=True,
                            photsystem=self.__photsystem)
                        phdata = np.squeeze(phot.phot)
                    else:
                        phdata = np.copy(phot.phot
                                         [:,self.__allfilters.index(filter_name)])
                    warn_nebular = True
                    warn_extinct = False
                    nebular = False
                    extinct = False

            elif self.__use_extinction:
                # Wanted extinction, no nebular
                if 'phot_ex' in phot._fields:
                    # Got it
                    if self.__lib is None:
                        phdata = np.squeeze(phot.phot_ex)
                    else:
                        phdata = np.copy(phot.phot_ex
                                         [:,self.__allfilters.index(filter_name)])
                    warn_nebular = False
                    warn_extinct = False
                    nebular = False
                    extinct = True
                else:
                    # Didn't get it; go to no nebular or extinction
                    if self.__lib is None:
                        phot = read_cluster_phot(
                            self.__libname, 
                            read_filters=filter_name,
                            read_nebular=False,
                            read_extinct=False,
                            phot_only=True,
                            photsystem=self.__photsystem)
                        phdata = np.squeeze(phot.phot)
                    else:
                        phdata = np.copy(phot.phot
                                         [:,self.__allfilters.index(filter_name)])
                    warn_nebular = False
                    warn_extinct = True
                    nebular = False
                    extinct = False
            elif self.__use_extinction is False and self.__use_nebular is False:
                # Didn't get it; go to no nebular or extinction
                if self.__lib is None:
                    phot = read_cluster_phot(
                        self.__libname, 
                        read_filters=filter_name,
                        read_nebular=False,
                        read_extinct=False,
                        phot_only=True,
                        photsystem=self.__photsystem)
                    phdata = np.squeeze(phot.phot)
                else:
                    phdata = np.copy(phot.phot
                                     [:,self.__allfilters.index(filter_name)])
                warn_nebular = False
                warn_extinct = False
                nebular = False
                extinct = False

            # Make sure data aren't NaN's. If they are, that means we
            # don't have extinction for this filter, so we need to use
            # non-extincted data.
            if np.isnan(phdata[0]):
                warn_extinct = True
                if nebular:
                    if self.__lib is None:
                        phot = read_cluster_phot(
                            self.__libname, 
                            read_filters=filter_name,
                            read_nebular=True,
                            read_extinct=False,
                            phot_only=True,
                            photsystem=self.__photsystem)
                        phdata = np.squeeze(phot.phot_neb)
                    else:
                        phdata = np.copy(phot.phot_neb
                                         [:,self.__allfilters.index(filter_name)])
                else:
                    if self.__lib is None:
                        phot = read_cluster_phot(
                            self.__libname, 
                            read_filters=filter_name,
                            read_nebular=False,
                            read_extinct=False,
                            phot_only=True,
                            photsystem=self.__photsystem)
                        phdata = np.squeeze(phot.phot)
                    else:
                        phdata = np.copy(phot.phot
                                         [:,self.__allfilters.index(filter_name)])
            else:
                warn_extinct = False

            # Issue warnings if necessary
            if warn_nebular:
                warnstr = ("cluster_slug: nebular data requested for "+ \
                           "filter {}, but is not available; using non-"+ \
                           "nebular data instead").format(filter_name)
                warnings.warn(warnstr)
            if warn_extinct:
                warnstr = ("cluster_slug: extincted data requested for "+ \
                           "filter {}, but is not available; using non-"+ \
                           "extinced data instead").format(filter_name)
                warnings.warn(warnstr)

        # Take the log of the photometric values if they're recorded
        # in a linear system; also set the bandwidth based on whether
        # we're in a magnitude system or not; we can override this
        # later if we want
        if 'mag' not in phot.filter_units[0]:
            phdata[phdata <= 0] = 1.0e-99
            phdata = np.log10(phdata)
            if bandwidth is None:
                if self.__bw_phot_default is None:
                    self.__photbw[filter_name] = 0.1
                else:
                    self.__photbw[filter_name] = self.__bw_phot_default
            else:
                self.__photbw[filter_name] = bandwidth
        else:
            if bandwidth is None:
                if self.__bw_phot_default is None:
                    self.__photbw[filter_name] = 0.25
                else:
                    self.__photbw[filter_name] = self.__bw_phot_default
            else:
                self.__photbw[filter_name] = bandwidth

        # Fix any inf's that were generated by photometric system
        # conversions or taking logs
        phdata[np.isinf(phdata)] = 99.0

        # Store the data
        self.__photdata[filter_name] = phdata

        # Restore the numpy error state
        np.seterr(divide=errstate['divide'], over=errstate['over'], 
                  under=errstate['under'], invalid=errstate['invalid'])


    ##################################################################
    # Method to return information on available filters
    ##################################################################
    def filters(self):
        """
        Returns list of all available filters

        Parameters:
           None

        Returns:
           filters : list of strings
              list of available filter names
        """
        return copy.deepcopy(self.__allfilters)

    def filter_units(self):
        """
        Returns list of all available filter units

        Parameters:
           None

        Returns:
           units : list of strings
              list of available filter units
        """
        return copy.deepcopy(self.__allunits)

    def filtersets(self):
        """
        Returns list of all currently-loaded filter sets

        Parameters:
           None

        Returns:
           filtersets : list of list of strings
              list of currently-loaded filter sets
        """
        return copy.deepcopy([f['filters'] for 
                              f in self.__filtersets])


    ##################################################################
    # Method to prepare to analyze a particular set of filters
    ##################################################################
    def add_filters(self, filters, bandwidth=None, pobs=None):
        """
        Add a set of filters to use for cluster property estimation

        Parameters
           filters : iterable of stringlike
              list of filter names to be used for inferenence
           bandwidth : None | 'auto' | float | array
              bandwidth for the photometric quantities; if set to
              None, the bandwidth is unchanged for an existing filter
              set, and for a newly-created one the default physical
              and photometric bandwidths are used; if set to 'auto',
              bandwidth is estimated automatically; if set to a float,
              this bandwidth is used for all physical and photometric
              dimensions; if set to an array, the array must have the
              same number of entries as nphys+len(filters)
           pobs : array, shape (N) | callable | 'equal' | None
              the probability that a particular object would be observed,
              which is used, like prior, to weight the library;
              interpretation depends on type. 'equal' means all objects are
              equally likely to be observed, array is an array giving the
              observation probability of each object in the library, and
              callable means must be a function that takes an array
              containing the photometry, of shape (N, nhpot), as an
              argument, and returns an array of shape (N) giving the
              probability of observation for that object. Finally,
              None leaves the observational probability unchanged

        Returns
           nothing
        """

        # If we already have this filter set in our dict, just set the
        # bandwidth and pobs and return
        for i, f in enumerate(self.__filtersets):
            if filters == f['filters']:
                if bandwidth is not None:
                    self.__filtersets[i]['bp'].bandwidth = bandwidth
                if pobs is not None:
                    if pobs is 'equal':
                        self.__filtersets[i]['bp'].pobs = None
                    elif not self.__pruning or hasattr(pobs, '__call__'):
                        self.__filtersets[i]['bp'].pobs = pobs
                    else:
                        pobs_tmp = pobs[f['keep']]
                        self.__filtersets[i]['bp'].pobs = pobs_tmp
                return

        # We're adding a new filter set, so save its name
        newfilter = { 'filters' : copy.deepcopy(filters) }

        # Construct data set to use with this filter combination, and
        # add the physical property data
        dataset = np.zeros((self.__ds_phys.shape[0], 
                                         self.__nphys+len(filters)))
        dataset[:,:self.__nphys] = self.__ds_phys
        
        # Loop over filters
        for i, f in enumerate(filters):

            # Do we have this filter loaded already? If not, read it.
            if f not in self.__photdata.keys():
                if bandwidth is None or type(bandwidth) is str:
                    self.load_data(f)
                else:
                    if hasattr(bandwidth, '__iter__'):
                        self.load_data(f, bandwidth=bandwidth[i])
                    else:
                        self.load_data(f, bandwidth=bandwidth)

            # Add data for this filter
            dataset[:,self.__nphys+i] = self.__photdata[f]

        # Make copies of the priors and sample_density for use by the
        # derived bayesphot class; we make copies because bp needs to
        # be able to sort these
        priors = deepcopy(self.__priors)
        sample_density = deepcopy(self.__sample_density)
        
        # Are we pruning?
        if self.__pruning:

            # Get pruning based on p_obs
            if hasattr(pobs, '__call__'):
                po = pobs(dataset[:,self.__nphys:])
                pobs_nonzero = po > 0.0
            elif pobs is not None:
                pobs_nonzero = pobs > 0.0
            else:
                pobs_nonzero = np.ones(len(dataset),
                                       dtype=np.bool)

            # Get list of clusters with p_obs > 0 and prior > 0; store
            # this because we will need it to change p_obs or prior
            # later
            newfilter['keep'] = np.logical_and(self.__prior_nonzero, 
                                               pobs_nonzero)

            # Prune the data set
            dataset = dataset[newfilter['keep']]

            # Make pruned versions of sample_density, priors, and pobs
            # if they are not callables; if they are callables, they
            # will be computed on the fly by the bp object, and we
            # don't have to worry about it
            if not hasattr(sample_density, '__call__') and \
               sample_density is not None:
                sample_density = sample_density[newfilter['keep']]
            if not hasattr(priors, '__call__') and priors is not None:
                priors = priors[newfilter['keep']]
            if not hasattr(pobs, '__call__') and pobs is not None:
                pobs_tmp = pobs[newfilter['keep']]
            else:
                pobs_tmp = pobs

        else:
            
            # Just set pointer to pobs if we're not doing any pruning
            pobs_tmp = pobs

        # Store the new data set
        newfilter['dataset'] = dataset
        
        # Set bandwidth
        if self.__bw_phys == 'auto' or bandwidth == 'auto':
            bw = 'auto'
        else:
            bw = np.zeros(self.__nphys+len(filters))
            bw[:self.__nphys] = self.__bw_phys
            if bandwidth is not None:
                bw[self.__nphys:] = bandwidth
            else:
                for i in range(len(filters)):
                    bw[self.__nphys+i] = self.__photbw[f]

        # Build the bp object
        newfilter['bp'] \
            = bp(newfilter['dataset'],
                 self.__nphys,
                 filters = filters,
                 bandwidth = bw,
                 ktype = self.__ktype,
                 priors = priors,
                 pobs = pobs_tmp,
                 sample_density = sample_density,
                 reltol = self.__reltol,
                 abstol = self.__abstol,
                 thread_safe = self.__thread_safe,
                 caching = self.__caching)

        # Save to the master filter list
        self.__filtersets.append(newfilter)


    ##################################################################
    # Method to delete a set of filters
    ##################################################################
    def del_filters(self, filters):
        """
        Remove a set of filters, freeing the memory associated with
        them. Note that this does not delete the underlying library
        data, just the data for the KD tree used internally.

        Parameters
           filters : iterable of stringlike
              list of filter names

        Returns
           Nothing

        Raises
           KeyError if the input set of filters is not loaded
        """

        # Figure out which filter set to delete
        idel = -1
        for i, f in enumerate(self.__filtersets):
            if filters == f['filters']:
                idel = i
                break
        if idel == -1:
            raise KeyError("filter set not found")

        # Delete filter set; destroy the bayesphot object explicitly
        # to ensure that it is remove from scope immediately, since it
        # can be rather large
        del self.__filtersets[idel]['bp']
        self.__filtersets.pop(idel)


    ##################################################################
    # Define the priors property. This just wraps around the
    # corresponding property defined for bp objects, but we do have to
    # be careful about pruning if it is turned on
    ##################################################################
    @property
    def priors(self):
        return self.__priors

    @priors.setter
    def priors(self, pr):
        self.__priors = pr
        if not self.__pruning or hasattr(pr, '__call__') or \
           pr is None:
            for f in self.__filtersets:
                f['bp'].priors = deepcopy(self.__priors)
        else:
            # If the data have been pruned before passing to the bp
            # object, we need to apply the same pruning to the new
            # priors
            for f in self.__filtersets:
                pr_tmp = pr[f['keep']]
                f['bp'].priors = pr_tmp

    ##################################################################
    # Define the pobs property. This wraps around the corresponding
    # property for bp objects, with the complication that if there's
    # more than one filter set is has to accept or return a list, and
    # that we have to handle pruning
    ##################################################################
    @property
    def pobs(self):
        po = []
        for f in self.__filtersets:
            po.append(f['bp'].pobs)
        return po

    @pobs.setter
    def pobs(self, po):
        # If po is an iterable and matches the number of filter sets,
        # we assume that this is an observation probability to be
        # assigned to each filter set. If not, we try assigning it to
        # each filter set
        if hasattr(po, '__iter__'):
            if len(po) == len(self.__filtersets):
                for f, p in zip(self.__filtersets, po):
                    if not self.__pruning or hasattr(p, '__call__'):
                        f['bp'].pobs = p
                    else:
                        p_tmp = p[f['keep']]
                        f['bp'].pobs = p_tmp
                return
        else:
            if not self.__pruning or hasattr(po, '__call__') \
               or po is None:
                for f in self.__filtersets:
                    f['bp'].pobs = deepcopy(po)
            else:
                for f in self.__filtersets:
                    po_tmp = deepcopy(po)[f['keep']]
                    f['bp'].pobs = po_tmp
                        
    ##################################################################
    # Define properties that update the current values for the
    # cluster_slug object, and also update all the child bp objects.
    ##################################################################
    @property
    def abstol(self):
        return self.__abstol

    @abstol.setter
    def abstol(self, newtol):
        self.__abstol = newtol
        for f in self.__filtersets:
            f['bp'].abstol = self.__abstol

    @property
    def reltol(self):
        return self.__reltol

    @reltol.setter
    def reltol(self, newtol):
        self.__reltol = newtol
        for f in self.__filtersets:
            f['bp'].reltol = self.__reltol

    @property
    def thread_safe(self):
        return self.__thread_safe

    @thread_safe.setter
    def thread_safe(self, new_thread_safe):
        self.__thread_safe = new_thread_safe
        for f in self.__filtersets:
            f['bp'].thread_safe = self.__thread_safe

    ##################################################################
    # Methods to make and destroy caches
    ##################################################################
    def make_cache(self, margindims, filters=None):
        """
        This method builds a cache to do faster calculation of PDFs
        where certain dimensions are marginalised out. If such caches
        exist, they are used automatically by all the computation
        methods.

        Parameters
           margindims : listlike of integers
              list of dimensions to be marginalised out; physical
              dimensions go from 0 - nphys-1, photometric dimensions
              from nphys to nphys + nphot - 1; note that the indexing
              depends on the filter set specified by filters
           filters : listlike of strings
              list of photometric filters to use; if left as None, and
              only 1 set of photometric filters has been defined for
              the cluster_slug object, that set will be used by
              default

        Returns
           Nothing
        """

        # Were we given a set of filters?
        if filters is None:

            # No filters given; if we have only a single filter set
            # stored, just use it
            if len(self.__filtersets) == 1:
                return self.__filtersets[0]['bp']. \
                    make_cache(margindims)
            else:
                raise ValueError("must specify a filter set")

        else:

            # We were given a filter set; add it if it doesn't exist
            self.add_filters(filters)

            # Find the bp object we should use
            for f in self.__filtersets:
                if f['filters'] == filters:
                    bp = f['bp']
                    break

            # Make cache
            bp.make_cache(margindims)

    def clear_cache(self, margindims=None):
        """
        This method deletes from the cache

        Parameters
           margindims : listlike of integers
              list of marginalised dimensions that should be removed
              from the cache, in the same format as make_cache; if
              left as None, the cache is completely emptied
           filters : listlike of strings
              list of photometric filters to use; if left as None, and
              only 1 set of photometric filters has been defined for
              the cluster_slug object, that set will be used by
              default

        Returns
           Nothing
        """
        # Were we given a set of filters?
        if filters is None:

            # No filters given; if we have only a single filter set
            # stored, just use it
            if len(self.__filtersets) == 1:
                return self.__filtersets[0]['bp']. \
                    clear_cache(margindims)
            else:
                raise ValueError("must specify a filter set")

        else:

            # We were given a filter set; add it if it doesn't exist
            self.add_filters(filters)

            # Find the bp object we should use
            for f in self.__filtersets:
                if f['filters'] == filters:
                    bp = f['bp']
                    break

            # Make cache
            bp.clear_cache(margindims)
            

    ##################################################################
    # The functions below just wrap around the bayesphot functions of
    # the same name. They all have in common that they accept an
    # additional keyword argument, filters, which specifies the filter
    # set to use.
    ##################################################################    
    def logL(self, physprop, photprop, photerr=None, filters=None,
             margindim=None):
        """
        This function returns the natural log of the likelihood
        function evaluated at a particular log mass, log age,
        extinction, and set of log luminosities

        Parameters:
           physprop : arraylike, shape (nhpys) or (..., nphys)
              array giving values of the log M, log T, and A_V; for a
              multidimensional array, the operation is vectorized over
              the leading dimensions; if created with use_extinct =
              False, the A_V dimension should be omitted.
              Will also include variable any variable parameters VPx if
              they are requested.
           photprop : arraylike, shape (nfilter) or (..., nfilter)
              array giving the photometric values; for a
              multidimensional array, the operation is vectorized over
              the leading dimensions
           photerr : arraylike, shape (nfilter) or (..., nfilter)
              array giving photometric errors; for a multidimensional
              array, the operation is vectorized over the leading
              dimensions
           filters : listlike of strings
              list of photometric filters to use; if left as None, and
              only 1 set of photometric filters has been defined for
              the cluster_slug object, that set will be used by
              default
           margindim : int | arraylike of ints | None
              The index or indices of the physical or photometric
              properties to be maginalized over, numbered from 0 -
              nphys-1 for physical properties and from nphys - nfilter +
              nfilter - 1 for photometric properties. If this keyword is
              set, then physprop and/or photprop should have fewer
              than nphys or nphot elements due to the omission of
              marginalised dimensions. If all physical or photometric
              dimensions are marginalised out, that corresponding
              argument for physprop or photprop should be set to None

        Returns:
           logL : float or arraylike
              natural log of the likelihood function
        """

        # Were we given a set of filters?
        if filters is None:

            # No filters given; if we have only a single filter set
            # stored, just use it
            if len(self.__filtersets) == 1:
                return self.__filtersets[0]['bp']. \
                    logL(physprop, photprop, photerr=photerr,
                         margindim=margindim)
            else:
                raise ValueError("must specify a filter set")

        else:

            # We were given a filter set; add it if it doesn't exist
            self.add_filters(filters)

            # Find the bp object we should use
            for f in self.__filtersets:
                if f['filters'] == filters:
                    bp = f['bp']
                    break

            # Call the logL method
            return bp.logL(physprop, photprop, photerr=photerr,
                           margindim=margindim)


    def mpdf(self, idx, photprop, photerr=None, ngrid=128,
             qmin=None, qmax=None, grid=None, norm=True,
             filters=None):
        """
        Returns the marginal probability for one or mode physical
        quantities for one or more input sets of photometric
        properties. Output quantities are computed on a grid of
        values, in the same style as meshgrid

        Parameters:
           idx : int or listlike containing ints
              index of the physical quantity whose PDF is to be
              computed; 0 = log M, 1 = log T, 2 = A_V, (2 or 3)+x = VPx; 
              if this is an iterable, the joint distribution of the indicated
              quantities is returned
           photprop : arraylike, shape (nfilter) or (..., nfilter)
              array giving the photometric values; for a
              multidimensional array, the operation is vectorized over
              the leading dimensions
           photerr : arraylike, shape (nfilter) or (..., nfilter)
              array giving photometric errors; for a multidimensional
              array, the operation is vectorized over the leading
              dimensions
           ngrid : int or listlike containing ints
              number of points in each dimension of the output grid;
              if this is an iterable, it must have the same number of
              elements as idx
           qmin : float or listlike
              minimum value in the output grid in each quantity; if
              left as None, defaults to the minimum value in the
              library; if this is an iterable, it must contain the
              same number of elements as idx
           qmax : float or listlike
              maximum value in the output grid in each quantity; if
              left as None, defaults to the maximum value in the
              library; if this is an iterable, it must contain the
              same number of elements as idx
           grid : listlike of arrays
              set of values defining the grid on which the PDF is to
              be evaluated, in the same format used by meshgrid
           norm : bool
              if True, returned pdf's will be normalized to integrate
              to 1
           filters : listlike of strings
              list of photometric filters to use; if left as None, and
              only 1 set of photometric filters has been defined for
              the cluster_slug object, that set will be used by
              default

        Returns:
           grid_out : array
              array of values at which the PDF is evaluated; contents
              are the same as returned by meshgrid
           pdf : array
              array of marginal posterior probabilities at each point
              of the output grid, for each input cluster; the leading
              dimensions match the leading dimensions produced by
              broadcasting the leading dimensions of photprop and
              photerr together, while the trailing dimensions match
              the dimensions of the output grid
        """

        # Were we given a set of filters?
        if filters is None:

            # No filters given; if we have only a single filter set
            # stored, just use it
            if len(self.__filtersets) == 1:
                return self.__filtersets[0]['bp']. \
                    mpdf(idx, photprop, photerr=photerr, 
                         ngrid=ngrid, qmin=qmin, qmax=qmax, 
                         grid=grid, norm=norm)
            else:
                raise ValueError("must specify a filter set")

        else:

            # We were given a filter set; add it if it doesn't exist
            self.add_filters(filters)

            # Find the bp object we should use
            for f in self.__filtersets:
                if f['filters'] == filters:
                    bp = f['bp']
                    break

            # Call the logL method
            return bp.mpdf(idx, photprop, photerr=photerr,
                           ngrid=ngrid, qmin=qmin, qmax=qmax,
                           grid=grid, norm=norm)


    def mpdf_phot(self, idx, physprop, ngrid=128,
                  qmin=None, qmax=None, grid=None, norm=True,
                  filters=None):
        """
        Returns the marginal probability for one or more photometric
        quantities corresponding to an input set or distribution of
        physical properties. Output quantities are computed on a grid of
        values, in the same style as meshgrid.

        Parameters:
           idx : int or listlike containing ints
              index of the photometric quantity whose PDF is to be
              computed, starting at 0; indices correspond to the order
              of elements in the filters argument; if this is an
              iterable, the joint distribution of the indicated
              quantities is returned
           physprop : arraylike, shape (nphys) or (..., nphys)
              physical properties to be used; if this is an array of
              nphys elements, these give the physical properties; if
              it is a multidimensional array, the operation is
              vectoried over the leading dimensions
              physical properties -- the function must take an array
              of (nphys) elements as an input, and return a floating
              point value representing the PDF evaluated at that set
              of physical properties as an output
           ngrid : int or listlike containing ints
              number of points in each dimension of the output grid;
              if this is an iterable, it must have the same number of
              elements as idx
           qmin : float or listlike
              minimum value in the output grid in each quantity; if
              left as None, defaults to the minimum value in the
              library; if this is an iterable, it must contain the
              same number of elements as idx
           qmax : float or listlike
              maximum value in the output grid in each quantity; if
              left as None, defaults to the maximum value in the
              library; if this is an iterable, it must contain the
              same number of elements as idx
           grid : listlike of arrays
              set of values defining the grid on which the PDF is to
              be evaluated, in the same format used by meshgrid
           norm : bool
              if True, returned pdf's will be normalized to integrate
              to 1
           filters : listlike of strings
              list of photometric filters to use; if left as None, and
              only 1 set of photometric filters has been defined for
              the cluster_slug object, that set will be used by
              default

        Returns:
           grid_out : array
              array of values at which the PDF is evaluated; contents
              are the same as returned by meshgrid
           pdf : array
              array of marginal posterior probabilities at each point
              of the output grid, for each input set of properties; the leading
              dimensions match the leading dimensions produced by
              broadcasting the leading dimensions of photprop and
              photerr together, while the trailing dimensions match
              the dimensions of the output grid
        """
        # Were we given a set of filters?
        if filters is None:

            # No filters given; if we have only a single filter set
            # stored, just use it
            if len(self.__filtersets) == 1:
                return self.__filtersets[0]['bp']. \
                    mpdf_phot(idx, physprop,
                              ngrid=ngrid, qmin=qmin, qmax=qmax, 
                              grid=grid, norm=norm)
            else:
                raise ValueError("must specify a filter set")

        else:

            # We were given a filter set; add it if it doesn't exist
            self.add_filters(filters)

            # Find the bp object we should use
            for f in self.__filtersets:
                if f['filters'] == filters:
                    bp = f['bp']
                    break

            # Call the logL method
            return bp.mpdf_phot(idx, physprop,
                                ngrid=ngrid, qmin=qmin, qmax=qmax,
                                grid=grid, norm=norm)

    def mpdf_gen(self, fixeddim, fixedprop, margindim, ngrid=128,
                 qmin=None, qmax=None, grid=None, norm=True,
                 filters=None):
        """
        Returns the marginal probability for one or more physical or
        photometric properties, keeping other properties fixed and
        marginalizing over other quantities. This is the most general
        marginal PDF routine provided.

        Parameters:
           fixeddim : int | arraylike of ints | None
              The index or indices of the physical or photometric
              properties to be held fixed; physical properties are
              numbered 0 ... nphys-1, and photometric ones are numbered
              nphys ... nphys + nphot - 1. This can also be set to
              None, in which case no properties are held fixed.
           fixedprop : array | None
              The values of the properties being held fixed; the size
              of the final dimension must be equal to the number of
              elements in fixeddim, and if fixeddim is None, this must
              be too
           margindim : int | arraylike of ints | None
              The index or indices of the physical or photometric
              properties to be maginalized over, numbered in the same
              way as with fixeddim; if set to None, no marginalization
              is performed
           ngrid : int or listlike containing ints
              number of points in each dimension of the output grid;
              if this is an iterable, it must have nphys + nphot -
              len(fixeddim) - len(margindim) elements
           qmin : float | arraylike
              minimum value in the output grid in each quantity; if
              left as None, defaults to the minimum value in the
              library; if this is an iterable, it must contain a
              number of elements equal to nphys + nphot -
              len(fixeddim) - len(margindim)
           qmax : float | arraylike
              maximum value in the output grid in each quantity; if
              left as None, defaults to the maximum value in the
              library; if this is an iterable, it must have the same
              number of elements as qmin
           grid : listlike of arrays
              set of values defining the grid on which the PDF is to
              be evaluated, in the same format used by meshgrid
           norm : bool
              if True, returned pdf's will be normalized to integrate
              to 1
           filters : listlike of strings
              list of photometric filters to use; if left as None, and
              only 1 set of photometric filters has been defined for
              the cluster_slug object, that set will be used by
              default

        Returns:
           grid_out : array
              array of values at which the PDF is evaluated; contents
              are the same as returned by meshgrid
           pdf : array
              array of marginal posterior probabilities at each point
              of the output grid, for each input set of properties; the leading
              dimensions match the leading dimensions produced by
              broadcasting the leading dimensions of photprop and
              photerr together, while the trailing dimensions match
              the dimensions of the output grid
        """
        # Were we given a set of filters?
        if filters is None:

            # No filters given; if we have only a single filter set
            # stored, just use it
            if len(self.__filtersets) == 1:
                return self.__filtersets[0]['bp']. \
                    mpdf_gen(fixeddim, fixedprop, margindim,
                             ngrid=ngrid, qmin=qmin, qmax=qmax,
                             grid=grid, norm=norm)
            else:
                raise ValueError("must specify a filter set")

        else:

            # We were given a filter set; add it if it doesn't exist
            self.add_filters(filters)

            # Find the bp object we should use
            for f in self.__filtersets:
                if f['filters'] == filters:
                    bp = f['bp']
                    break

            # Call the logL method
            return bp.mpdf_gen(fixeddim, fixedprop, margindim,
                               ngrid=ngrid, qmin=qmin, qmax=qmax,
                               grid=grid, norm=norm)
        
        
    def mcmc(self, photprop, photerr=None, mc_walkers=100,
             mc_steps=500, mc_burn_in=50, filters=None):
        """
        This function returns a sample of MCMC walkers for cluster
        mass, age, and extinction 

        Parameters:
           photprop : arraylike, shape (nfilter) or (..., nfilter)
              array giving the photometric values; for a
              multidimensional array, the operation is vectorized over
              the leading dimensions
           photerr : arraylike, shape (nfilter) or (..., nfilter)
              array giving photometric errors; for a multidimensional
              array, the operation is vectorized over the leading
              dimensions
           mc_walkers : int
              number of walkers to use in the MCMC
           mc_steps : int
              number of steps in the MCMC
           mc_burn_in : int
              number of steps to consider "burn-in" and discard
           filters : listlike of strings
              list of photometric filters to use; if left as None, and
              only 1 set of photometric filters has been defined for
              the cluster_slug object, that set will be used by
              default

        Returns
           samples : array
              array of sample points returned by the MCMC
        """

        # Were we given a set of filters?
        if filters is None:

            # No filters given; if we have only a single filter set
            # stored, just use it
            if len(self.__filtersets) == 1:
                return self.__filtersets[0]['bp']. \
                    mcmc(photprop, photerr, mc_walkers, mc_steps, 
                         mc_burn_in)
            else:
                raise ValueError("must specify a filter set")

        else:

            # We were given a filter set; add it if it doesn't exist
            self.add_filters(filters)

            # Find the bp object we should use
            for f in self.__filtersets:
                if f['filters'] == filters:
                    bp = f['bp']
                    break

            # Call the mcmc method
            return bp.mcmc(photprop, photerr, mc_walkers, mc_steps, 
                           mc_burn_in)


    def bestmatch(self, phot, photerr=None, nmatch=1, 
                  bandwidth_units=False, filters=None):
        """
        Searches through the simulation library and returns the closest
        matches to an input set of photometry.

        Parameters:
           phot : arraylike, shape (nfilter) or (..., nfilter)
              array giving the photometric values; for a
              multidimensional array, the operation is vectorized over
              the leading dimensions
           photerr : arraylike, shape (nfilter) or (..., nfilter)
              array giving photometric errors, which must have the
              same shape as phot; if this is not None,
              then distances will be measured in units of the
              photometric error if bandwidth_units is False, or in
              units of the bandwidth added in quadrature with the
              errors if it is True
           nmatch : int
              number of matches to return; returned matches will be
              ordered by distance from the input
           bandwidth_units : bool
              if False, distances are computed based on the
              logarithmic difference in luminosity; if True, they are
              measured in units of the bandwidth
           filters : listlike of strings
              list of photometric filters to use; if left as None, and
              only 1 set of photometric filters has been defined for
              the cluster_slug object, that set will be used by
              default

        Returns:
           matches : array, shape (..., nmatch, nphys + nfilter)
              best matches to the input photometry; shape in the
              leading dimensions will be the same as for phot, and if
              nmatch == 1 then that dimension will be omitted; in the
              final dimension, the first 3 elements give log M, log T,
              and A_V, while the last nfilter give the photometric
              values; if created with use_extinct = False, the A_V
              dimension is omitted
           dist : array, shape (..., nmatch)
              distances between the matches and the input photometry
        """

        # Were we given a set of filters?
        if filters is None:

            # No filters given; if we have only a single filter set
            # stored, just use it
            if len(self.__filtersets) == 1:
                return self.__filtersets[0]['bp']. \
                    bestmatch(phot, photerr, nmatch, bandwidth_units)
            else:
                raise ValueError("must specify a filter set")

        else:

            # We were given a filter set; add it if it doesn't exist
            self.add_filters(filters)

            # Find the bp object we should use
            for f in self.__filtersets:
                if f['filters'] == filters:
                    bp = f['bp']
                    break

            # Call the bestmatch method
            return bp.bestmatch(phot, photerr, nmatch, bandwidth_units)


    def bestmatch_phys(self, phys, nmatch=1, bandwidth_units=False,
                       filters=None):
        """
        Searches through the simulation library and returns the closest
        matches to an input set of photometry.

        Parameters:
           phys : arraylike, shape (nphys) or (..., nphys)
              array giving the physical values; for a
              multidimensional array, the operation is vectorized over
              the leading dimensions
           nmatch : int
              number of matches to return; returned matches will be
              ordered by distance from the input
           bandwidth_units : bool
              if False, distances are computed based on the
              logarithmic difference in physical properties; if True,
              they are measured in units of the bandwidth
           filters : listlike of strings
              list of photometric filters to use; if left as None, and
              only 1 set of photometric filters has been defined for
              the cluster_slug object, that set will be used by
              default

        Returns:
           matches : array, shape (..., nmatch, nphys + nfilter)
              best matches to the input properties; shape in the
              leading dimensions will be the same as for phot, and if
              nmatch == 1 then that dimension will be omitted
           dist : array, shape (..., nmatch)
              distances between the matches and the input physical
              properties
        """
        # Were we given a set of filters?
        if filters is None:

            # No filters given; if we have only a single filter set
            # stored, just use it
            if len(self.__filtersets) == 1:
                return self.__filtersets[0]['bp']. \
                    bestmatch_phys(phys, nmatch, bandwidth_units)
            else:
                raise ValueError("must specify a filter set")

        else:

            # We were given a filter set; add it if it doesn't exist
            self.add_filters(filters)

            # Find the bp object we should use
            for f in self.__filtersets:
                if f['filters'] == filters:
                    bp = f['bp']
                    break

            # Call the bestmatch method
            return bp.bestmatch_phys(phys, nmatch, bandwidth_units)


    def make_approx_phot(self, phys, squeeze=True, filter_ignore=None,
                         filters=None):
        """
        Returns an object that can be used for a fast approximation of
        the PDF of photometric properties that corresponds to a set of
        physical properties. The PDF produced by summing over the
        points returned is guaranteed to account for at least 1-reltol
        of the marginal photometric probability, and to represent the
        shape of the PDF in photometric space within a local accuracy
        of reltol as well.

        Parameters:
           phys : arraylike, shape (nphys) or (N, nphys)
              the set or sets of physical properties for which the
              approximation is to be generated
           squeeze : bool
              if True, the representation returned will be squeezed to
              minimize the number of points included, using reltol as
              the error tolerance
           filter_ignore : None or listlike of bool
              if None, the kernel density representation returned
              covers all filters; otherwise this must be a listlike of
              bool, one entry per filter, with a value of False
              indicating that filter should be excluded from the
              values returned; suppressing filters can allow for more
              efficient representations
           filters : listlike of strings
              list of photometric filters to use; if left as None, and
              only 1 set of photometric filters has been defined for
              the cluster_slug object, that set will be used by
              default

        Returns:
           x : array, shape (M, nphot), or a list of such arrays
              an array containing the list of points to be used for
              the approximation
           wgts : array, shape (M), or a list of such arrays
              an array containing the weights of the points

        Notes:
           if the requested relative tolerance cannot be reached for
           numerical reasons (usually because the input point is too
           far from the library to allow accurate computation), x and
           wgts will be return as None, and a warning will be issued
        """

        # Were we given a set of filters?
        if filters is None:

            # No filters given; if we have only a single filter set
            # stored, just use it
            if len(self.__filtersets) == 1:
                return self.__filtersets[0]['bp']. \
                    make_approx_phot(phys, squeeze=squeeze,
                                     filter_ignore=filter_ignore)
            else:
                raise ValueError("must specify a filter set")

        else:

            # We were given a filter set; add it if it doesn't exist
            self.add_filters(filters)

            # Find the bp object we should use
            for f in self.__filtersets:
                if f['filters'] == filters:
                    bp = f['bp']
                    break

            # Call the method
            return bp.make_approx_phot(phys, squeeze=squeeze,
                                       filter_ignore=filter_ignore)


    def make_approx_phys(self, phot, photerr=None, squeeze=True, 
                         phys_ignore=None, filters=None, tol=None):
        """
        Returns an object that can be used for a fast approximation of
        the PDF of physical properties that corresponds to a set of
        photometric properties. The PDF produced by summing over the
        points returned is guaranteed to account for at least 1-reltol
        of the marginal photometric probability, and to represent the
        shape of the PDF in photometric space within a local accuracy
        of reltol as well.

        Parameters:
           phot : arraylike, shape (nfilter) or (N, nfilter)
              the set or sets of photometric properties for which the
              approximation is to be generated
           photerr : arraylike, shape (nfilter) or (N, nfilter)
              array giving photometric errors; the number of elements
              in the output lists will be the size that results from
              broadcasting together the leading dimensions of phot and
              photerr
           squeeze : bool
              if True, the representation returned will be squeezed to
              minimize the number of points included, using reltol as
              the error tolerance
           phys_ignore : None or listlike of bool
              if None, the kernel density representation returned
              covers all physical properties; otherwise this must be a
              listlike of bool, one entry per physical dimension, with
              a value of False indicating that dimension should be
              excluded from the values returned; suppressing
              dimensions can allow for more efficient representations
           filters : listlike of strings
              list of photometric filters to use; if left as None, and
              only 1 set of photometric filters has been defined for
              the cluster_slug object, that set will be used by
              default
           tol : float
              if set, this tolerance overrides the value of reltol

        Returns:
           x : array, shape (M, nphys), or a list of such arrays
              an array containing the list of points to be used for
              the approximation, where nphys is the number of
              physical dimensions being returned
           wgts : array, shape (M), or a list of such arrays
              an array containing the weights of the points

        Notes:
           if the requested relative tolerance cannot be reached for
           numerical reasons (usually because the input point is too
           far from the library to allow accurate computation), x and
           wgts will be return as None, and a warning will be issued
        """

        # Were we given a set of filters?
        if filters is None:

            # No filters given; if we have only a single filter set
            # stored, just use it
            if len(self.__filtersets) == 1:
                return self.__filtersets[0]['bp']. \
                    make_approx_phys(phot, photerr=photerr,
                                     squeeze=squeeze,
                                     phys_ignore=phys_ignore,
                                     tol=tol)
            else:
                raise ValueError("must specify a filter set")

        else:

            # We were given a filter set; add it if it doesn't exist
            self.add_filters(filters)

            # Find the bp object we should use
            for f in self.__filtersets:
                if f['filters'] == filters:
                    bp = f['bp']
                    break

            # Call the method
            return bp.make_approx_phys(phot, photerr=photerr,
                                       squeeze=squeeze,
                                       phys_ignore=phys_ignore,
                                       tol=tol)


    def squeeze_rep(self, x, wgts, dims=None, filters=None):
        """
        Takes an input array of positions and weights that form a
        kernel density representation and approximates them using
        fewer points, using an error tolerance of reltol

        Parameters:
           x : array, shape (N, ndim)
              an array of points forming a kernel density
              representation; on exit, x will be resized to (M, ndim)
              with M <= N
           wgts : array, shape (N)
              an array of weights for the kernel density
              representation; on exit, wgts will be resized to (M),
              with M <= N
           dims : array, shape (ndim)
              array specifying which dimensions in the kernel density
              representation the coordinates in x correspond to; if
              left as None, they are assumed to correspond to the
              first ndim dimensions in the data set
           filters : listlike of strings
              list of photometric filters to use; if left as None, and
              only 1 set of photometric filters has been defined for
              the cluster_slug object, that set will be used by
              default

        Returns:
           Nothing
        """
        # Were we given a set of filters?
        if filters is None:

            # No filters given; if we have only a single filter set
            # stored, just use it, then return
            if len(self.__filtersets) == 1:
                self.__filtersets[0]['bp']. \
                    squeeze_rep(x, wgts, dims=dims)
                return
            else:
                raise ValueError("must specify a filter set")

        else:

            # We were given a filter set; add it if it doesn't exist
            self.add_filters(filters)

            # Find the bp object we should use
            for f in self.__filtersets:
                if f['filters'] == filters:
                    bp = f['bp']
                    break

            # Call the method and return
            bp.squeeze_rep(x, wgts, dims=dims)
            return


    def mpdf_approx(self, x, wgts, dims='phys', dims_return=None,
                    ngrid=64, qmin='all', qmax='all', grid=None,
                    norm=True, filters=None):
        """
        Returns the marginal posterior PDF computed from a kernel
        density approximation returned by make_approx_phys or
        make_approx_phot. Outputs are computed on a grid of values, in
        the same style as meshgrid.

        Parameters:
           x : array, shape (M, ndim), or a list of such arrays
              array of points retured by make_approx_phot or
              make_approx_phys
           wgts : array, shape (M) or a list of such arrays
              array of weights returned by make_approx_phot or
              make_approx_phys
           dims : 'phys' | 'phot' | arraylike of ints
              dimensions covered by x and wgts; the strings 'phys' or
              'phot' indicate that they cover all physical or
              photometric dimensions, and correspond to the defaults
              returned by make_approx_phys and make_approx_phot,
              respectively; if dims is an array of ints, these specify
              the dimensions covered by x and wgts, where the
              physical dimensions are numbered 0, 1, ... nphys-1, and
              the photometric ones are nphys, nphys+1,
              ... nphys+nphot-1
           dims_return : None or arraylike of ints
              if None, the output PDF has the same dimensions as
              specified in dims; if not, then dims_return must be a
              subset of dims, and a marginal PDF in certain dimensions
              will be generated
           ngrid : int or listlike containing ints
              number of points in each dimension of the output grid;
              if this is an iterable, it must have the same number of
              elements as idx
           qmin : float | listlike | 'zoom' | 'all' 
              minimum value in the output grid in each quantity; if
              this a float, it is applied to each dimension; if it is
              an iterable, it must contain the same number of elements
              as the number of dimensions being returned, as gives the
              minimum in each dimension; if it is 'zoom' or 'all', the
              minimum is chosen automatically, with 'zoom' focusing on
              a region encompassing the probability maximum, and 'all'
              encompassing all the points in the representation
           qmax : float | listlike | 'zoom' | 'all'
              same as qmin, but for the maximum of the output grid
           grid : listlike of arrays
              set of values defining the grid on which the PDF is to
              be evaluated, in the same format used by meshgrid
           norm : bool
              if True, returned pdf's will be normalized to integrate
              to 1
           filters : listlike of strings
              list of photometric filters to use; if left as None, and
              only 1 set of photometric filters has been defined for
              the cluster_slug object, that set will be used by
              default

        Returns:
           grid_out : array
              array of values at which the PDF is evaluated; contents
              are the same as returned by meshgrid
           pdf : array
              array of marginal posterior probabilities at each point
              of the output grid, for each input cluster; the leading
              dimensions match the leading dimensions produced by
              broadcasting the leading dimensions of photprop and
              photerr together, while the trailing dimensions match
              the dimensions of the output grid
        """

        # Were we given a set of filters?
        if filters is None:

            # No filters given; if we have only a single filter set
            # stored, just use it
            if len(self.__filtersets) == 1:
                return self.__filtersets[0]['bp']. \
                    mpdf_approx(x, wgts, dims=dims, 
                                dims_return=dims_return,
                                ngrid=ngrid, qmin=qmin,
                                qmax=qmax, grid=grid, norm=norm)

            else:
                raise ValueError("must specify a filter set")

        else:

            # We were given a filter set; add it if it doesn't exist
            self.add_filters(filters)

            # Find the bp object we should use
            for f in self.__filtersets:
                if f['filters'] == filters:
                    bp = f['bp']
                    break

            # Call the method
            return bp.mpdf_approx(x, wgts, dims=dims, 
                                  dims_return=dims_return,
                                  ngrid=ngrid, qmin=qmin,
                                  qmax=qmax, grid=grid, norm=norm)

    ##################################################################
    # Methods to draw samples
    ##################################################################
    def draw_sample(self,  photerr=None, nsample=1, filters=None):
        """
        Returns a randomly-drawn sample of clusters for a given filter
        set

        Parameters:
           photerr : arraylike, shape (nphot)
              photometric errors to apply to the output photometry;
              these are added in quadrature with the kernel density
              estimation bandwidth
           nsample : int
              number of random samples to draw for each set of
              physical properties; must be positive
           filters : listlike of strings
              list of photometric filters to use; if left as None, and
              only 1 set of photometric filters has been defined for
              the cluster_slug object, that set will be used by
              default

        Returns:
           samples : array, shape (nsample, nphys+nphot)
              random sample drawn from the kernel density object; for
              the final dimension in the output, the first nphys
              elements are the physical quantities, the next nphot are
              the photometric quantities
        """
        
        # Were we given a set of filters?
        if filters is None:

            # No filters given; if we have only a single filter set
            # stored, just use it
            if len(self.__filtersets) == 1:
                return self.__filtersets[0]['bp']. \
                    draw_sample(photerr=photerr, nsample=nsample)

            else:
                raise ValueError("must specify a filter set")

        else:

            # We were given a filter set; add it if it doesn't exist
            self.add_filters(filters)

            # Find the bp object we should use
            for f in self.__filtersets:
                if f['filters'] == filters:
                    bp = f['bp']
                    break

            # Call the method
            return bp.draw_sample(photerr=photerr, nsample=nsample)

    def draw_phot(self, physprop, physidx=None, photerr=None,
                  nsample=1, filters=None):
        """
        Returns a randomly-drawn sample of clusters for a given filter
        set

        Parameters:
           physprop : arraylike
              physical properties to be used; the final dimension of
              the input must have len(physidx) indices, or nphys
              indicates if physidx is None; if the input is a
              multidimensional array, the operation is vectorized over
              the leading dimensions physical properties
           physidx : arraylike
              indices of the physical quantities being constrained; if
              left as None, all physical properties are set, and
              physprop must have a trailing dimension of size equal to
              nphys; otherwise this must be an arraylike of <= nphys
              positive integers, each unique and in the range [0,
              nphys), specying which physical dimensions are
              constrained
           photerr : arraylike, shape (nphot)
              photometric errors to apply to the output photometry;
              these are added in quadrature with the kernel density
              estimation bandwidth
           nsample : int
              number of random samples to draw for each set of
              physical properties; must be positive
           filters : listlike of strings
              list of photometric filters to use; if left as None, and
              only 1 set of photometric filters has been defined for
              the cluster_slug object, that set will be used by
              default

        Returns:
           samples : array, shape (nsample, nphys+nphot)
              random sample drawn from the kernel density object; for
              the final dimension in the output, the first nphys
              elements are the physical quantities, the next nphot are
              the photometric quantities
        """
        
        # Were we given a set of filters?
        if filters is None:

            # No filters given; if we have only a single filter set
            # stored, just use it
            if len(self.__filtersets) == 1:
                return self.__filtersets[0]['bp']. \
                    draw_phot(physprop, physidx=physidx,
                              photerr=photerr, nsample=nsample)

            else:
                raise ValueError("must specify a filter set")

        else:

            # We were given a filter set; add it if it doesn't exist
            self.add_filters(filters)

            # Find the bp object we should use
            for f in self.__filtersets:
                if f['filters'] == filters:
                    bp = f['bp']
                    break

            # Call the method
            return bp.draw_phot(physprop, physidx=physidx,
                                photerr=photerr, nsample=nsample)
