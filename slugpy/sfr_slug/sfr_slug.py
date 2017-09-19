"""
This defines a class that can be used to estimate the PDF of true star
formation rate from a set of input point mass estimates of the star
formation rate.
"""

import numpy as np
import copy
import os
import os.path as osp
import warnings
import errno
try:
    # Python 3
    from urllib.request import urlopen
except ImportError:
    # Python 2
    from urllib2 import urlopen

# Import the data reading and Bayesian inference stuff we need
from ..bayesphot import bp
from ..read_integrated_prop import read_integrated_prop
from ..read_integrated_phot import read_integrated_phot

# Functions for pre-defined prior probability distributions on the SFR
def dndlogsfr_flat(logsfr):
    return 1.0+logsfr*0.0

def dndlogsfr_schechter(logsfr):
    sfr = np.exp(logsfr)
    return sfr**(-0.51) * np.exp(-sfr/9.2)


# Main sfr_slug class
class sfr_slug(object):
    """
    A class that can be used to estimate the PDF of true star
    formation rate from a set of input point mass estimates of the
    star formation rate.

    Properties
       priors : array, shape (N) | callable | 'flat' | 'schechter' | None
          prior probability on each data point; interpretation
          depends on the type passed; array, shape (N): values are
          interpreted as the prior probability of each data point;
          callable: the callable must take as an argument an array
          of shape (N, nphys), and return an array of shape (N)
          giving the prior probability at each data point; None:
          all data points have equal prior probability; the values
          'flat' and 'schechter' use priors p(log SFR) ~ constant and
          p(log SFR) ~ SFR^alpha exp(-SFR/SFR_*), respectively, where
          alpha = -0.51 and SFR_* = 9.2 Msun/yr are the values
          measured by Bothwell et al. (2011)
       bandwidth : 'auto' | array, shape (M)
          bandwidth for kernel density estimation; if set to
          'auto', the bandwidth will be estimated automatically
    """

    ##################################################################
    # Initializer method
    ##################################################################
    def __init__(self, libname=None, detname=None, filters=None, 
                 bandwidth=0.1, ktype='gaussian', priors=None, 
                 sample_density='read', reltol=1.0e-3,
                 abstol=1.0e-10, leafsize=16):
        """
        Initialize an sfr_slug object.

       Parameters
           libname : string
              name of the SLUG model to load; if left as None, the default
              is $SLUG_DIR/sfr_slug/SFR_SLUG
           detname : string
              name of a SLUG model run with the same parameters but no
              stochasticity; used to establish the non-stochastic
              photometry to SFR conversions; if left as None, the default
              is libname_DET
           filters : iterable of stringlike
              list of filter names to be used for inferenence
           bandwidth : 'auto' | float | array, shape (M)
              bandwidth for kernel density estimation; if set to
              'auto', the bandwidth will be estimated automatically;
              if set to a float, the same bandwidth is used in all
              dimensions
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
           sample_density : array, shape (N) | callable | 'auto' | 'read' | None
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
              this analytically if it is known; 'read': the sample
              density is to be read from a numpy save file whose name
              matches that of the library, with the extension _density.npy
              added; None: data are assumed to be uniformly sampled
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

        Returns
           Nothing

        Raises
           IOError, if the library cannot be found
        """

        # If using the default library, assign the library name
        if libname is None:
            self.__libname = osp.join('sfr_slug', 'SFR_SLUG')
            if 'SLUG_DIR' in os.environ:
                self.__libname = osp.join(os.environ['SLUG_DIR'], 
                                          self.__libname)
        else:
            self.__libname = libname

        # Load the data
        try:
            prop = read_integrated_prop(self.__libname)
            phot = read_integrated_phot(self.__libname)
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
                              "failed to read default sfr_slug " +
                              "library sfr_slug/SFR_SLUG due to missing " +
                              "astropy.io.fits; install astropy " +
                              "or specify a library in a non-FITS " +
                              "format")

            # If we're here, we couldn't open the default library, and
            # it's not because we don't have FITS capability. The file
            # must not exist, or must be damaged. Check if we're in
            # interactive mode. If not, just raise an error and
            # suggest the user to go get the library file.
            errstr = "Unable to open default sfr_slug " + \
                     "library file sfr_slug/SFR_SLUG."
            import __main__ as main
            if hasattr(main, '__file__'):
                # We're not interactive; just raise an error
                raise IOError(errno.EIO,
                              errstr + " " +
                              "Try downloading it from " +
                              "https://sites.google.com/site/runslug/data")

            # If we're here, we don't have hte library file, but we
            # are in interactive mode. Thus offer the user an option
            # to go download the file now.
            usr_response \
                = raw_input(errstr + " Would you like to download it "
                            "now (total size 160 MB)? [y/n] ").\
                lower().strip()
            if not usr_response in ['yes', 'y', 'ye']:
                # User didn't say yes, so raise error
                raise IOError(errno.EIO, "Unable to proceeed")

            # If we're here, download the files
            print("Fetching SFR_SLUG_integrated_prop.fits " +
                  "(this may take a while)...")
            url = urlopen(
                'https://dl.dropboxusercontent.com/s/7la1b5h986rdz29/SFR_SLUG_integrated_prop.fits')
            rawdata = url.read()
            url.close()
            fp = open(osp.join(osp.dirname(self.__libname),
                               'SFR_SLUG_integrated_prop.fits'), 'wb')
            fp.write(rawdata)
            fp.close()
            print("Fetching SFR_SLUG_integrated_phot.fits " +
                  "(this make take a while)...")
            url = urlopen(
                'https://dl.dropboxusercontent.com/s/ra8qf5raqutcf50/SFR_SLUG_integrated_phot.fits')
            rawdata = url.read()
            url.close()
            fp = open(osp.join(osp.dirname(self.__libname),
                               'SFR_SLUG_integrated_phot.fits'), 'wb')
            fp.write(rawdata)
            fp.close()

            # Now try reading the data
            try:
                prop = read_integrated_prop(self.__libname)
                phot = read_integrated_phot(self.__libname)
            except IOError:
                raise IOError(errno.EIO,
                              "still unable to open default library")

        # Load the determinstic run
        if detname is None:
            self.__detname = self.__libname + '_DET'
        else:
            self.__detname = detname
        propdet = read_integrated_prop(self.__detname)
        photdet = read_integrated_phot(self.__detname)

        # If we have been told to read the sample density, do so
        try:
            if sample_density == 'read':
                self.__sample_density = np.load(self.__libname + 
                                                '_density.npy')
        except IOError:
            warnings.warn("unable to load requested sample "+
                          "density file " + 
                          self.__libname + "_density.npy; setting " +
                          "sample_density = 'auto' instead")
            sample_density = 'auto'

        # Store filters
        self.__allfilters = phot.filter_names

        # Get conversions from photometry to SFR for the deterministic run
        sfrdet = propdet.target_mass[0,0]/propdet.time[0]
        self.__conversions = photdet.phot[:,0,0] / sfrdet

        # Deduce SFR from target mass and time
        logsfr = np.log10(
            np.transpose(np.transpose(prop.target_mass)/prop.time[0]))

        # Convert the photometry from the stochastic runs to SFRs
        # estimated using the point mass approximation
        sfrphot = np.transpose(np.transpose(phot.phot) / 
                               self.__conversions)
        sfrphotclip = np.clip(sfrphot, 1e-100, np.inf)
        logsfrphot = np.log10(sfrphotclip)

        # Build array of log SFR, log SFR_phot
        self.__ds = np.zeros((logsfrphot.shape[2],
                                 logsfrphot.shape[0]*
                                 logsfrphot.shape[1]+1))
        self.__ds[:,0] = logsfr
        self.__ds[:,1:] = np.transpose(logsfrphot.reshape(
            (logsfrphot.shape[0]*logsfrphot.shape[1],
             logsfrphot.shape[2])))

        # Record other stuff that we'll use to construct bp objects
        # later
        if type(bandwidth) is not np.ndarray:
            bandwidth = np.array([bandwidth]*(1+len(phot.filter_names)))
        self.__bandwidth = bandwidth
        self.__ktype = ktype
        self.__reltol = reltol
        self.__abstol = abstol
        if sample_density != 'read':
            self.__sample_density = sample_density

        # Initialize empty dict containing filter sets
        self.__filtersets = []

        # Set priors
        self.priors = priors

        # If we have been given a filter list, create the data set to
        # go with it
        if filters is not None:
            self.add_filters(filters)

    ##################################################################
    # Method to return list of available filters
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


    ##################################################################
    # Method to prepare to analyze a particular set of filters
    ##################################################################
    def add_filters(self, filters):
        """
        Add a set of filters to use for cluster property estimation

        Parameters
           filters : iterable of stringlike
              list of filter names to be used for inferenence

        Returns
           nothing
        """

        # Handle the case where we're given just a string for one filter
        if type(filters) is str:
            filters = [filters]

        # If we already have this filter set in our dict, do nothing
        for f in self.__filtersets:
            if filters == f['filters']:
                return

        # We're adding a new filter set, so save its name
        newfilter = { 'filters' : filters}

        # Construct data set to use with this filter combination, and
        # fill it in
        newfilter['dataset'] = np.zeros((self.__ds.shape[0], 
                                         1+len(filters)))
        newfilter['dataset'][:,:1] = self.__ds[:,:1]
        newfilter['idx'] = np.zeros(1+len(filters), dtype=int)
        newfilter['idx'][0] = 0
        for i, f in enumerate(filters):
            if f not in self.__allfilters:
                raise ValueError("unknown filter "+str(f))
            idx = self.__allfilters.index(f)
            newfilter['dataset'][:,1+i] = self.__ds[:,1+idx]
            newfilter['idx'][i+1] = 1+idx

        # Build a bp object to go with this data set
        newfilter['bp'] = bp(newfilter['dataset'], 1,
                             filters=filters,
                             bandwidth = self.__bandwidth,
                             ktype = self.__ktype,
                             priors = self.__priors,
                             sample_density = self.__sample_density,
                             reltol = self.__reltol,
                             abstol = self.__abstol)

        # Save to the master filter list
        self.__filtersets.append(newfilter)


    ##################################################################
    # Define the priors property. This just wraps around the
    # corresponding property defined for bp objects.
    ##################################################################
    @property
    def priors(self):
        return self.__priors

    @priors.setter
    def priors(self, pr):
        if pr == 'flat':
            self.priors = dndlogsfr_flat
        elif pr == 'schechter':
            self.priors = dndlogsfr_schechter
        else:
            self.__priors = pr
            for f in self.__filtersets:
                f['bp'].priors = self.__priors


    ##################################################################
    # Define the bandwidth property. This just wraps around the
    # corresponding property defined for bp objects.
    ##################################################################
    @property
    def bandwidth(self):
        return self.__bandwidth

    @bandwidth.setter
    def bandwidth(self, bandwidth):
        self.__bandwidth = bandwidth
        for f in self.__filtersets:
            f['bp'].bandwidth = np.array(bandwidth)[f['idx']]



    ##################################################################
    # Wrappers around the bp logL, mpdf, and mcmc functions
    ##################################################################
    def logL(self, logSFR, logSFRphot, logSFRphoterr=None, 
             filters=None):
        """
        This function returns the natural log of the likelihood
        function evaluated at a particular log SFR and set of log
        luminosities

        Parameters:
           logSFR : float or arraylike
              float or array giving values of the log SFR; for an
              array, the operation is vectorized
           logSFRphot : float or arraylike, shape (nfilter) or (..., nfilter)
              float or array giving the SFR inferred from photometry using a
              deterministic conversion; for an array, the operation is
              vectorized over the leading dimensions
           logSFRphoterr : float arraylike, shape (nfilter) or (..., nfilter)
              float or array giving photometric SFR errors; for a
              multidimensional array, the operation is vectorized over
              the leading dimensions
           filters : listlike of strings
              list of photometric filters used for the SFR estimation;
              if left as None, and only 1 set of photometric filters
              has been defined for the sfr_slug object, that set will
              be used by default

        Returns:
           logL : float or arraylike
              natural log of the likelihood function
        """

        # If we were given a floats, make them an arrays of shape 1 for
        # the purposes of communicating with c
        if not type(logSFR) is np.ndarray:
            physprop = np.array(logSFR)
            if physprop.ndim == 0:
                physprop = np.reshape(physprop, (1,))
        else:
            physprop = logSFR
        if not type(logSFRphot) is np.ndarray:
            photprop = np.array(logSFRphot)
            if photprop.ndim == 0:
                photprop = np.reshape(photprop, (1,))
        else:
            photprop = logSFRphot
        if logSFRphoterr is not None and \
           type(logSFRphoterr) is not np.ndarray:
            photerr = np.array(logSFRphoterr)
            if photerr.ndim == 0:
                photerr = np.reshape(photerr, (1,))
        else:
            photerr = logSFRphoterr

        # Were we given a set of filters?
        if filters is None:

            # No filters given; if we have only a single filter set
            # stored, just use it
            if len(self.__filtersets) == 1:
                return self.__filtersets[0]['bp']. \
                    logL(physprop, photprop, photerr)
            else:
                raise ValueError("must specify a filter set")

        else:

            # We were given a filter set; add it if it doesn't exist
            self.add_filters(filters)

            # Handle the case where we're given just a string for one filter
            if type(filters) is str:
                filters = [filters]

            # Find the bp object we should use
            for f in self.__filtersets:
                if f['filters'] == filters:
                    bp = f['bp']
                    break

            # Call the logL method
            return bp.logL(physprop, photprop, photerr)


    def mpdf(self, logSFRphot, logSFRphoterr=None, ngrid=128,
             qmin=None, qmax=None, grid=None, norm=True,
             filters=None):
        """
        Returns the marginal probability of log SFR for one or more
        input sets of photometric properties. Output quantities are
        computed on a grid of values, in the same style as meshgrid

        Parameters:
           logSFRphot : float or arraylike
              float or array giving the log SFR inferred from
              photometry using a deterministic conversion; if the
              argument is an array, the operation is vectorized over
              it
           logSFRphoterr : arraylike, shape (nfilter) or (..., nfilter)
              array giving photometric errors; for a multidimensional
              array, the operation is vectorized over the leading
              dimensions
           ngrid : int
              number of points in the output log SFR grid
           qmin : float
              minimum value in the output log SFR grid
           qmax : float
              maximum value in the output log SFR grid
           grid : array
              set of values defining the grid of SFR values at which
              to evaluate; if set, overrides ngrid, qmin, and qmax
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
              array of log SFR values at which the PDF is evaluated
           pdf : array
              array of marginal posterior probabilities at each point
              of the output grid, for each input photometric value;
              the leading dimensions match the leading dimensions
              produced by broadcasting the leading dimensions of
              photprop and photerr together, while the trailing
              dimensions match the dimensions of the output grid
        """

        # If we were given a floats, make them an arrays of shape 1 for
        # the purposes of communicating with c
        if not type(logSFRphot) is np.ndarray:
            photprop = np.array(logSFRphot)
            if photprop.ndim == 0:
                photprop = np.reshape(photprop, (1,))
        else:
            photprop = logSFRphot
        if logSFRphoterr is not None and \
           type(logSFRphoterr) is not np.ndarray:
            photerr = np.array(logSFRphoterr)
            if photerr.ndim == 0:
                photerr = np.reshape(photerr, (1,))
        else:
            photerr = logSFRphoterr

        # Were we given a set of filters?
        if filters is None:

            # No filters given; if we have only a single filter set
            # stored, just use it
            if len(self.__filtersets) == 1:
                return self.__filtersets[0]['bp']. \
                    mpdf(0, photprop, photerr, ngrid,
                         qmin, qmax, grid, norm)
            else:
                raise ValueError("must specify a filter set")

        else:

            # We were given a filter set; add it if it doesn't exist
            self.add_filters(filters)

            # Handle the case where we're given just a string for one filter
            if type(filters) is str:
                filters = [filters]

            # Find the bp object we should use
            for f in self.__filtersets:
                if f['filters'] == filters:
                    bp = f['bp']
                    break

            # Call the logL method
            return bp.mpdf(0, photprop, photerr, ngrid,
                           qmin, qmax, grid, norm)


    def mcmc(self, photprop, photerr=None, mc_walkers=100,
             mc_steps=500, mc_burn_in=50, filters=None):
        """
        This function returns a sample of MCMC walkers for log SFR

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

        # If we were given a floats, make them an arrays of shape 1 for
        # the purposes of communicating with c
        if not type(logSFRphot) is np.ndarray:
            photprop = np.array([logSFRphot])
            if photprop.ndim == 0:
                photprop = np.reshape(photprop, (1,))
        else:
            photprop = logSFRphot
        if logSFRphoterr is not None and \
           type(logSFRphoterr) is not np.ndarray:
            photerr = np.array([logSFRphoterr])
            if photerr.ndim == 0:
                photerr = np.reshape(photerr, (1,))
        else:
            photerr = logSFRphoterr
              
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

            # Handle the case where we're given just a string for one filter
            if type(filters) is str:
                filters = [filters]

            # Call the logL method
            return bp.mcmc(photprop, photerr, mc_walkers, mc_steps, 
                           mc_burn_in)

