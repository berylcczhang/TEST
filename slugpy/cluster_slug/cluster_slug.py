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

# Import the data reading and Bayesian inference stuff we need
from ..bayesphot import bp
from ..read_cluster_prop import read_cluster_prop
from ..read_cluster_phot import read_cluster_phot

##################################################################
# Define the cluster_slug class                                  #
##################################################################

class cluster_slug(object):
    """
    A class that can be used to estimate the PDF of star cluster
    properties (mass, age, extinction) from a set of input photometry
    in various bands.
    """

    ##################################################################
    # Initializer method
    ##################################################################
    def __init__(self, libname=None, filters=None, bandwidth='auto',
                 ktype='gaussian', priors=None, sample_density=None,
                 reltol=1.0e-3, abstol=1.0e-10, leafsize=16):
        """
        Initialize a cluster_slug object.

        Parameters
           libname : string
              name of the SLUG model to load; if left as None, the default
              is $SLUG_DIR/cluster_slug/CLUSTER_SLUG
           filters : iterable of stringlike
              list of filter names to be used for inferenence
           bandwidth : 'auto' | array, shape (M)
              bandwidth for kernel density estimation; if set to
              'auto', the bandwidth will be estimated automatically
           ktype : string
              type of kernel to be used in densty estimation; allowed
              values are 'gaussian' (default), 'epanechnikov', and
              'tophat'; only Gaussian can be used with error bars
           priors : array, shape (N) | callable | None
              prior probability on each data point; interpretation
              depends on the type passed:
                 array, shape (N) : 
                    values are interpreted as the prior probability of
                    each data point
                 callable : 
                    the callable must take as an argument an array of
                    shape (N, nphys), and return an array of shape (N)
                    giving the prior probability at each data point
                  None :
                    all data points have equal prior probability
           sample_density : array, shape (N) | callable | 'auto' | None
              the density of the data samples at each data point; this
              need not match the prior density; interpretation depends
              on the type passed:
                 array, shape (N) : 
                    values are interpreted as the density of data
                    sampling at each sample point
                 callable : 
                    the callable must take as an argument an array of
                    shape (N, nphys), and return an array of shape (N)
                    giving the sampling density at each point
                 'auto' :
                    the sample density will be computed directly from
                    the data set; note that this can be quite slow for
                    large data sets, so it is preferable to specify
                    this analytically if it is known
                  None :
                    data are assumed to be uniformly sampled
           reltol : float
              relative error tolerance; errors on all returned
              probabilities p will satisfy either
              |p_est - p_true| <= reltol * p_est   OR
              |p_est - p_true| <= abstol,
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

        # Load the cluster data
        if libname is None:
            self.__libname = osp.join('cluster_slug', 'CLUSTERSLUG')
            if 'SLUG_DIR' in os.environ:
                self.__libname = osp.join(os.environ['SLUG_DIR'], 
                                          self.__libname)
        else:
            self.libname = libname
        prop = read_cluster_prop(self.__libname)
        phot = read_cluster_phot(self.__libname)

        # Record available filters
        self.__allfilters = phot.filter_names

        # Record other stuff that we'll use to construct bp objects
        # later
        self.__bandwidth = bandwidth
        self.__ktype = ktype
        self.__priors = priors
        self.__sample_density = sample_density
        self.__reltol = reltol
        self.__abstol = abstol

        # Only keep cases where the luminosity is non-zero; zero
        # luminosity cases correspond to clusters where the most
        # massive star produced was below the minimum mass included in
        # our evolutionary tracks
        idx = np.where(np.amin(phot.phot, axis=1) > 0)[0]

        # Build dataset array; 1st 3 dimensions are log mass, log age,
        # and A_V remaining ones are photometric values
        self.__ds = np.zeros((len(idx),
                                   3+len(phot.filter_names)))
        self.__ds[:,0] = np.log10(prop.actual_mass[idx])
        self.__ds[:,1] = np.log10(prop.time[idx] - prop.form_time[idx])
        self.__ds[:,2] = prop.A_V[idx]
        self.__ds[:,3:] = phot.phot_ex[idx,:]

        # Fill in any cases where extincted photometry is not
        # available
        for i in np.where(np.isnan(phot.phot_ex[idx[0],:]))[0]:
            warnstr = "cluster_slug: extincted photometry " + \
                      "unavailable for filter " + \
                      phot.filter_names[i] + \
                      ", probably because extinction curve " + \
                      "does not cover the required wavelength " + \
                      "range; using unextincted values instead"
            warnings.warn(warnstr)
            self.__ds[:,3+i] = phot.phot[idx,i]

        # Take log of photometric values if they are recorded in a
        # linear system
        for i, f in enumerate(phot.filter_units):
            if 'mag' not in f:
                self.__ds[:,3+i] = np.log10(self.__ds[:,3+i])

        # Initialize empty dict containing filter sets
        self.__filtersets = []

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
        Parameters
           filters : iterable of stringlike
              list of filter names to be used for inferenence

        Returns
           nothing
        """

        # If we already have this filter set in our dict, do nothing
        for f in self.__filtersets:
            if filters == f['filters']:
                return

        # We're adding a new filter set, so save its name
        newfilter = { 'filters' : filters}

        # Construct data set to use with this filter combination, and
        # fill it in
        newfilter['dataset'] = np.zeros((self.__ds.shape[0], 
                                         3+len(filters)))
        newfilter['dataset'][:,:3] = self.__ds[:,:3]
        for i, f in enumerate(filters):
            if f not in self.__allfilters:
                raise ValueError("unknown filter "+str(f))
            idx = self.__allfilters.index(f)
            newfilter['dataset'][:,3+i] = self.__ds[:,3+idx]

        # Build a bp object to go with this data set
        newfilter['bp'] = bp(newfilter['dataset'], 3,
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
        self.__priors = pr
        for f in self.__filterset:
            f['bp'].priors = self.__priors


    ##################################################################
    # Define the bandwidth property. This just wraps around the
    # corresponding property defined for bp objects.
    ##################################################################
    @property
    def bandwidth(self):
        return self.__bandwidth

    @priors.setter
    def bandwidth(self, bandwidth):
        self.__bandwidth = bandwidth
        for f in self.__filterset:
            f['bp'].bandwidth = bandwidth


    ##################################################################
    # Wrappers around the bp logL, mpdf, and mcmc functions
    ##################################################################
    def logL(self, physprop, photprop, photerr=None, filters=None):
        """
        This function returns the natural log of the likelihood
        function evaluated at a particular log mass, log age,
        extinction, and set of log luminosities

        Parameters:
           physprop : arraylike, shape (3) or (..., 3)
              array giving values of the log M, log T, and A_V; for a
              multidimensional array, the operation is vectorized over
              the leading dimensions
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
                    logL(physprop, photprop, photerr)
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
            return bp.logL(physprop, photprop, photerr)


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
              computed; 0 = log M, 1 = log T, 2 = A_V; if this is an
              iterable, the joint distribution of the indicated
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

           (only returned if the c library was compiled in DIAGNOSTIC mode)
           nodecheck : int or array, shape(N)
              Number of nodes examined during the evaluation
           leafcheck : int or array, shape(N)
              Number of leaves examined during the evaluation
           termcheck : int or array, shape(N)
              Number of nodes examined during the evaluation for which
              no children were examined
        """

        # Were we given a set of filters?
        if filters is None:

            # No filters given; if we have only a single filter set
            # stored, just use it
            if len(self.__filtersets) == 1:
                return self.__filtersets[0]['bp']. \
                    mpdf(idx, photprop, photerr, ngrid,
                         qmin, qmax, grid, norm)
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
            return bp.mpdf(idx, photprop, photerr, ngrid,
                           qmin, qmax, grid, norm)


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

            # Call the logL method
            return bp.mcmc(photprop, photerr, mc_walkers, mc_steps, 
                           mc_burn_in)
