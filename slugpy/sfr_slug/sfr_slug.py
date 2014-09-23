"""
This defines a class that can be used to estimate the PDF of true star
formation rate from a set of input point mass estimates of the star
formation rate.
"""

from sklearn.neighbors.kde import KernelDensity
import numpy as np
import os
import os.path as osp
from copy import deepcopy
from collections import namedtuple
from ..read_integrated_phot import read_integrated_phot
from ..read_integrated_prop import read_integrated_prop
from gkcons import gk

class sfr_slug(object):
    """
    A class that can be used to estimate the PDF of true star
    formation rate from a set of input point mass estimates of the
    star formation rate.

    Attributes
       dataset : array
          the training dataset to be used for KDE estimation
       dataset_filters : list of string
          filters represented in the training data set
       conversions : array
          conversions between luminosity and SFR in the point mass
          estimate
       kde : KernelDensity object
          a KernelDensity estimator constructed from the dataset
       kde_filters : list of string
          filters represented in the KDE object
       kde_limits : array
          range over which KDE is non-zero
       libname : string
          name of the SLUG model from which the data set was read
       detname : string
          name of a SLUG model run with the same parameters as
          libname, but with no stochasticity

    Methods
       __init__ : initializer
       __call__ : returns the PDF of SFR from input point-mass SFR
                  estimates
       pdf : returns the joint PDF of the data set in at a specified
             set of data points
       pdfgrid : returns the joint PDF of the data set evaluated on a
                 grid
    """

    ##################################################################
    # Initializer method
    ##################################################################
    def __init__(self, libname=None, detname=None, bandwidth=0.1):
        """
        Initialize an sfr_pdf_est object.

        Parameters
           libname : string
              name of the SLUG model to load; if left as None, the default
              is $SLUG_DIR/sfr_slug/SFR_SLUG
           detname : string
              name of a SLUG model run with the same parameters but no
              stochasticity; used to establish the non-stochastic
              photometry to SFR conversions; if left as None, the default
              is libname_DET
           bandwidth : float
             bandwidth of the kernel to use in density estimates

        Returns
           Nothing

        Raises
           IOError, if the library cannot be found
        """

        # Load the data
        if libname is None:
            self.libname = osp.join('sfr_slug', 'SFR_SLUG')
            if 'SLUG_DIR' in os.environ:
                self.libname = osp.join(os.environ['SLUG_DIR'], 
                                        self.libname)
        else:
            self.libname = libname
        prop = read_integrated_prop(self.libname)
        phot = read_integrated_phot(self.libname)

        # Load the determinstic run
        if detname is None:
            self.detname = self.libname + '_DET'
        else:
            self.detname = detname
        propdet = read_integrated_prop(self.detname)
        photdet = read_integrated_phot(self.detname)

        # Store filters
        self.dataset_filters = ['SFR'] + phot.filter_names

        # Get conversions from photometry to SFR for the deterministic run
        sfrdet = propdet.target_mass[0,0]/propdet.time[0]
        self.conversions = photdet.phot[:,0,0] / sfrdet

        # Deduce SFR from target mass and time
        logsfr = np.log10(
            np.transpose(np.transpose(prop.target_mass)/prop.time[0]))

        # Convert the photometry from the stochastic runs to SFRs
        # estimated using the point mass approximation
        sfrphot = np.transpose(np.transpose(phot.phot) / 
                               self.conversions)
        sfrphotclip = np.clip(sfrphot, 1e-100, np.inf)
        logsfrphot = np.log10(sfrphotclip)

        # Build array of log SFR, log SFR_phot
        self.dataset = np.zeros((logsfrphot.shape[2],
                                 logsfrphot.shape[0]*
                                 logsfrphot.shape[1]+1))
        self.dataset[:,0] = logsfr
        self.dataset[:,1:] = np.transpose(logsfrphot.reshape(
            (logsfrphot.shape[0]*logsfrphot.shape[1],
             logsfrphot.shape[2])))

        # Purge data set of runs where SFR came out to exactly 0
        idx = np.where(self.dataset[:,2].flatten() > -99.99)[0]
        self.dataset = self.dataset[idx,:]

        # Record the bandwidth
        self.bandwidth = bandwidth

        # Initialize to no KDE filters set
        self.kde_filters = None



    ##################################################################
    # Method to set up a KDE
    ##################################################################
    def set_kde(self, filter_name):
        """
        Builds a KernelDensity object to estimate the joint PDF of the
        dataset in one or more filters

        Parameters
           filter_name : string or iterable
              name or names of filters used for photometric SFR
              estimates, or 'SFR' for the true SFR; if left as None,
              the currently-stored filters are used

        Returns
           Nothing

        Raises
           ValueError, if filter_name is not one of the filters
           available in dataset_filters
        """

        # Check if we already have a kde build for this set of
        # filters; if not, build it
        if filter_name != self.kde_filters:

            # Figure out which parts of the data set to use
            idx = []
            if hasattr(filter_name, '__iter__'):
                # List of filters given
                for f in filter_name:
                    if not f in self.dataset_filters:
                        raise ValueError("filter "+str(f) + 
                                         " not available")
                    idx.append(self.dataset_filters.index(f))
            else:
                if not filter_name in self.dataset_filters:
                    raise ValueError("filter "+str(filter_name) + 
                                     " not available")
                idx.append(self.dataset_filters.index(filter_name))

            # Grab the parts of the data set we want, and record the
            # filters they correspond to
            data_sub = self.dataset[:,idx]
            self.kde_filters = filter_name

            # Construct the KDE object
            self.kde \
                = KernelDensity(bandwidth=self.bandwidth, 
                                kernel='epanechnikov').fit(data_sub)

            # Compute the limits on the KDE
            self.kde_limits = np.zeros((len(idx),2))
            for i in range(len(idx)):
                self.kde_limits[i,0] \
                    = np.amin(data_sub[:,i]) - self.bandwidth
                self.kde_limits[i,1] \
                    = np.amax(data_sub[:,i]) + self.bandwidth



    ##################################################################
    # Method to compute a KDE
    ##################################################################
    def get_kde(self, filter_name):
        """
        Builds a KernelDensity object to estimate the joint PDF of the
        dataset in one or more filters, then returns it. Same as
        set_kde, but the kde is returned rather than being stored
        interally.

        Parameters
           filter_name : string or iterable
              name or names of filters; the string 'SFR' corresponds
              to the true star formation rate

        Returns
           kde : KernelDensity object
              the computed KernelDensity object
           kde_limts : array, shape (N_dim, 2)
              range over which the kde is non-zero; element [:,0]
              gives the minimum in each dimension, and [:,1] gives the
              maximum

        Raises
           ValueError, if filter_name is not one of the filters
           available in dataset_filters
        """

        # Check if we already have a kde build for this set of
        # filters; if not, build it
        if filter_name == self.kde_filters:

            return deepcopy(self.kde)

        else:

            # Figure out which parts of the data set to use
            idx = []
            if hasattr(filter_name, '__iter__'):
                # List of filters given
                for f in filter_name:
                    if not f in self.dataset_filters:
                        raise ValueError("filter "+str(f) + 
                                         " not available")
                    idx.append(self.dataset_filters.index(f))
                idx = np.array(idx)
            else:
                if not filter_name in self.dataset_filters:
                    raise ValueError("filter "+str(filter_name) + 
                                     " not available")
                idx.append(self.dataset_filters.index(filter_name))

            # Grab the parts of the data set we want
            data_sub = self.dataset[:,idx]

            # Construct the KDE object
            kde = \
                  KernelDensity(bandwidth=self.bandwidth, 
                                kernel='epanechnikov').fit(data_sub)

            # Compute the limits on the KDE
            kde_limits = np.zeros((len(idx),2))
            for i in range(len(idx)):
                kde_limits[i,0] \
                    = np.amin(data_sub[:,i]) - self.bandwidth
                kde_limits[i,1] \
                    = np.amax(data_sub[:,i]) + self.bandwidth

            # Return
            return kde, kde_limits



    ##################################################################
    # Method to compute the joint PDF in one or more filters
    ##################################################################
    def pdf(self, x, filter_name=None, nosave=False):
        """
        Return the PDF of the data set in one or more filters.

        Parameters
           x : array
              point or points at which the PDF is to be evaluated; the
              trailing dimension of x must match the number of
              elements in filter_name
           filter_name : string or iterable
              name or names of filters used for photometric SFR
              estimates, or 'SFR' for the true SFR; if left as None,
              the currently-stored filters are used
           nosave : bool
              if True, the KDE constructed as part of the computation
              is not saved; if False, it is saved and overwrites the
              existing KDE

        Returns
           logpdf : array
              log of the value of the PDF evaluated at each of the
              input points

        Raises
           ValueError, if filter_name is left as None and no filters
           are set
        """

        # Build the KDE if needed
        if filter_name is None:
            if self.kde_filters is None:
                raise ValueError("must set filter_name if none "+
                                 "are stored")
            kde = self.kde
            kde_limits = self.kde_limits
        else:
            if nosave:
                kde, kde_limits = self.get_kde(filter_name)
            else:
                self.set_kde(filter_name)
                kde = self.kde
                kde_limits = self.kde_limits

        # Evalute the PDF
        pdf = kde.score_samples(x)

        # Return
        return pdf


    ##################################################################
    # Method to compute the joint PDF in one or more filters on a grid
    ##################################################################
    def pdfgrid(self, filter_name=None, nmesh=50, lim=None,
                nosave=False):
        """
        Return the PDF of the data set in one or more filters,
        evaluated on a uniformly-spaced grid over the data.

        Parameters
           filter_name : string or iterable
              name or names of filters used for photometric SFR
              estimates; if left as None, the currently-stored filters
              are used
           nmesh : int
              number of sample points per dimension to use in
              constructing a sampling grid
           lim : array, shape (Ndim, 2)
              limits of the sampling grid; element [i,0] gives the
              lower limit in the dimension corresponding to
              filter_name[i], and element [i,1] gives the upper
              limit; if left as None, limits are chosen automatically
              to fit the data set
           nosave : bool
              if True, the KDE constructed as part of the computation
              is not saved; if False, it is saved and overwrites the
              existing KDE

        Returns
           xlim : array, shape (Ndim, 2)
              limits of the grid over which the data was evaluated;
              Ndim is the number of dimensions of the data, which is
              equal to the number of elements in filter_name
           logpdf : ndarray
              log of the value of the PDF evaluated at each of the
              input points

        Raises
           ValueError, if filter_name is left as None and no filters
           are set
        """

        # Build the KDE if needed
        if filter_name is None:
            if self.kde_filters is None:
                raise ValueError("must set filter_name if none "+
                                 "are stored")
            kde = self.kde
            kde_limits = self.kde_limits
        else:
            if nosave:
                kde, kde_limits = self.get_kde(filter_name)
            else:
                self.set_kde(filter_name)
                kde = self.kde
                kde_limits = self.kde_limits

        # Build the grid of points at which to evaluate the data
        ndim = kde_limits.shape[0]
        if lim is None:
            lim = np.copy(kde_limits)
        if ndim == 1:
            x = np.linspace(lim[0,0],
                            lim[0,1],
                            nmesh)[:, np.newaxis]
        else:
            vec = []
            for i in range(ndim):
                vec.append(np.linspace(lim[i,0], lim[i,1],
                                       nmesh))
            mg = np.meshgrid(*vec)
            x = mg[0].flatten()[:,np.newaxis]
            for i in range(1,ndim):
                x = np.append(x, mg[i].flatten()[:,np.newaxis], 
                              axis=1)

        # Evalute the PDF
        pdf = kde.score_samples(x)

        # Reshape the PDF into the desired shape
        pdf = pdf.reshape((nmesh,)*ndim)

        # Return
        return lim, pdf



    ##################################################################
    # Method to compute PDF of SFR given a point mass estimate
    ##################################################################
    def __call__(self, logsfr_est, logsfr_err=None, filter_name=None,
                 nmesh=100, logsfr_lim=None, prior=None,
                 error_est=False, gkorder='61'):
        """
        Return an estimate of the PDF of true star formation rate for
        one or more point mass estimates of the SFR using a particular
        photometric filter.

        Parameters
           logsfr_est : float or array
              a point mass estimate of log_10 SFR; can be a float, 1D
              array or 2D array. For a 1D array the data are taken to
              be a series of point mass estimates of log_10 SFR. For a
              2D array, the trailing dimension must match the length
              of filter_name, and the data are taken to represent one
              or more log_10 SFR point mass estimates using different
              filters
           logsfr_err : float or array
              error on logsfr_est; must be the same shape as
              logsfr_est; if left as None, data are assumed to have
              negligible errors and are treated as delta functions
           filter_name : string or iterable
              name or names of filters used for photometric SFR
              estimates; if left as None, stored values are used
           nmesh : int
              number of mesh points at which to evaluate the SFR; note
              that accurate normalization requires that this not be
              too small
           logsfr_lim : array, shape (2)
              limit the log SFR values considered to lie in the range
              logsfr_lim[0] to logsfr_lim[1]
           prior : callable
              a callable that returns the prior probability
              distribution of log SFR; if set to None, the prior
              simply matches the distribution of input models
           error_est: bool
              if True, an estimate of the error in the numerical
              convolution of the observaitonal uncertainties with the
              model is returned; ignored if logsfr_err is None
           gkorder : string
              order Gauss-Kronrod quadrature; allowed values are '15',
              '21', '31', '41', '51', '61'; only used if logsfr_err is
              not None

        Returns
           A namedtuple consisting of:
           logsfr : array, shape (nmesh)
              an array of log_10 SFR values at which the PDF is evaluated
           sfrpdf : array, shape (nmesh) or shape (nmesh, ndata)
              the PDF of log_10 SFR evaluate at the points in the
              logsfr array; if logsfr_est is a float, this will be a
              1D array, while if logsr_est is an array whose leading
              dimension has ndata elements, it will be a 2D array
              where entry [:,M] gives the PDF for the Mth input data
              value
           pdf_err : array, shape (nmesh) or shape (nmesh, ndata)
              an estimate of the numerical error in pdf; returned only
              if error_est is True

        Raises
           ValueError, if filter_name is not one of the filters
           available in dataset_filters
        """

        # Build the KDE for the joint PDF of SFR and filter-based SFR
        if filter_name is None:
            if self.kde_filters is None:
                raise ValueError("must set filter_name if none "+
                                 "are stored")
            filter_name = self.kde_filters
        else:
            if hasattr(filter_name, '__iter__'):
                fname = ['SFR'] + filter_name
            else:
                fname = ['SFR', filter_name]
            self.set_kde(fname)

        # Set up the grid of log SFR values on which the answer will
        # be returned
        if logsfr_lim is None:
            logsfr_lim = self.kde_limits[0,:]
        logsfr = np.linspace(logsfr_lim[0], logsfr_lim[1], nmesh)
        dlogsfr = logsfr[1] - logsfr[0]

        # Evalute the ratio p(log SFR)/p_M(log SFR)
        if prior != None:
            sfr_kde, sfr_kdelim = self.get_kde('SFR')
            priorfac = prior(logsfr) \
                       / (1.0e-100 + 
                          np.exp(sfr_kde.score_samples(logsfr[:,np.newaxis])))
        else:
            priorfac = 1.0

        # Get KDE to estimate PDF of photometrically-inferred SFRs
        phot_kde, phot_kdelim = self.get_kde(filter_name)

        # Were we given observational errors?
        if logsfr_err is None:

            # Case without errors: in this case we just treat the
            # input photometric SFR estimate as a delta function, and
            # evaluate all quantities at that value

            # Make new array containing output SFRs to be evaluated
            # meshed within SFR estimates
            if not hasattr(logsfr_est, '__iter__'):
                # Input log SFR estimnate is a single number
                ndata = 1
                kdedata = np.zeros((nmesh,2))
                kdedata[:,0] = logsfr
                kdedata[:,1] = logsfr_est
            elif logsfr_est.ndim == 1:
                # Input is a 1D array giving log SFR estimated from a
                # single filter
                ndata = len(logsfr_est)
                logsfr_grd, logsfr_est_grd = np.meshgrid(logsfr, logsfr_est)
                kdedata = np.append(logsfr_grd.ravel()[:,np.newaxis], 
                                    logsfr_est_grd.ravel()[:,np.newaxis], 
                                    axis=1)
            else:
                # Input is a 2D array giving log SFR estimate from
                # multiple filters
                ndata = logsfr_est.shape[0]
                nphot = logsfr_est.shape[1]
                logsfr_grd, logsfr_est_grd \
                    = np.meshgrid(logsfr, logsfr_est[:,0])
                kdedata = np.append(logsfr_grd.ravel()[:,np.newaxis], 
                                    logsfr_est_grd.ravel()[:,np.newaxis], 
                                    axis=1)
                for i in range(nphot-1):
                    kdedata = np.append(kdedata,
                                        np.repeat(logsfr_est[:,i+1],
                                                  nmesh)[:,np.newaxis],
                                        axis=1)

            # Compute the kernel probability at each data point
            kdeprob = np.exp(self.kde.score_samples(kdedata))

            # If given multiple input data, reshape to 2D array
            if ndata > 1:
                sfrpdf = kdeprob.reshape(ndata, nmesh)
            else:
                sfrpdf = kdeprob

        else:

            # Case with errors: in this case for every input data
            # point we need to convolve the PDF with the luminosity /
            # photometric SFR distribution implied by the errors

            xgkvec = gk[gkorder]['xgkvec']
            wgkvec = gk[gkorder]['wgkvec']
            wgvec = gk[gkorder]['wgvec']
            gkord = gk[gkorder]['gkorder']

            # Step 1: figure out how many data points and photometric
            # values per data point we have
            if not hasattr(logsfr_est, '__iter__'):
                # Input log SFR estimate is a single number
                ndata = 1
                nphot = 1
            elif logsfr_est.ndim == 1:
                # Input is a 1D array giving log SFR estimated from a
                # single filter
                ndata = len(logsfr_est)
                nphot = 1
            else:
                ndata = logsfr_est.shape[0]
                nphot = logsfr_est.shape[1]

            # Step 2: get p(log SFR_est); this will be our prior on
            # the photometry-estimated SFR
            logsfr_est_kde, lim = self.get_kde(filter_name)

            # Step 3: for each input data point, generate a grid of
            # luminosity values from -6 to +6 sigma around it in each
            # direction. The placement of grid points is set up to do
            # Gauss-Kronrod integration. At each point on this grid,
            # evaluate p(L | data) / p(L).
            if nphot == 1:
                gkgrid = xgkvec
                logsfr_est_grid = logsfr_est + \
                                  np.multiply.outer(gkgrid,
                                                    5*logsfr_err)
                probdata = np.exp(-(6.0*gkgrid)**2/2.0)
                priordata = np.exp(logsfr_est_kde.score_samples(
                    logsfr_est_grid.reshape(ndata*gkord,1)). \
                                reshape(gkord, ndata))+1e-100
                pdata = probdata / np.transpose(priordata)
            else:
                veclist = [xgkvec]*nphot
                gkgrid = np.array(np.meshgrid(*veclist))

            # Step 4: generate a grid of points at each combination of
            # SFR and luminosity, and evaluate p(log SFR, log L) at
            # each of those points
            if nphot == 1:
                logsfr_grd, full_grid \
                    = np.meshgrid(logsfr, np.ravel(logsfr_est_grid))
                full_grid = full_grid.reshape(gkord, ndata, nmesh)
                logsfr_grd = logsfr_grd.reshape(gkord, ndata, nmesh)
                kdedata = np.append(
                    logsfr_grd.reshape(logsfr_grd.shape+(1,)),
                    full_grid.reshape(full_grid.shape+(1,)),
                    axis=logsfr_grd.ndim)
                kdeprob = np.exp(self.kde.score_samples(kdedata))
            else:
                pass

            # Step 5: multiply by the prior on log L to get the full
            # integrand
            integrand = np.transpose(kdeprob)*pdata

            # Step 6: last step: perform a weighted sum of the
            # integrand at each input photometric SFR and each output
            # SFR to get the quadrature estimate of the integral
            sfrpdf = np.transpose(np.sum(integrand*wgkvec, axis=2))
            if error_est:
                sfrpdf_gauss = np.transpose(
                    np.sum(integrand[:,:,1::2]*wgvec, axis=2))
            if not hasattr(logsfr_est, '__iter__'):
                sfrpdf=np.ravel(sfrpdf)
                if error_est:
                    sfrpdf_gauss=np.ravel(sfrpdf_gauss)

        # Adjust prior on SFR
        sfrpdf = sfrpdf * priorfac

        # Normalize so sum over bins is unity
        norm = np.sum(sfrpdf, axis=sfrpdf.ndim-1)*dlogsfr+1.0e-100
        if ndata > 1:
            sfrpdf = np.transpose(sfrpdf)/norm
        else:
            sfrpdf = sfrpdf/norm

        # Repeat for Gaussian estimate if we're returning an error
        if error_est and logsfr_err is not None:
            sfrpdf_gauss = sfrpdf_gauss * priorfac
            norm = np.sum(sfrpdf_gauss, 
                          axis=sfrpdf_gauss.ndim-1)*dlogsfr+1.0e-100
            if ndata > 1:
                sfrpdf_gauss = np.transpose(sfrpdf_gauss)/norm
            else:
                sfrpdf_gauss = sfrpdf_gauss/norm


        # Return
        if error_est and logsfr_err is not None:
            out_type = namedtuple('SFR_estimate',
                                  ['logsfr', 'sfrpdf', 'sfrpdf_err'])
            out = out_type(logsfr, sfrpdf, np.abs(sfrpdf-sfrpdf_gauss))
        else:
            out_type = namedtuple('SFR_estimate',
                                  ['logsfr', 'sfrpdf'])
            out = out_type(logsfr, sfrpdf)
        return out
