"""
This defines a class that can be used to estimate the PDF of physical
quantities from a set of input photometry in various bands, together
with a training data set.
"""

import numpy as np
import os
import os.path as osp
import ctypes
from ctypes import POINTER
from ctypes import c_void_p
from ctypes import c_int
from ctypes import c_uint
from ctypes import c_double
from ctypes import c_bool
import numpy.ctypeslib as npct
import random
from copy import deepcopy
try:
    import emcee
    mc_avail = True
except:
    mc_avail = False
    pass

##################################################################
# Define some types for use later                                #
##################################################################
array_1d_double = npct.ndpointer(dtype=np.double, ndim=1,
                                 flags="CONTIGUOUS")
array_1d_uint = npct.ndpointer(dtype=np.uintc, ndim=1,
                               flags="CONTIGUOUS")

##################################################################
# Define the cluster_slug class                                  #
##################################################################

class bp(object):
    """
    A class that can be used to estimate the PDF of the physical
    properties of stellar population from a training set plus a set of
    measured photometric values.
    """

    ##################################################################
    # Initializer method
    ##################################################################
    def __init__(self, dataset, nphys, filters=None, bandwidth='auto',
                 ktype='gaussian', priors=None, sample_density=None,
                 reltol=1.0e-3, abstol=1.0e-10, leafsize=16):
        """
        Initialize a cluster_slug object.

        Parameters
           dataset : array, shape (N, M)
              training data set; this is a set of N sample stellar
              populations, having M properties each; the first npys
              represent physical properties (e.g., log mass, log age),
              while the next M - nphys represent photometric
              properties
           npys : int
              number of physical properties in dataset
           filters : listlike of strings, nphys elements
              names of photometric filters; if left as None, the
              filters can be referred to by index, but not by name
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
           IOError, if the bayesphot c library cannot be found

        """

        # Load the c library
        self.__clib = npct.load_library("bayesphot", 
                                        osp.realpath(__file__))

        # Check for diagnostic mode
        self.__clib.diagnostic_mode.restype = c_bool
        self.__clib.diagnostic_mode.argtypes = None
        self.__diag_mode = bool(self.__clib.diagnostic_mode())

        # Define interfaces to all the c library functions
        self.__clib.build_kd.restype = c_void_p
        self.__clib.build_kd.argtypes \
            = [ array_1d_double,   # x
                c_uint,            # ndim
                c_uint,            # npt
                ctypes.
                POINTER(c_double), # wgt
                c_uint,            # leafsize
                array_1d_double,   # bandwidth
                c_int ]            # ktype

        self.__clib.free_kd.restype = None
        self.__clib.free_kd.argtypes = [ c_void_p ]

        self.__clib.kd_change_wgt.restype = None
        self.__clib.kd_change_wgt.argtypes \
            = [ array_1d_double,   # wgt
                c_void_p ]         # kd

        self.__clib.kd_change_bandwidth.restype = None
        self.__clib.kd_change_bandwidth.argtypes \
            = [ array_1d_double,   # bandwidth
                c_void_p ]         # kd

        self.__clib.kd_neighbors.restype = None
        self.__clib.kd_neighbors.argtypes \
            = [ c_void_p,          # kd
                array_1d_double,   # xpt
                POINTER(c_uint),   # dims
                c_uint,            # ndim
                c_uint,            # nneighbor
                c_bool,            # bandwidth_units
                array_1d_double,   # pos
                POINTER(c_double), # dptr
                array_1d_double ]  # d2
        self.__clib.kd_neighbors_vec.restype = None
        self.__clib.kd_neighbors_vec.argtypes \
            = [ c_void_p,          # kd
                array_1d_double,   # xpt
                POINTER(c_uint),   # dims
                c_uint,            # ndim
                c_uint,            # npt
                c_uint,            # nneighbor
                c_bool,            # bandwidth_units
                array_1d_double,   # pos
                POINTER(c_double), # dptr
                array_1d_double ]  # d2

        self.__clib.kd_neighbors_point.restype = None
        self.__clib.kd_neighbors_point.argtypes \
            = [ c_void_p,          # kd
                c_uint,            # idxpt
                c_uint,            # nneighbor
                c_bool,            # bandwidth_units
                array_1d_uint,     # idx
                array_1d_double ]  # d2
        self.__clib.kd_neighbors_point_vec.restype = None
        self.__clib.kd_neighbors_point_vec.argtypes \
            = [ c_void_p,          # kd
                array_1d_uint,     # idxpt
                c_uint,            # npt
                c_uint,            # nneighbor
                c_bool,            # bandwidth_units
                array_1d_uint,     # idx
                array_1d_double ]  # d2

        self.__clib.kd_neighbors_all.restype = None
        self.__clib.kd_neighbors_all.argtypes \
            = [ c_void_p,          # kd
                c_uint,            # nneighbor
                c_bool,            # bandwidth_units
                array_1d_uint,     # idx
                array_1d_double ]  # d2

        self.__clib.kd_pdf.restype = c_double
        if self.__diag_mode:
            self.__clib.kd_pdf.argtypes \
                = [ c_void_p,          # kd
                    array_1d_double,   # x
                    c_double,          # reltol
                    c_double,          # abstol
                    ctypes.POINTER(
                        c_uint),       # nodecheck
                    ctypes.POINTER(
                        c_uint),       # leafcheck
                    ctypes.POINTER(
                        c_uint) ]      # termcheck
        else:
            self.__clib.kd_pdf.argtypes \
                = [ c_void_p,          # kd
                    array_1d_double,   # x
                    c_double,          # reltol
                    c_double ]         # abstol

        self.__clib.kd_pdf_vec.restype = None
        if self.__diag_mode:
            self.__clib.kd_pdf_vec.argtypes \
                = [ c_void_p,          # kd
                    array_1d_double,   # x
                    c_uint,            # npt
                    c_double,          # reltol
                    c_double,          # abstol
                    array_1d_double,   # pdf
                    array_1d_uint,     # nodecheck
                    array_1d_uint,     # leafcheck
                    array_1d_uint ]    # termcheck
        else:
            self.__clib.kd_pdf_vec.argtypes \
                = [ c_void_p,          # kd
                    array_1d_double,   # x
                    c_uint,            # npt
                    c_double,          # reltol
                    c_double,          # abstol
                    array_1d_double ]  # pdf
        self.__clib.kd_pdf_int.restype = c_double
        if self.__diag_mode:
            self.__clib.kd_pdf_int.argtypes \
                = [ c_void_p,          # kd
                    array_1d_double,   # x
                    array_1d_uint,     # dims
                    c_uint,            # ndim
                    c_double,          # reltol
                    c_double,          # abstol
                    array_1d_uint,     # nodecheck
                    array_1d_uint,     # leafcheck
                    array_1d_uint ]    # termcheck
        else:
            self.__clib.kd_pdf_int.argtypes \
                = [ c_void_p,          # kd
                    array_1d_double,   # x
                    array_1d_uint,     # dims
                    c_uint,            # ndim
                    c_double,          # reltol
                    c_double ]         # abstol
        self.__clib.kd_pdf_int_vec.restype = None
        if self.__diag_mode:
            self.__clib.kd_pdf_int_vec.argtypes \
                = [ c_void_p,          # kd
                    array_1d_double,   # x
                    array_1d_uint,     # dims
                    c_uint,            # ndim
                    c_uint,            # npt
                    c_double,          # reltol
                    c_double,          # abstol
                    array_1d_double,   # pdf
                    array_1d_uint,     # nodecheck
                    array_1d_uint,     # leafcheck
                    array_1d_uint ]    # termcheck
        else:
            self.__clib.kd_pdf_int_vec.argtypes \
                = [ c_void_p,          # kd
                    array_1d_double,   # x
                    array_1d_uint,     # dims
                    c_uint,            # ndim
                    c_uint,            # npt
                    c_double,          # reltol
                    c_double,          # abstol
                    array_1d_double ]  # pdf

        # Record some of the input parameters
        if (ktype == 'gaussian'):
            self.__ktype = 2
        elif (ktype == 'epanechnikov'):
            self.__ktype = 0
        elif (ktype == 'tophat'):
            self.__ktype = 1
        self.leafsize = leafsize
        self.abstol = abstol
        self.reltol = reltol
        self.__sden = sample_density
        self.__sample_density = None

        # Store data set
        self.__dataset = np.ascontiguousarray(dataset)

        # Store list of available filters
        self.__filters = deepcopy(filters)

        # Initialize internal data
        self.__ndata = self.__dataset.shape[0]
        self.__nphys = nphys
        self.__nphot = self.__dataset.shape[1] - self.__nphys
        self.__auto_bw = None
        self.__auto_bw_set = False
        self.__priors = None
        self.__kd_phys = None

        # Build the initial kernel density estimation object, using a
        # dummy bandwidth
        self.__bandwidth = np.ones(self.__nphys + self.__nphot)
        self.__kd = self.__clib.build_kd(np.ravel(self.__dataset),
                                         self.__dataset.shape[1],
                                         self.__ndata, None, 
                                         leafsize, self.__bandwidth,
                                         self.__ktype)

        # Initialize the bandwidth
        self.bandwidth = bandwidth

        # Set priors
        self.priors = priors


    ##################################################################
    # De-allocation method
    ##################################################################
    def __del__(self):
        if self.__kd is not None:
            self.__clib.free_kd(self.__kd)
        if self.__kd_phys is not None:
            self.__clib.free_kd(self.__kd_phys)


    ##################################################################
    # Return a copy of the list of available filters
    ##################################################################
    def filters(self):
        return deepcopy(self.__filters)


    ##################################################################
    # Define the priors property
    ##################################################################

    @property
    def priors(self):
        return self.__priors

    @priors.setter
    def priors(self, pr):
        """
        This function sets the prior probabilities to use

        Parameters:
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

        Returns
           Nothing
        """

        # If the prior is unchanged, do nothing
        if (type(pr) == np.ndarray) and \
           (type(self.__priors) == np.ndarray):
            if np.array_equal(pr, self.__priors):
                return
        elif pr == self.__priors:
            return

        # If priors is None, just remove all weighting
        if pr is None:
            self.__clib.kd_change_wgt(None, self.__kd)
            self.__priors = None
            return

        else:
            # If we're here, we have a non-trival prior

            # Evaluate the raw sample density at each point if we have
            # not previously done so
            if self.__sample_density is None:

                # Choose computation method
                if self.__sden is None:

                    # None means uniform sampling
                    self.__sample_density = np.ones(self.__ndata)

                elif hasttr(self.__sden, '__call__'):

                    # Callable, so pass the physical data to the
                    # callable and store the result
                    self.__sample_density \
                        = self.__sden(self.__dataset[:,:self.__nphys])

                elif type(self.__sden) is np.ndarray:

                    # Array, so treat treat this as the data
                    self.__sample_densty = self.__sden

                elif self.__sden == 'auto':

                    # We've been asked to calculate the sample density
                    # ourselves, so do so

                    # Create unweighted kernel density object for just
                    # the physical parameters if we have not done so
                    # already
                    if self.__kd_phys is None:
                        self.__dataset_phys \
                            = np.copy(self.__dataset[:,:self.__nphys])
                        self.__kd_phys \
                            = self.__clib.build_kd(
                                np.ravel(self.__dataset_phys), 
                                self.__nphys, self.__ndata,
                                None, leafsize, self.__bandwidth,
                                self.__ktype)

                        # Use the unweighted kernel density object to
                        # evaluate the raw sample density near each data
                        # point
                        self.__sample_density = np.zeros(self.__ndata)
                        if not self.__diag_mode:
                            self.__clib.kd_pdf_vec(
                                self.__kd_phys, np.ravel(self.__dataset_phys),
                                self.__ndata, self.reltol, self.abstol,
                                self.__sample_density)
                        else:
                            nodecheck = np.zeros(self.__ndata, dtype=c_uint)
                            leafcheck = np.zeros(self.__ndata, dtype=c_uint)
                            termcheck = np.zeros(self.__ndata, dtype=c_uint)
                            self.__clib.kd_pdf_vec(
                                self.__kd_phys, np.ravel(self.__dataset_tmp),
                                self.__ndata, self.reltol, self.abstol,
                                self.__sample_density, nodecheck, leafcheck, 
                                termcheck)

            # We now have the sample density. Record the new prior,
            # and, if our prior is a callable, call it; otherwise
            # just record the input data
            self.__priors = pr
            if hasattr(self.__priors, '__call__'):
                prior_data = self.__priors(self.__dataset[:,:self.__nphys])
            else:
                prior_data = self.__priors

            # Compute the weights from the ratio of the prior to
            # the sample density, then adjust the weights in the kd
            self.__wgt = prior_data / self.__sample_density
            self.__clib.kd_change_wgt(self.__wgt, self.__kd)


    ##################################################################
    # Define the bandwidth property
    ##################################################################
 
    @property
    def bandwidth(self):
        return deepcopy(self.__bandwidth)

    @bandwidth.setter
    def bandwidth(self, bw):

        if np.array_equal(self.__bandwidth, bw):
            # If new bandwidth equals old bandwidth, do nothing
            return

        elif bw != 'auto':
            # If we've been given a specified bandwidth, set to that
            self.__bandwidth = np.copy(bw)
            self.__auto_bw_set = False

        else:
            # Automatic bandwidth setting

            # Are we already set on auto? If so, just return
            if self.__auto_bw_set:
                return

            # Do we have a stored value for the automatic bandwidth?
            # If not, we need to compute it.
            if self.__auto_bw is None:

                # Find 10th nearest neighbors
                nneighbor=10
                if self.__ndata > 5000:
                    # For data sets with > 5000 samples, just use a
                    # sub-sample of 5,000, which is more than enough
                    # to get a reasonable estimate of the distribution
                    idxpt = np.array(
                        random.sample(np.arange(self.__ndata), 
                                      5000), dtype=np.uintc)
                    neighbors = np.zeros(nneighbor*5000, 
                                              dtype=np.uintc)
                    d2 = np.zeros(nneighbor*5000)
                    self.__clib.kd_neighbors_point_vec(
                        self.__kd, idxpt, 5000, nneighbor, False,
                        neighbors, d2)

                else:
                    # For smaller data sets, use it all
                    neighbors = np.zeros(nneighbor*self.__ndata, 
                                         dtype=np.uintc)
                    d2 = np.zeros(nneighbor*self.__ndata)
                    idxpt = np.arange(self.__ndata, dtype=np.uintc)
                    self.__clib.kd_neighbors_all(self.__kd, nneighbor, 
                                                 False, neighbors, d2)

                # Take the bandwidth in each dimension to be the 90th
                # percentile of the 10th nearest neighbor distance
                offset = np.abs(self.__dataset[idxpt,:] -
                                self.__dataset[
                                    neighbors[nneighbor-1::nneighbor],:])
                self.__auto_bw = np.zeros(self.__nphys+self.__nphot)
                for i in range(self.__nphys+self.__nphot):
                    self.__auto_bw[i] = np.percentile(offset[:,i], 95)

            # Set to the auto bandwidth
            self.__bandwidth = np.copy(self.__auto_bw)
            self.__auto_bw_set = True

        # Set the new bandwidth
        self.__clib.kd_change_bandwidth(self.__bandwidth, self.__kd)

        # If we have stored priors, and we're using automatic sample
        # density setting, we need to recompute them for the
        # new bandwidth; zero out the stored sample density, and
        # adjust the kernel density estimator for the physical
        # parameters to the new bandwidth before doing so
        if self.__sden == 'auto':
            self.__sample_density = None
            if self.__kd_phys is not None:
                self.__clib.kd_change_bandwidth(
                    self.__bandwidth[:self.__nphys], self.__kd_phys)
            pr = self.priors
            self.priors = None
            self.priors = pr


    ##################################################################
    # Method to compute the log likelihood function for a particular
    # set of physical properties given a particular set of photometric
    # properties
    ##################################################################
    def logL(self, physprop, photprop, photerr=None):
        """
        This function returns the natural log of the likelihood
        function evaluated at a particular log mass, log age,
        extinction, and set of log luminosities

        Parameters:
           physprop : arraylike, shape (nphys) or (..., nphys)
              array giving values of the physical properties; for a
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

        Returns:
           logL : float or arraylike
              natural log of the likelihood function

           (only returned if the c library was compiled in DIAGNOSTIC mode)
           nodecheck : int or array, shape(N)
              Number of nodes examined during the evaluation
           leafcheck : int or array, shape(N)
              Number of leaves examined during the evaluation
           termcheck : int or array, shape(N)
              Number of nodes examined during the evaluation for which
              no children were examined
        """

        # Safety check
        if np.array(physprop).shape[-1] != self.__nphys:
            raise ValueError("need " + str(self.__nphys) + 
                             " physical properties!")
        if np.array(photprop).shape[-1] != self.__nphot:
            raise ValueError("need " + str(self.__nphot) +
                             " photometric properties!")
        if photerr is not None:
            if np.array(photerr).shape[-1] != self.__nphot:
                raise ValueError("need " + str(self.__nphot) +
                                 " photometric errors!")

        # Figure out number of distinct input sets of physical and
        # photometric properties
        nphys_in = np.array(physprop).size/self.__nphys
        nphot_in = np.array(photprop).size/self.__nphot

        # Figure out how many sets of photometric errors we have. We
        # unfortunately have to iterate over these, because each
        # distinct set requires changing the bandwidth of the kernel
        # density estimation.
        nphot_err = np.array(photerr).size/self.__nphot

        # Allocate an array to hold the results
        if photerr is not None:
            pdf = np.zeros(
                np.broadcast(np.array(physprop)[..., 0],
                             np.array(photprop)[..., 0],
                             np.array(photerr)[..., 0]).shape)
        else:
            pdf = np.zeros(
                np.broadcast(np.array(physprop)[..., 0],
                             np.array(photprop)[..., 0]).shape)
        if self.__diag_mode:
            nodecheck = np.zeros(pdf.shape, dtype=c_uint)
            leafcheck = np.zeros(pdf.shape, dtype=c_uint)
            termcheck = np.zeros(pdf.shape, dtype=c_uint)

        # Make an array suitable for passing data to c routines
        cdata = np.zeros((pdf.size, self.__nphys+self.__nphot))
        cdata[..., :self.__nphys] \
            = np.vstack((physprop,) * (pdf.size/nphys_in))
        cdata[..., self.__nphys:] \
            = np.vstack((photprop,) * (pdf.size/nphot_in))

        # Separate cases with single / no photometric errors from
        # cases with multiple sets of photometric errors
        if nphot_err <= 1:

            # Case with at most one set of photometric errors

            # Set the bandwidth based on the photometric errors if we
            # were given some
            if photerr is not None:
                err = np.zeros(self.__bandwidth.size)
                err[self.__nphys:] = photerr
                bandwidth = np.sqrt(self.__bandwidth**2+err**2)
                self.__clib.kd_change_bandwidth(bandwidth, self.__kd)

            # Call the PDF computation routine
            if not self.__diag_mode:
                self.__clib.kd_pdf_vec(
                    self.__kd, np.ravel(cdata), pdf.size, 
                    self.reltol, self.abstol, np.ravel(pdf))
            else:
                self.__clib.kd_pdf_vec(
                    self.__kd, np.ravel(cdata), pdf.size, 
                    self.reltol, self.abstol, np.ravel(pdf),
                    np.ravel(nodecheck), np.ravel(leafcheck),
                    np.ravel(termcheck))

            # Set the bandwidth back to its default if necessary
            if photerr is not None:
                self.__clib.kd_change_bandwidth(self.__bandwidth, 
                                                self.__kd)

            # Return
            if not self.__diag_mode:
                return pdf
            else:
                return pdf, nodecheck, leafcheck, termcheck

        else:

            # Case with multiple sets of photometric errors

            # Loop over photometric errors
            for i in np.ndindex(*photerr.shape[:-1]):

                # Set bandwidth based on photometric error for this
                # iteration
                err = np.zeros(self.__bandwidth.size)
                err[self.__nphys:] = photerr[i]
                bandwidth = np.sqrt(self.bandwidth**2+err**2)
                self.__clib.kd_change_bandwidth(bandwidth, self.__kd)

                # Grab the corresponding portions of the arrays going
                # to and from the c code
                cdata_sub = cdata[i]
                pdf_sub = np.zeros(np.array(pdf[i]).shape)
                if self.__diag_mode:
                    nodecheck_sub = np.zeros(np.array(nodecheck[i]).shape)
                    leafcheck_sub = np.zeros(np.array(leafcheck[i]).shape)
                    termcheck_sub = np.zeros(np.array(termcheck[i]).shape)

                # Call kernel density estimate with this bandwidth
                if not self.__diag_mode:
                    self.__clib.kd_pdf_vec(
                        self.__kd, np.ravel(cdata_sub), pdf_sub.size, 
                        self.reltol, self.abstol, np.ravel(pdf_sub))
                else:
                    self.__clib.kd_pdf_vec(
                        self.__kd, np.ravel(cdata_sub), pdf_sub.size, 
                        self.reltol, self.abstol, np.ravel(pdf_sub),
                        np.ravel(nodecheck_sub), np.ravel(leafcheck_sub),
                        np.ravel(termcheck_sub))
                pdf[i] = pdf_sub
                if self.__diag_mode:
                    nodecheck[i] = nodecheck_sub
                    leafcheck[i] = leafcheck_sub
                    termcheck[i] = termcheck_sub

            # Restore the bandwidth
            self.__clib.kd_change_bandwidth(self.__bandwidth, 
                                            self.__kd)

            # Return
            if not self.__diag_mode:
                return pdf
            else:
                return pdf, nodecheck, leafcheck, termcheck


    ##################################################################
    # Function to return the marginal distribution of one of the
    # physical properties for a specified set of photometric
    # properties
    ##################################################################
    def mpdf(self, idx, photprop, photerr=None, ngrid=128,
             qmin=None, qmax=None, grid=None, norm=True):
        """
        Returns the marginal probability for one or mode physical
        quantities for one or more input sets of photometric
        properties. Output quantities are computed on a grid of
        values, in the same style as meshgrid

        Parameters:
           idx : int or listlike containing ints
              index of the physical quantity whose PDF is to be
              computed; if this is an iterable, the joint distribution of
              the indicated quantities is returned
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

        # Safety check
        if np.array(photprop).shape[-1] != self.__nphot:
            raise ValueError("need " + str(self.__nphot) +
                             " photometric properties!")
        if photerr is not None:
            if np.array(photerr).shape[-1] != self.__nphot:
                raise ValueError("need " + str(self.__nphot) +
                                 " photometric errors!")

        # Set up the grid of outputs
        if grid is not None:
            grid_out = grid
        else:
            if qmin is None:
                qmin = np.amin(self.__dataset[:,idx], axis=0)
            if qmax is None:
                qmax = np.amax(self.__dataset[:,idx], axis=0)
            griddims = []
            if hasattr(idx, '__len__'):
                # Case for multiple indices
                griddims = []
                if hasattr(ngrid, '__len__'):
                    ngrid_tmp = ngrid
                else:
                    ngrid_tmp = [ngrid]*len(idx)
                for i in range(len(idx)):
                    griddims.append(qmin[i] + np.arange(ngrid_tmp[i]) * 
                                    float(qmax[i]-qmin[i])/(ngrid_tmp[i]-1))
                grid_out = np.array(np.meshgrid(*griddims,
                                                indexing='ij'))
                out_shape = grid_out[0, ...].shape
            else:
                # Case for a single index
                grid_out = qmin + \
                           np.arange(ngrid) * \
                           float(qmax-qmin)/(ngrid-1)
                out_shape = grid_out.shape

        # Figure out how many distinct photometric values we've been
        # given, and how many sets of photometric errors
        nphot = np.array(photprop).size/self.__nphot
        nphot_err = np.array(photerr).size/self.__nphot

        # Set up a grid to hold the outputs
        if photerr is not None:
            pdf = np.zeros(
                np.broadcast(np.array(photprop)[..., 0],
                             np.array(photerr)[..., 0]).shape +
                out_shape)
        else:
            pdf = np.zeros(np.array(photprop)[..., 0].shape +
                           out_shape)
        if self.__diag_mode:
            nodecheck = np.zeros(pdf.shape, dtype=c_uint)
            leafcheck = np.zeros(pdf.shape, dtype=c_uint)
            termcheck = np.zeros(pdf.shape, dtype=c_uint)

        # Prepare data for c library
        if hasattr(idx, '__len__'):
            nidx = len(idx)
            ng = grid_out[nidx-1:].size / max(nphot_err, 1)
        else:
            nidx = 1
            ng = grid_out.size / max(nphot_err, 1)
        dims = np.zeros(nidx+self.__nphot, dtype=np.uintc)
        dims[:nidx] = idx
        dims[nidx:] = self.__nphys+np.arange(self.__nphot, dtype='int')
        ndim = np.uintc(nidx + self.__nphot)
        npt = np.uintc(ng)
        cdata = np.zeros(pdf.shape+(ndim,))
        if nidx == 1:
            cdata[..., 0] = np.vstack((grid_out,) * 
                                      photprop[...,0].size). \
                reshape(cdata[..., 0].shape)
            cdata[..., nidx:] = np.hstack((photprop,) * grid_out.size). \
                                reshape(cdata[..., nidx:].shape)
        else:
            for i in range(nidx):
                cdata[..., i] = np.vstack((grid_out[i,...],) * 
                                          photprop[...,0].size). \
                    reshape(cdata[..., i].shape)
            cdata[..., nidx:] = np.hstack((photprop,) * 
                                          grid_out[0,...].size). \
                                reshape(cdata[..., nidx:].shape)

        # Separate cases with single / no photometric errors from
        # cases with multiple sets of photometric errors
        if nphot_err <= 1:

            # Case with at most one set of photometric errors

            # Set the bandwidth based on the photometric errors if we
            # were given some
            if photerr is not None:
                err = np.zeros(self.__bandwidth.size)
                err[self.__nphys:] = photerr
                bandwidth = np.sqrt(self.__bandwidth**2+err**2)
                self.__clib.kd_change_bandwidth(bandwidth, self.__kd)

            # Call the PDF computation routine
            if not self.__diag_mode:
                self.__clib.kd_pdf_int_vec(self.__kd, np.ravel(cdata),
                                           dims, ndim, pdf.size,
                                           self.reltol, self.abstol,
                                           np.ravel(pdf))
            else:
                self.__clib.kd_pdf_int_vec(self.__kd, np.ravel(cdata),
                                           dims, ndim, pdf.size,
                                           self.reltol, self.abstol,
                                           np.ravel(pdf), 
                                           np.ravel(nodecheck), 
                                           np.ravel(leafcheck),
                                           np.ravel(termcheck))

        else:

            # Case with multiple sets of photometric errors

            # Loop over photometric errors
            for i in np.ndindex(*photerr.shape[:-1]):

                # Set bandwidth based on photometric error for this
                # iteration
                err = np.zeros(self.__bandwidth.size)
                err[self.__nphys:] = photerr[i]
                bandwidth = np.sqrt(self.bandwidth**2+err**2)
                self.__clib.kd_change_bandwidth(bandwidth, self.__kd)

                # Grab the corresponding portions of the arrays going
                # to and from the c code
                cdata_sub = cdata[i]
                pdf_sub = np.zeros(np.array(pdf[i]).shape)
                if self.__diag_mode:
                    nodecheck_sub = np.zeros(np.array(nodecheck[i]).shape)
                    leafcheck_sub = np.zeros(np.array(leafcheck[i]).shape)
                    termcheck_sub = np.zeros(np.array(termcheck[i]).shape)

                # Call kernel density estimate with this bandwidth
                if not self.__diag_mode:
                    self.__clib.kd_pdf_int_vec(
                        self.__kd, np.ravel(cdata_sub), dims, ndim,
                        pdf_sub.size, self.reltol, self.abstol, 
                        np.ravel(pdf_sub))
                else:
                    self.__clib.kd_pdf_int_vec(
                        self.__kd, np.ravel(cdata_sub), dims, ndim,
                        pdf_sub.size, self.reltol, self.abstol, 
                        np.ravel(pdf_sub), np.ravel(nodecheck_sub), 
                        np.ravel(leafcheck_sub), np.ravel(termcheck_sub))
                pdf[i] = pdf_sub
                if self.__diag_mode:
                    nodecheck[i] = nodecheck_sub
                    leafcheck[i] = leafcheck_sub
                    termcheck[i] = termcheck_sub


        # Set the bandwidth back to its default if necessary
        if photerr is not None:
            self.__clib.kd_change_bandwidth(self.__bandwidth, 
                                                self.__kd)

        # Normalize if requested
        if norm:

            # Compute the sizes of the output cells
            if nidx == 1:
                cellsize = np.zeros(grid_out.size)
                cellsize[1:-1] = 0.5*(grid_out[2:]-grid_out[:-2])
                cellsize[0] = grid_out[1] - grid_out[0]
                cellsize[-1] = grid_out[-1] - grid_out[-2]
            else:
                # Get the cell sizes in each dimension
                csize = []
                for i in range(nidx):
                    vec = grid_out[(i,)+i*(0,)+(slice(None),) + 
                                   (grid_out.shape[0]-i-1)*(0,)]
                    csize.append(np.zeros(vec.size))
                    csize[i][1:-1] = 0.5*(vec[2:]-vec[:-2])
                    csize[i][0] = vec[1] - vec[0]
                    csize[i][-1] = vec[-1] - vec[-2]
                # Take outer product to get grid of sizes
                cellsize = np.multiply.outer(csize[0], csize[1])
                for i in range(2, nidx):
                    cellsize = np.multiply.outer(cellsize, csize[i])

            # Compute integral
            normfac = np.sum(pdf*cellsize, axis = 
                             tuple(range(photprop.ndim-1, pdf.ndim)))

            # Normalize
            pdf = np.transpose(np.transpose(pdf)/normfac)

        # Return
        if not self.__diag_mode:
            return grid_out, pdf
        else:
            return grid_out, pdf, nodecheck, leafcheck, termcheck


    ##################################################################
    # Method to return log likelihood function at a specified set of
    # physical properties for a particular set of photometric
    # variables; this is set up for use by emcee, and is not intended
    # for use by humans
    ##################################################################
    def __logL(self, *args):
        x = np.zeros(self.__nphys+self.__nphot)
        x[:self.__nphys] = args[0]
        x[self.__nphys:] = args[1:]
        if not self.__diag_mode:
            return np.log(self.__clib.kd_pdf(self.__kd, x, self.reltol,
                                             self.abstol))
        else:
            nodecheck = c_uint(0)
            leafcheck = c_uint(0)
            termcheck = c_uint(0)
            return np.log(self.__clib.kd_pdf(self.__kd, x, self.reltol,
                                             self.abstol, nodecheck,
                                             leafcheck, termcheck))


    ##################################################################
    # Function to compute an MCMC sampler for a particular set of
    # photometric values
    ##################################################################
    def mcmc(self, photprop, photerr=None, mc_walkers=100,
             mc_steps=500, mc_burn_in=50):
        """
        This function returns a sample of MCMC walkers sampling the
        physical parameters at a specified set of photometric values.

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

        Returns
           samples : array
              array of sample points returned by the MCMC
        """

        # See if we have emcee
        if not mc_avail:
            raise ImportError("unable to import emcee")

        # Safety check
        if np.array(photprop).shape[-1] != self.__nphot:
            raise ValueError("need " + str(self.__nphot) +
                             " photometric properties!")
        if photerr is not None:
            if np.array(photerr).shape[-1] != self.__nphot:
                raise ValueError("need " + str(self.__nphot) +
                             " photometric properties!")
 
        # Prepare storage for output
        samples = []

        # Make dummy photometric errors if necessary
        if photerr is None:
            photerr = np.zeros(self.__nphot)

        # Loop over photometric errors
        for i in np.ndindex(*np.array(photerr).shape[:-1]):

            # Set bandwidth based on photometric error for this
            # iteration
            err = np.zeros(self.__bandwidth.size)
            err[self.__nphys:] = photerr[i]
            bandwidth = np.sqrt(self.bandwidth**2+err**2)
            self.__clib.kd_change_bandwidth(bandwidth, self.__kd)

            # Grab the clusters that go with this photometric error
            if photprop.ndim < photerr.ndim:
                ph = photprop
            else:
                ph = photprop[i]

            # Loop over clusters
            for j in np.ndindex(*np.array(ph).shape[:-1]):

                # Grab photometric values for this cluster
                ph_tmp = ph[j]

                # Search the data set for the sample cluster closest to
                # the observed luminosities; this will be the starting
                # point for the MCMC
                dims = np.arange(self.__nphys,
                                 self.__nphys+self.__nphot, 
                                 dtype=np.uint32)
                nearpt = np.zeros(self.__nphys+self.__nphot)
                wgt = np.zeros(1)
                dist2 = np.zeros(1)
                self.__clib. \
                    kd_neighbors(self.__kd, ph_tmp,
                                 dims.ctypes.data_as(POINTER(c_uint)),
                                 self.__nphot, 1, True, nearpt, 
                                 wgt.ctypes.data_as(POINTER(c_double)),
                                 dist2)

                # Generate a set of starting points by scattering walkers
                # around the starting position
                pos = [nearpt[:self.__nphys] + 
                       self.bandwidth[:self.__nphys] * 
                       np.random.randn(self.__nphys) 
                       for i in range(mc_walkers)]

                # Run the MCMC
                sampler=emcee.EnsembleSampler(mc_walkers, self.__nphys, 
                                              self.__logL, 
                                              args=ph_tmp)
                sampler.run_mcmc(pos, mc_steps)

                # Store the result
                samples.append(sampler.chain[:,mc_burn_in:,:].
                               reshape((-1,self.__nphys)))


        # Set bandwidth back to default if necessary
        if photerr is not None:
            self.__clib.kd_change_bandwidth(self.__bandwidth, 
                                                self.__kd)

        # Reshape the samples
        samples = np.squeeze(
            np.array(samples).reshape(
                np.broadcast(np.array(photprop)[..., 0],
                             np.array(photerr)[..., 0]).shape +
                samples[0].shape))

        # Return
        return samples
