"""
This defines a class that can be used to estimate the PDF of physical
quantities from a set of input photometry in various bands, together
with a training data set.
"""

import numpy as np
import scipy.interpolate as interp
import os
import os.path as osp
import ctypes
from ctypes import POINTER
from ctypes import c_void_p
from ctypes import c_int
from ctypes import c_uint
from ctypes import c_ulong
from ctypes import c_double
from ctypes import c_bool
import numpy.ctypeslib as npct
import random
from copy import deepcopy
from warnings import warn
try:
    import emcee
    mc_avail = True
except:
    mc_avail = False
    pass

##################################################################
# Define some types for use later                                #
##################################################################
array_1d_double = npct.ndpointer(dtype=c_double, ndim=1,
                                 flags="CONTIGUOUS")
array_1d_ulong = npct.ndpointer(dtype=c_ulong, ndim=1,
                               flags="CONTIGUOUS")
array_1d_int = npct.ndpointer(dtype=c_int, ndim=1,
                              flags="CONTIGUOUS")

##################################################################
# Define the cluster_slug class                                  #
##################################################################

class bp(object):
    """
    A class that can be used to estimate the PDF of the physical
    properties of stellar population from a training set plus a set of
    measured photometric values.

    Properties
       priors : array, shape (N) | callable | None
          prior probability on each data point; interpretation
          depends on the type passed; array, shape (N): values are
          interpreted as the prior probability of each data point;
          callable: the callable must take as an argument an array
          of shape (N, nphys), and return an array of shape (N)
          giving the prior probability at each data point; None:
          no reweighting is performed, so all data points in the library
          have equal prior probability, i.e. the prior is the same as
          the sampling of points in the library
       pobs : array, shape (N) | callable | None
          the probability that a particular object would be observed,
          which is used, like prior, to weight the library;
          interpretation depends on type. None means all objects are
          equally likely to be observed, array is an array giving the
          observation probability of each object in the library, and
          callable means must be a function that takes an array
          containing the photometry, of shape (N, nhpot), as an
          argument, and returns an array of shape (N) giving the
          probability of observation for that object
       bandwidth : 'auto' | float | array, shape (M)
          bandwidth for kernel density estimation; if set to
          'auto', the bandwidth will be estimated automatically; if
          set to a scalar quantity, the same bandwidth is used for all
          dimensions
       nphys : int
          number of physical properties in the library
       nphot : int
          number of photometric properties in the library
       ndim : int
          nphys + nphot
    """

    ##################################################################
    # Initializer method
    ##################################################################
    def __init__(self, dataset, nphys, filters=None, bandwidth='auto',
                 ktype='gaussian', priors=None, pobs=None,
                 sample_density=None, reltol=1.0e-2, abstol=1.0e-10,
                 leafsize=16, nosort=None, thread_safe=True,
                 caching='none'):
        """
        Initialize a bp object.

        Parameters
           dataset : array, shape (N, M)
              training data set; this is a set of N sample stellar
              populations, having M properties each; the first nphys
              represent physical properties (e.g., log mass, log age),
              while the next M - nphys represent photometric
              properties
           nphys : int
              number of physical properties in dataset
           filters : listlike of strings
              names of photometric filters; not used, but can be
              stored for convenience
           bandwidth : 'auto' | float | array, shape (M)
              bandwidth for kernel density estimation; if set to
              'auto', the bandwidth will be estimated automatically; if
              set to a scalar quantity, the same bandwidth is used for all
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
           pobs : array, shape (N) | callable | None
              the probability that a particular object would be observed,
              which is used, like prior, to weight the library;
              interpretation depends on type. None means all objects are
              equally likely to be observed, array is an array giving the
              observation probability of each object in the library, and
              callable means must be a function that takes an array
              containing the photometry, of shape (N, nhpot), as an
              argument, and returns an array of shape (N) giving the
              probability of observation for that object
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
              to be uniformly sampled
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
           nosort : arraylike of bool, shape (N) | None
              if specified, this keyword causes the KD tree not to be
              sorted along the dimensions for which nosort is True
           thread_safe : bool
              if True, bayesphot will make extra copies of internals
              as needed to ensure thread safety when the computation
              routines (logL, mpdf, mcmc, bestmatch, make_approx_phot,
              make_approx_phys, mpdf_approx) are used with
              multiprocessing; this incurs a minor performance
              penalty, and can be disabled by setting to False if the
              code will not be run with the multiprocessing module
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

        Returns
           Nothing

        Raises
           IOError, if the bayesphot c library cannot be found

        Notes
           Because the data sets passed in may be large, this class
           does not make copies of any of its arguments, and instead
           modifies them in place. However, this means it is the
           responsibility of the user not to alter the any of the
           arguments once they are passed to this class; for example,
           dataset must not be modified after it is passed. Altering the
           arguments, except through this class's methods, may cause
           incorrect results to be generated.

           Caching with 'aggressive' or 'lazy' does create significant
           memory overhead.
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
                c_ulong,           # ndim
                c_ulong,           # npt
                ctypes.
                POINTER(c_double), # wgt
                c_ulong,           # leafsize
                array_1d_double,   # bandwidth
                c_int,             # ktype
                c_ulong,           # minsplit
                POINTER(c_ulong) ] # sortmap

        self.__clib.build_kd_sortdims.restype = c_void_p
        self.__clib.build_kd_sortdims.argtypes \
            = [ array_1d_double,   # x
                c_ulong,           # ndim
                c_ulong,           # npt
                ctypes.
                POINTER(c_double), # wgt
                c_ulong,           # leafsize
                array_1d_double,   # bandwidth
                c_int,             # ktype
                array_1d_int,      # nosort
                POINTER(c_ulong) ] # sortmap

        self.__clib.copy_kd.restype = c_void_p
        self.__clib.copy_kd.argtypes = [ c_void_p ] # kd

        self.__clib.free_kd.restype = None
        self.__clib.free_kd.argtypes = [ c_void_p ]

        self.__clib.free_kd_copy.restype = None
        self.__clib.free_kd_copy.argtypes = [ c_void_p ]

        self.__clib.kd_change_wgt.restype = None
        self.__clib.kd_change_wgt.argtypes \
            = [ POINTER(c_double), # wgt
                c_void_p ]         # kd

        self.__clib.kd_change_bandwidth.restype = None
        self.__clib.kd_change_bandwidth.argtypes \
            = [ array_1d_double,   # bandwidth
                c_void_p ]         # kd

        self.__clib.kd_neighbors.restype = None
        self.__clib.kd_neighbors.argtypes \
            = [ c_void_p,          # kd
                array_1d_double,   # xpt
                POINTER(c_ulong),  # dims
                c_ulong,           # ndim
                c_ulong,           # nneighbor
                c_bool,            # bandwidth_units
                array_1d_double,   # pos
                POINTER(c_double), # dptr
                array_1d_double ]  # d2
        self.__clib.kd_neighbors_vec.restype = None
        self.__clib.kd_neighbors_vec.argtypes \
            = [ c_void_p,          # kd
                array_1d_double,   # xpt
                POINTER(c_ulong),  # dims
                c_ulong,           # ndim
                c_ulong,           # npt
                c_ulong,           # nneighbor
                c_bool,            # bandwidth_units
                array_1d_double,   # pos
                POINTER(c_double), # dptr
                array_1d_double ]  # d2

        self.__clib.kd_neighbors_point.restype = None
        self.__clib.kd_neighbors_point.argtypes \
            = [ c_void_p,          # kd
                c_ulong,           # idxpt
                c_ulong,           # nneighbor
                c_bool,            # bandwidth_units
                array_1d_ulong,     # idx
                array_1d_double ]  # d2
        self.__clib.kd_neighbors_point_vec.restype = None
        self.__clib.kd_neighbors_point_vec.argtypes \
            = [ c_void_p,          # kd
                array_1d_ulong,     # idxpt
                c_ulong,           # npt
                c_ulong,           # nneighbor
                c_bool,            # bandwidth_units
                array_1d_ulong,     # idx
                array_1d_double ]  # d2

        self.__clib.kd_neighbors_all.restype = None
        self.__clib.kd_neighbors_all.argtypes \
            = [ c_void_p,          # kd
                c_ulong,           # nneighbor
                c_bool,            # bandwidth_units
                array_1d_ulong,     # idx
                array_1d_double ]  # d2

        self.__clib.kd_pdf.restype = c_double
        if self.__diag_mode:
            self.__clib.kd_pdf.argtypes \
                = [ c_void_p,          # kd
                    array_1d_double,   # x
                    c_double,          # reltol
                    c_double,          # abstol
                    POINTER(c_ulong),  # nodecheck
                    POINTER(c_ulong),  # leafcheck
                    POINTER(c_ulong) ] # termcheck
        else:
            self.__clib.kd_pdf.argtypes \
                = [ c_void_p,          # kd
                    array_1d_double,   # x
                    c_double,          # reltol
                    c_double ]         # abstol

        self.__clib.kd_pdf_grid.restype = None
        self.__clib.kd_pdf_grid.argtypes \
            = [ c_void_p,              # kd
                array_1d_double,       # xfixed
                array_1d_ulong,         # dimfixed
                c_ulong,               # ndimfixed
                c_ulong,               # nfixed
                array_1d_double,       # xgrid
                array_1d_ulong,         # dimgrid
                c_ulong,               # ndimgrid
                c_ulong,               # ngrid
                c_double,              # reltol
                c_double,              # abstol
                array_1d_double ]      # pdf
        self.__clib.kd_pdf_reggrid.restype = None
        self.__clib.kd_pdf_reggrid.argtypes \
            = [ c_void_p,              # kd
                array_1d_double,       # xfixed
                array_1d_ulong,         # dimfixed
                c_ulong,               # ndimfixed
                c_ulong,               # nfixed
                array_1d_double,       # xgridlo
                array_1d_double,       # xgridhi
                array_1d_ulong,         # ngrid
                array_1d_ulong,         # dimgrid
                c_ulong,               # ndimgrid
                c_double,              # reltol
                c_double,              # abstol
                array_1d_double ]      # pdf

        self.__clib.kd_pdf_vec.restype = None
        if self.__diag_mode:
            self.__clib.kd_pdf_vec.argtypes \
                = [ c_void_p,          # kd
                    array_1d_double,   # x
                    POINTER(c_double), # bandwidth
                    c_ulong,           # npt
                    c_double,          # reltol
                    c_double,          # abstol
                    array_1d_double,   # pdf
                    array_1d_ulong,     # nodecheck
                    array_1d_ulong,     # leafcheck
                    array_1d_ulong ]    # termcheck
        else:
            self.__clib.kd_pdf_vec.argtypes \
                = [ c_void_p,          # kd
                    array_1d_double,   # x
                    POINTER(c_double), # bandwidth
                    c_ulong,           # npt
                    c_double,          # reltol
                    c_double,          # abstol
                    array_1d_double ]  # pdf
        self.__clib.kd_pdf_int.restype = c_double
        if self.__diag_mode:
            self.__clib.kd_pdf_int.argtypes \
                = [ c_void_p,          # kd
                    array_1d_double,   # x
                    array_1d_ulong,    # dims
                    c_ulong,           # ndim
                    c_double,          # reltol
                    c_double,          # abstol
                    array_1d_ulong,    # nodecheck
                    array_1d_ulong,    # leafcheck
                    array_1d_ulong ]   # termcheck
        else:
            self.__clib.kd_pdf_int.argtypes \
                = [ c_void_p,          # kd
                    array_1d_double,   # x
                    array_1d_ulong,    # dims
                    c_ulong,           # ndim
                    c_double,          # reltol
                    c_double ]         # abstol
        self.__clib.kd_pdf_int_vec.restype = None
        if self.__diag_mode:
            self.__clib.kd_pdf_int_vec.argtypes \
                = [ c_void_p,          # kd
                    array_1d_double,   # x
                    POINTER(c_double), # bandwidth
                    array_1d_ulong,    # dims
                    c_ulong,           # ndim
                    c_ulong,           # npt
                    c_double,          # reltol
                    c_double,          # abstol
                    array_1d_double,   # pdf
                    array_1d_ulong,    # nodecheck
                    array_1d_ulong,    # leafcheck
                    array_1d_ulong ]   # termcheck
        else:
            self.__clib.kd_pdf_int_vec.argtypes \
                = [ c_void_p,          # kd
                    array_1d_double,   # x
                    POINTER(c_double), # bandwidth
                    array_1d_ulong,    # dims
                    c_ulong,           # ndim
                    c_ulong,           # npt
                    c_double,          # reltol
                    c_double,          # abstol
                    array_1d_double ]  # pdf
        self.__clib.kd_pdf_int_grid.restype = None
        self.__clib.kd_pdf_int_grid.argtypes \
            = [ c_void_p,              # kd
                array_1d_double,       # xfixed
                array_1d_ulong,        # dimfixed
                c_ulong,               # ndimfixed
                c_ulong,               # nfixed
                array_1d_double,       # xgrid
                array_1d_ulong,        # dimgrid
                c_ulong,               # ndimgrid
                c_ulong,               # ngrid
                c_double,              # reltol
                c_double,              # abstol
                array_1d_double ]      # pdf
        self.__clib.kd_pdf_int_reggrid.restype = None
        self.__clib.kd_pdf_int_reggrid.argtypes \
            = [ c_void_p,              # kd
                array_1d_double,       # xfixed
                array_1d_ulong,        # dimfixed
                c_ulong,               # ndimfixed
                c_ulong,               # nfixed
                array_1d_double,       # xgridlo
                array_1d_double,       # xgridhi
                array_1d_ulong,        # ngrid
                array_1d_ulong,        # dimgrid
                c_ulong,               # ndimgrid
                c_double,              # reltol
                c_double,              # abstol
                array_1d_double ]      # pdf                

        self.__clib.kd_rep.restype = c_ulong
        self.__clib.kd_rep.argtypes \
            = [ c_void_p,              # kd
                array_1d_double,       # x
                array_1d_ulong,        # dims
                c_ulong,               # ndim
                c_double,              # reltol
                POINTER(c_ulong),      # dim_return
                c_ulong,               # ndim_return
                POINTER(POINTER(
                    c_double)),        # xpt
                POINTER(POINTER(
                    c_double)) ]       # wgts
        self.__clib.free_kd_rep.restype = None
        self.__clib.free_kd_rep.argtypes \
            = [ POINTER(POINTER(
                c_double)),            # xpt
               POINTER(POINTER(
                   c_double)) ]        # wgts
        self.__clib.squeeze_rep.restype = c_ulong
        self.__clib.squeeze_rep.argtypes \
            = [ c_ulong,               # npts
                c_ulong,               # ndim
                array_1d_double,       # h
                c_double,              # tol
                POINTER(POINTER(
                    c_double)),        # xpt
                POINTER(POINTER(
                    c_double)) ]       # wgts

        self.__clib.kd_pdf_draw.restype = None
        self.__clib.kd_pdf_draw.argtypes \
            = [ c_void_p,              # kd
                POINTER(c_double),     # x
                POINTER(c_ulong),      # dims
                c_ulong,               # ndim
                c_ulong,               # nsample
                c_int,                 # draw_method
                c_void_p,              # rng
                array_1d_double ]      # out
        self.__clib.rng_init.restype = c_void_p
        self.__clib.rng_init.argtypes \
            = [ c_ulong ]              # seed
        self.__clib.rng_free.restype = None
        self.__clib.rng_free.argtypes \
            = [ c_void_p ]             # r

        # Record some of the input parameters
        if (ktype == 'gaussian'):
            self.__ktype = 2
        elif (ktype == 'epanechnikov'):
            self.__ktype = 0
        elif (ktype == 'tophat'):
            self.__ktype = 1
        self.leafsize = leafsize
        self.__abstol = abstol
        self.__reltol = reltol
        self.thread_safe = thread_safe

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
        self.__pobs = None
        self.__prior_data = None
        self.__pobs_data = None
        self.__kd_phys = None
        self.__cache = []

        # Build the initial kernel density estimation object, using a
        # dummy bandwidth; record mapping from indices of input data
        # to sorted indices, in case we need it to interpret arrays of
        # sample density, observation probability, or prior
        self.__bandwidth = np.ones(self.__nphys + self.__nphot)
        self.__idxmap = np.zeros(self.__ndata, dtype=c_ulong)
        self.__idxmap_inv = np.argsort(self.__idxmap)
        if nosort is None:
            self.__kd = self.__clib.build_kd(
                np.ravel(self.__dataset), self.__dataset.shape[1],
                self.__ndata, None, self.leafsize, self.__bandwidth,
                self.__ktype, 0,
                self.__idxmap.ctypes.data_as(POINTER(c_ulong)))
        else:
            nosort_c = np.zeros(nosort.shape, dtype=np.intc)
            nosort_c[:] = nosort == True
            self.__kd = self.__clib.build_kd_sortdims(
                np.ravel(self.__dataset), self.__dataset.shape[1],
                self.__ndata, None, self.leafsize, self.__bandwidth,
                self.__ktype, nosort_c,
                self.__idxmap.ctypes.data_as(POINTER(c_ulong)))

        # Store sample density
        if hasattr(sample_density, '__iter__'):
            self.__sden = np.array(sample_density)[self.__idxmap]
        else:
            self.__sden = sample_density
        self.__sample_density = None
        dummy = self.sample_density

        # Initialize the bandwidth
        self.bandwidth = bandwidth

        # Set priors and observation probability
        self.priors = priors
        self.pobs = pobs

        # Set the caching policy, and, if aggressive, construct the
        # data caches
        self.caching = caching
        if caching == 'lazy' and thread_safe:
            warn(
                "bp: cannot use lazy caching in thread safe"
                " mode; caching policy changed to 'none'")
            self.caching == 'none'
        if self.caching == 'aggressive':
            # Construct cached data sets for calculation of PDFs over
            # 1 dimension with all others marginalised out
            for d in range(self.ndim):
                margindims = list(range(0,d)) + \
                             list(range(d+1, self.ndim))
                self.make_cache(margindims)
            # Construct cached data sets for calculation of 1d
            # marginal PDFs of physical variables over photometric
            # ones
            if self.__nphys > 1:
                for d in range(self.__nphys):
                    margindims = list(range(0,d)) + \
                                 list(range(d+1, self.__nphys))
                    self.make_cache(margindims)


    ##################################################################
    # De-allocation method
    ##################################################################
    def __del__(self):
        if hasattr(self, '__kd'):
            if self.__kd is not None:
                self.__clib.free_kd(self.__kd)
        if hasattr(self, '__kd_phys'):
            if self.__kd_phys is not None:
                self.__clib.free_kd(self.__kd_phys)
        if hasattr(self, '__rng'):
            self.__clib.rng_free(self.__rng)


    ##################################################################
    # Return a copy of the list of available filters
    ##################################################################
    def filters(self):
        return deepcopy(self.__filters)
    
    ##################################################################
    # Property to get the sampling density
    ##################################################################

    @property
    def sample_density(self):
        """
        The density with which the library was sampled, evaluated for
        each simulation in the library
        """
        
        if self.__sample_density is None:

            # Choose computation method
            if hasattr(self.__sden, '__call__'):

                # Callable, so pass the physical data to the
                # callable and store the result
                self.__sample_density \
                    = self.__sden(self.__dataset[:,:self.__nphys])

            elif type(self.__sden) is np.ndarray or self.__sden is None:

                # Array, so treat treat this as the data
                self.__sample_density = self.__sden

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
                            None, self.leafsize, self.__bandwidth,
                            self.__ktype, 0)

                    # Use the unweighted kernel density object to
                    # evaluate the raw sample density near each data
                    # point, or near a sub-sample which we can
                    # interpolate from
                    self.__sample_density = np.zeros(self.__ndata)
                    nsamp = 500
                    if self.__ndata < nsamp:
                        # Few data points, just use them all; note
                        # that we cannot pass self.__dataset_phys,
                        # because it is not in the same order as
                        # the full data set anymore
                        pts = np.ravel(self.__dataset[:,:self.__nphys])
                        if not self.__diag_mode:
                            self.__clib.kd_pdf_vec(
                                self.__kd_phys, pts, None,
                                self.__ndata, self.reltol, self.abstol,
                                self.__sample_density)
                        else:
                            nodecheck = np.zeros(self.__ndata, dtype=c_ulong)
                            leafcheck = np.zeros(self.__ndata, dtype=c_ulong)
                            termcheck = np.zeros(self.__ndata, dtype=c_ulong)
                            self.__clib.kd_pdf_vec(
                                self.__kd_phys, pts, None,
                                self.__ndata, self.reltol, self.abstol,
                                self.__sample_density_sorted, nodecheck, 
                                leafcheck, termcheck)
                    else:
                        # Many data points, so choose a sample at random
                        idxpt = np.array(
                            random.sample(np.arange(self.__ndata), 
                                          nsamp), dtype=c_ulong)
                        pos = np.copy(self.__dataset_phys[idxpt,:])
                        # Add points at the edges of the dataset
                        # to ensure we enclose all the points
                        lowlim = np.amin(self.__dataset_phys,
                                         axis=0)
                        pos = np.append(pos, lowlim)
                        hilim = np.amax(self.__dataset_phys,
                                        axis=0)
                        pos = np.append(pos, hilim)
                        # Compute density at selected points
                        sample_density = np.zeros(nsamp+2)
                        if not self.__diag_mode:
                            self.__clib.kd_pdf_vec(
                                self.__kd_phys, np.ravel(pos), None,
                                nsamp+2, self.reltol, self.abstol,
                                sample_density)
                        else:
                            nodecheck = np.zeros(nsamp+2, dtype=c_ulong)
                            leafcheck = np.zeros(nsamp+2, dtype=c_ulong)
                            termcheck = np.zeros(nsamp+2, dtype=c_ulong)
                            self.__clib.kd_pdf_vec(
                                self.__kd_phys, np.ravel(pos), None,
                                nsamp+2, self.reltol, self.abstol,
                                sample_density, nodecheck, 
                                leafcheck, termcheck)
                        # Now interpolate the sample points to all
                        # points in the data set
                        pts = np.ravel(self.__dataset[:,:self.__nphys])
                        self.__sample_density \
                            = np.exp(
                                interp.griddata(pos, 
                                                np.log(sample_density),
                                                pts,
                                                method='linear')).flatten()
                        
        # Return result, sorted back into the order of the original
        # data
        if self.__sample_density is not None:
            return self.__sample_density[self.__idxmap_inv]
        else:
            return None
    

    @sample_density.setter
    def sample_density(self, sden):
        """
        This function sets the sampling density for the library

        Parameters
           sden : array, shape (N) | callable | None
              sample density each data point; interpretation depends
              on the type passed:
                 array, shape (N) : 
                    values are interpreted as the sampling density at
                    each data point in the dataset used to construct
                    the bp object
                 callable : 
                    the callable must take as an argument an array of
                    shape (N, nphys), and return an array of shape (N)
                    giving the prior probability at each data point
                 None :
                    all data points have equal prior probability

        Returns
           Nothing
        """
        # Erase current sample density information, then force
        # recomputation
        self.__sample_density = None
        if hasattr(sden, '__iter__'):
            # If passed an array, assume it is ordered as was the
            # original data set, not our re-ordered version, so
            # re-order before storing
            self.__sden = np.array(sden)[self.__idxmap]
        else:
            self.__sden = sden
        dummy = self.sample_density

        # Update any cached objects as well
        for c in self.__cache:
            c['bp'].sample_density = dummy
            

    ##################################################################
    # Define the priors and pobs properties
    ##################################################################

    @property
    def priors(self):
        """
        The current set of prior probabilities for every
        simulation in the library; data returned are in the same order
        as the data originally used to construct the library
        """
        if self.__priors is None:
            return None
        else:
            return self.__prior_data[self.__idxmap_inv]

    @priors.setter
    def priors(self, prior):
        """
        This function sets the prior probabilities to use

        Parameters
           prior : array, shape (N) | callable | None
              prior probability on each data point; interpretation
              depends on the type passed:
                 array, shape (N) : 
                    values are interpreted as the prior probability of
                    each data point in the dataset used to construct
                    the bp object
                 callable : 
                    the callable must take as an argument an array of
                    shape (N, nphys), and return an array of shape (N)
                    giving the prior probability at each data point
                 None :
                    all data points have equal prior probability

        Returns
           Nothing
        """

        # If prior is an arraylike, sort it using our index map
        if hasattr(prior, '__iter__'):
            pr = np.array(prior)[self.__idxmap]
        else:
            pr = prior

        # If the prior is unchanged, do nothing
        if (type(pr) == np.ndarray) and \
           (type(self.__priors) == np.ndarray):
            if np.array_equal(pr, self.__priors):
                return
        elif (type(pr) == type(self.__priors)) and \
             (pr == self.__priors):
            return

        # If priors and pobs are both None, just remove all weighting
        if pr is None and self.pobs is None:
            self.__clib.kd_change_wgt(None, self.__kd)
            self.__priors = None
            self.__prior_data = None
            for c in self.__cache:
                c['bp'].priors = None
            return

        else:
            # If we're here, we have a non-trival prior or pobs;
            # record the new prior, and, if it is a callable, call it
            self.__priors = pr
            if hasattr(self.__priors, '__call__'):
                self.__prior_data \
                    = self.__priors(self.__dataset[:,:self.__nphys]) \
                          .flatten()
            else:
                self.__prior_data = self.__priors

            # Compute the weights from the ratio of the prior times
            # the observation probability to the sample density, then
            # adjust the weights in the kd; note that a bit of
            # trickery is required, because the sample_density property
            # returns the sample density sorted in the original data
            # order, not in the sorted order, which is what we need
            dummy = self.sample_density    # Force re-computation
            wgt = np.ones(self.__ndata)
            if self.__sample_density is not None:
                wgt *= 1.0 / self.__sample_density
            if self.__prior_data is not None:
                wgt *= self.__prior_data
            if self.__pobs_data is not None:
                wgt *= self.__pobs_data
            self.__wgt = wgt
            self.__clib.kd_change_wgt(self.__wgt.ctypes.data_as(
                POINTER(c_double)), self.__kd)

            # Updated cached bp objects
            for c in self.__cache:
                c['bp'].priors = self.__prior_data


    @property
    def pobs(self):
        """
        The current set of observation probabilities for every
        simulation in the library; data returned are in the same order
        as the data originally used to construct the library
        """
        if self.__pobs_data is None:
            return None
        else:
            return self.__pobs_data[self.__idxmap_inv]

    @pobs.setter
    def pobs(self, p_obs):
        """
        This function sets the observation probabilities to use

        Parameters
           p_obs : array, shape (N) | callable | None
              probability of observation for each data point;
              interpretation depends on the type passed:
                 array, shape (N) : 
                    values are interpreted as the probability of
                    observation each data point in the dataset used to
                    construct the bp object
                 callable : 
                    the callable must take as an argument an array of
                    shape (N, nphot), and return an array of shape (N)
                    giving the probability of observation for each
                    data point
                 None :
                    all data points have equal probability of
                    observation

        Returns
           Nothing
        """

        # If prior is an arraylike, sort it using our index map
        if hasattr(p_obs, '__iter__'):
            po = np.array(p_obs)[self.__idxmap]
        else:
            po = p_obs

        # If the observation probability is unchanged, do nothing
        if (type(po) == np.ndarray) and \
           (type(self.__pobs) == np.ndarray):
            if np.array_equal(po, self.__pobs):
                return
        elif (type(po) == type(self.__pobs)) and \
             (po == self.__pobs):
            return

        # If pobs and prior are both None, just remove all weighting
        if po is None and self.priors is None:
            self.__clib.kd_change_wgt(None, self.__kd)
            self.__pobs = None
            self.__pobs_data = None
            for c in self.__cache:
                c['bp'].pobs = None
            return

        else:

            # If we're here, we have a non-trivial pobs or prior; if
            # the new pobs is a callable, call it and save the result
            self.__pobs = po
            if hasattr(self.__pobs, '__call__'):
                self.__pobs_data \
                    = self.__pobs(self.__dataset[:,self.__nphys:]).\
                    flatten()
            else:
                self.__pobs_data = self.__pobs

            # Compute the weights on the points, and change the kd
            # object appropriately
            dummy = self.sample_density    # Force re-computation
            wgt = np.ones(self.__ndata)
            if self.__sample_density is not None:
                wgt *= 1.0 / self.__sample_density
            if self.__prior_data is not None:
                wgt *= self.__prior_data
            if self.__pobs_data is not None:
                wgt *= self.__pobs_data
            self.__wgt = wgt
            self.__clib.kd_change_wgt(self.__wgt.ctypes.data_as(
                POINTER(c_double)), self.__kd)

            # Update cached objects
            for c in self.__cache:
                c['bp'].pobs = None


    ##################################################################
    # Define the bandwidth property
    ##################################################################
 
    @property
    def bandwidth(self):
        """
        The current bandwidth
        """
        return deepcopy(self.__bandwidth)

    @bandwidth.setter
    def bandwidth(self, bw):

        if np.array_equal(self.__bandwidth, bw):
            # If new bandwidth equals old bandwidth, do nothing
            return

        elif bw is not 'auto':
            # If we've been given a specified bandwidth, set to that
            if hasattr(bw, '__iter__'):
                self.__bandwidth = np.copy(bw)
            else:
                self.__bandwidth \
                    = np.zeros(self.__nphys + self.__nphot) + bw
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
                                      5000), dtype=c_ulong)
                    neighbors = np.zeros(nneighbor*5000, 
                                              dtype=c_ulong)
                    d2 = np.zeros(nneighbor*5000)
                    self.__clib.kd_neighbors_point_vec(
                        self.__kd, idxpt, 5000, nneighbor, False,
                        neighbors, d2)

                else:
                    # For smaller data sets, use it all
                    neighbors = np.zeros(nneighbor*self.__ndata, 
                                         dtype=c_ulong)
                    d2 = np.zeros(nneighbor*self.__ndata)
                    idxpt = np.arange(self.__ndata, dtype=c_ulong)
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
        if type(self.__sden) is str:
            if self.__sden == 'auto':
                self.__sample_density = None
                if self.__kd_phys is not None:
                    self.__clib.kd_change_bandwidth(
                        self.__bandwidth[:self.__nphys], self.__kd_phys)
                pr = self.priors
                self.priors = None
                self.priors = pr

        # Update cached objects
        for c in self.__cache:
            c['bp'] = bw

    ##################################################################
    # Utility methods to change the tolerances. We need properties for
    # these because, if we have cached data sets, we need to update
    # their tolerances too.
    ##################################################################
    @property
    def abstol(self):
        return self.__abstol

    @abstol.setter
    def abstol(self, newtol):
        self.__abstol = newtol
        for c in self.__cache:
            c['bp'].abstol = self.__abstol

    @property
    def reltol(self):
        return self.__reltol

    @reltol.setter
    def reltol(self, newtol):
        self.__reltol = newtol
        for c in self.__cache:
            c['bp'].reltol = self.__reltol

    ##################################################################
    # Utitlity methods to return the number of physical and
    # photometric quantities, and the total number of dimensions,
    # without allowing them to be changed
    ##################################################################

    @property
    def nphys(self):
        """
        This is the number of physical properties for the bayesphot
        object.
        """
        return self.__nphys

    @nphys.setter
    def nphys(self, val):
        raise ValueError("can't set bayesphot.nphys!")
    
    @property
    def nphot(self):
        """
        This is the number of photometric properties for the bayesphot
        object.
        """
        return self.__nphot

    @nphot.setter
    def nphot(self, val):
        raise ValueError("can't set bayesphot.nphot!")

    @property
    def ndim(self):
        """
        This is the number of physical plus photometric properties for
        the bayesphot object.
        """
        return (self.__nphys + self.__nphot)

    @ndim.setter
    def ndim(self, val):
        raise ValueError("can't set bayesphot.ndim!")

    ##################################################################
    # Utility methods to set a new bandwidth to account for
    # photometric errors, to reset to the original bandwidth, and to
    # free a temporary kd_tree holder, all in a thread-safe way. These
    # are intended for internal use.
    ##################################################################

    # Broaden the bandwidth by the photometric error, and return a KD
    # object with the new bandwidth; if we are thread-safe, this is a
    # newly-created KD object, which will need to be freed later; if a
    # value for kd_cur is passed in, the bandwidth will be changed for
    # it, and no allocation will be done regardless of thread safety
    def __change_bw_err(self, photerr, kd_cur=None, err_only=False):
        if err_only:
            bandwidth = np.copy(self.__bandwidth)
            bandwidth[self.__nphys:] = photerr
        else:
            err = np.zeros(self.__bandwidth.size)
            err[self.__nphys:] = photerr
            bandwidth = np.sqrt(self.__bandwidth**2+err**2)
        if kd_cur is not None:
            kd_tmp = kd_cur
        elif not self.thread_safe: 
            kd_tmp = self.__kd
        else: 
            kd_tmp = self.__clib.copy_kd(self.__kd)
        self.__clib.kd_change_bandwidth(bandwidth, kd_tmp)
        return kd_tmp

    # Free a temporary KD object
    def __free_bw_err(self, kd_tmp):
        if self.thread_safe:
            self.__clib.free_kd_copy(kd_tmp)

    # Restore bandwidth to the default, de-allocating the temporary KD
    # object if needed
    def __restore_bw_err(self, kd_tmp):
        if not self.thread_safe:
            self.__clib.kd_change_bandwidth(
                self.__bandwidth, self.__kd)
        else:
            self.__clib.free_kd_copy(kd_tmp)

    ##################################################################
    # Methods to create and destroy caches
    ##################################################################
    def make_cache(self, margindims):
        """
        This method builds a cache to do faster calculation of PDFs
        where certain dimensions are marginalised out. If such caches
        exist, they are used automatically by all the computation
        methods.

        Parameters
           margindims : listlike of integers
              list of dimensions to be marginalised out; physical
              dimensions go from 0 - nphys-1, photometric dimensions
              from nphys to nphys + nphot - 1

        Returns
           Nothing
        """

        # Safety check on inputs
        mdims = np.asarray(margindims, dtype='int')
        if np.amin(mdims) < 0 or np.amin(mdims) > self.ndim-1 or \
           len(mdims) != len(np.unique(mdims)):
            raise ValueError("bp.make_cache: margindims must be "
                             " unique integers in the range "
                             " 0, ndim")

        # If we already have a cache for this set of marginalised
        # dimensions, do nothing
        for c in self.__cache:
            if c['margindims'] == list(margindims):
                return

        # Get list of dimensions that remain
        keepdims = np.array(list(set(range(self.ndim)) - set(mdims)),
                            dtype=int)
        nphys = np.sum(keepdims < self.nphys)

        # Extract the dimensions we're keeping from the data set
        data_cache = np.copy(self.__dataset[:,keepdims])
        bw = self.bandwidth[keepdims]
        if self.__prior_data is not None:
            prior_cache = np.copy(self.__prior_data)
        else:
            prior_cache = None
        if self.__pobs_data is not None:
            pobs_cache = np.copy(self.__pobs_data)
        else:
            pobs_cache = None

        # Build a new bp object from the dimensionally-reduced data
        if self.__ktype == 2:
            ktype = 'gaussian'
        elif self.__ktype == 1:
            ktype = 'tophat'
        elif self.__ktype == 0:
            ktype = 'epanechnikov'
        bp_cache = bp(data_cache, nphys, bandwidth=bw, ktype=ktype,
                      priors=prior_cache, pobs=pobs_cache,
                      sample_density=self.__sample_density,
                      reltol=self.reltol,
                      abstol=self.abstol, leafsize=self.leafsize,
                      thread_safe=self.thread_safe,
                      caching='none')

        # Store the new object
        self.__cache.append(
            { 'margindims' : list(margindims),
              'keepdims'   : list(keepdims),
              'data'       : data_cache,
              'prior'      : prior_cache,
              'pobs'       : pobs_cache,
              'bp'         : bp_cache } )

    def clear_cache(self, margindims=None):
        """
        This method deletes from the cache

        Parameters
           margindims : listlike of integers
              list of marginalised dimensions that should be removed
              from the cache, in the same format as make_cache; if
              left as None, the cache is completely emptied

        Returns
           Nothing
        """
        if margindims is None:
            self.__cache = []
        else:
            for i in range(len(cache)):
                if cache[i]['margindims'] == list(margindims):
                    del cache[i]
                    return
            raise ValueError(
                "bp.clear_cache: no cached data found for "+
                repr(margindims))
        

    ##################################################################
    # Method to compute the log likelihood function for a particular
    # set of physical properties given a particular set of photometric
    # properties; some dimensions can be marginalised over if desired
    ##################################################################
    def logL(self, physprop, photprop, photerr=None,
             margindim=None):
        """
        This function returns the natural log of the likelihood
        function evaluated at a particular log mass, log age,
        extinction, and set of log luminosities

        Parameters
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

        Returns
           logL : float or arraylike
              natural log of the likelihood function
        """

        # If we're marginalising, see if we should use a cached data
        # set
        if margindim is not None:

            # If we're in lazy mode, add a cache for this set of
            # marginalised dimensions if it doesn't exist
            if self.caching == 'lazy':
                self.make_cache(margindim)

            # Look for this set of dimensions in our existing cache;
            # if we find a match, use it
            for c in self.__cache:
                if list(margindim) == c['margindims']:
                    return c['bp'].logL(physprop, photprop,
                                        photerr=photerr)

        # Safety check: make sure all the dimensions of the input
        # arrays match what we expect. If we are marginalizing over
        # some dimensions, extract the dimensions we want to keep
        if margindim is not None:
            nphys_margin = np.sum(np.asarray(margindim) <
                                  self.__nphys)
            nphot_margin = len(margindim) - nphys_margin
            nphys = self.__nphys - nphys_margin
            nphot = self.__nphot - nphot_margin
            keepdims = np.array(
                list(set(range(self.ndim)) - set(margindim)),
                dtype=c_ulong)
            bw = self.__bandwidth[keepdims]
        else:
            nphys = self.__nphys
            nphot = self.__nphot
            bw = self.__bandwidth
        if nphys > 1:
            if (np.asarray(physprop).shape[-1] != nphys):
                raise ValueError("need " + str(nphys) + 
                                 " physical properties!")
        if nphys == 0 and physprop is not None:
            raise ValueError(
                "input includes physical properties, but "
                "all physical properties are marginalised "
                "out!")
                
        if nphot > 1:
            if (np.asarray(photprop).shape[-1] != nphot):
                raise ValueError("need " + str(nphot) +
                                 " photometric properties!")
            if photerr is not None:
                if (np.asarray(photerr).shape[-1] != nphot):
                    raise ValueError("need " + str(nphot) +
                                     " photometric errors!")
        if nphot == 0 and photprop is not None:
            raise ValueError(
                "input includes photometric properties, but "
                "all photometric properties are marginalised "
                "out!")
        if nphot == 0 and photerr is not None:
            raise ValueError(
                "input includes photometric errors, but "
                "all photometric properties are marginalised "
                "out!")

        # Handle special case of a single photometric or physical
        # quantity; in this case make sure we have a trailing
        # dimension of 1 so that all the fancy array manipulation
        # below works correctly
        if nphys == 1 and (np.asarray(physprop).shape[-1] != 1
                           or np.asarray(physprop).ndim == 1):
            physprop_ = np.array(physprop).\
                        reshape(np.asarray(physprop).shape+(1,))
        else:
            physprop_ = physprop
        if nphot == 1 and (np.asarray(photprop).shape[-1] != 1
                           or np.asarray(photprop).ndim == 1):
            photprop_ = np.array(photprop).\
                        reshape(np.asarray(photprop).shape+(1,))
        else:
            photprop_ = photprop
        if nphot == 1 and photerr is not None:
            if np.asarray(photerr).shape[-1] != 1 \
               or np.asarray(photerr).ndim == 1:
                photerr_ = np.array(photerr).\
                        reshape(np.asarray(photerr).shape+(1,))
        else:
            photerr_ = photerr

        # Figure out number of distinct input sets of physical and
        # photometric properties in the input data
        if nphys > 0:
            nphys_in = np.asarray(physprop_).size // nphys
        else:
            nphys_in = 0
        if nphot > 0:
            nphot_in = np.asarray(photprop_).size // nphot
        else:
            nphot_in = 0

        # Allocate an array to hold the results. The shape is a bit
        # tricky to figure out, since one or more of the inputs
        # (physprop, photprop, and photerr) can be None. We therefore
        # need to consider several possible cases.
        if physprop_ is not None and photprop_ is not None and \
           photerr_ is not None:
            pdf = np.zeros(
                np.broadcast(np.atleast_2d(physprop_)[..., 0],
                             np.atleast_2d(photprop_)[..., 0],
                             np.atleast_2d(photerr_)[..., 0]).shape,
                dtype=c_double)
        elif physprop_ is not None and photprop_ is not None:
            pdf = np.zeros(
                np.broadcast(np.atleast_2d(physprop_)[..., 0],
                             np.atleast_2d(photprop_)[..., 0]).shape,
                dtype=c_double)
        elif physprop_ is not None:
            pdf = np.zeros(np.atleast_2d(physprop_)[..., 0].shape,
                           dtype=c_double)
        elif photprop_ is not None and photerr_ is not None:
            pdf = np.zeros(
                np.broadcast(np.atleast_2d(photprop_)[..., 0],
                             np.atleast_2d(photerr_)[..., 0]).shape,
                dtype=c_double)
        elif photprop_ is not None:
            pdf = np.zeros(np.atleast_2d(photprop_)[..., 0].shape,
                           dtype=c_double)            
        if self.__diag_mode:
            nodecheck = np.zeros(pdf.shape, dtype=c_ulong)
            leafcheck = np.zeros(pdf.shape, dtype=c_ulong)
            termcheck = np.zeros(pdf.shape, dtype=c_ulong)

        # Make an array suitable for passing data to c routines
        cdata = np.zeros(pdf.shape + (nphys+nphot,),
                         dtype=c_double)
        if nphys != 0:
            cdata[..., :nphys] \
                = np.vstack((physprop_,) * (pdf.size//nphys_in)). \
                reshape(cdata[..., :nphys].shape)
        if nphot != 0:
            cdata[..., nphys:] \
                = np.vstack((photprop_,) * (pdf.size//nphot_in)). \
                reshape(cdata[..., nphys:].shape)
        if photerr_ is not None:
            cbandwidth = np.zeros(pdf.shape + (nphot,),
                                  dtype=c_double)                
            cbandwidth[...] = np.sqrt(bw**2 + photerr_**2)
            cbw_ptr = cbandwidth.ctypes.data_as(POINTER(c_double))
        else:
            cbw_ptr = None
            
        # Call the PDF computation routine
        if margindim is None:
            # No marginalization, so call kd_pdf

            if not self.__diag_mode:
                self.__clib.kd_pdf_vec(
                    self.__kd, np.ravel(cdata), cbw_ptr, pdf.size, 
                    self.reltol, self.abstol, np.ravel(pdf))
            else:
                self.__clib.kd_pdf_vec(
                    self.__kd, np.ravel(cdata), cbw_ptr, pdf.size, 
                    self.reltol, self.abstol, np.ravel(pdf),
                    np.ravel(nodecheck), np.ravel(leafcheck),
                    np.ravel(termcheck))
        else:
            # We are marginalizing, so call kd_pdf_int
            if not self.__diag_mode:
                self.__clib.kd_pdf_int_vec(
                    self.__kd, np.ravel(cdata), cbw_ptr, keepdims,
                    keepdims.size, pdf.size, self.reltol,
                    self.abstol, np.ravel(pdf))
            else:
                self.__clib.kd_pdf_int_vec(
                    self.__kd, np.ravel(cdata), cbw_ptr, keepdims,
                    keepdims.size, pdf.size, self.reltol,
                    self.abstol, np.ravel(pdf),
                    np.ravel(nodecheck), np.ravel(leafcheck),
                    np.ravel(termcheck))

        # Return
        if not self.__diag_mode:
            return np.log(pdf)
        else:
            return np.log(pdf), nodecheck, leafcheck, termcheck

    ##################################################################
    # Function to return the marginal distribution of one or more of
    # the physical properties for a specified set of photometric
    # properties
    ##################################################################
    def mpdf(self, idx, photprop, photerr=None, ngrid=128,
             qmin=None, qmax=None, grid=None, norm=True):
        """
        Returns the marginal probability for one or mode physical
        quantities for one or more input sets of photometric
        properties. Output quantities are computed on a grid of
        values, in the same style as meshgrid.

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
              of the output grid, for each input set of photometry; the leading
              dimensions match the leading dimensions produced by
              broadcasting the leading dimensions of photprop and
              photerr together, while the trailing dimensions match
              the dimensions of the output grid
        """

        # Check if we're marginalising over any dimensions, and if so
        # whether we have or should create a cache for those
        # dimensions
        if len(np.atleast_1d(idx)) < self.nphys:

            # Figure out which dimensions we're marginalising over
            keepdim = list(np.atleast_1d(idx))
            margindim = list(set(range(self.nphys)) - set(keepdim))

            # If we're in lazy mode, add a cache for this set of
            # marginalised dimensions
            if self.caching == 'lazy':
                self.make_cache(margindim)
            
            # Look for this set of dimensions in our existing cache
            for c in self.__cache:
                if list(margindim) == c['margindims']:

                    # Found a match; call the bp object that matches
                    # and return its output
                    return c['bp'].mpdf(keepdim, photprop,
                                        photerr=photerr,
                                        ngrid=ngrid,
                                        qmin=qmin, qmax=qmax,
                                        grid=grid, norm=norm)
                
        # Safety check
        if (np.array(photprop).shape[-1] != self.__nphot) and \
           (self.__nphot > 1):
            raise ValueError("need " + str(self.__nphot) +
                             " photometric properties!")
        if photerr is not None:
            if (np.array(photerr).shape[-1] != self.__nphot) and \
               (self.__nphot > 1):
                raise ValueError("need " + str(self.__nphot) +
                                 " photometric errors!")
        if (np.amax(idx) > self.nphys) or (np.amin(idx) < 0) or \
           (not np.array_equal(np.squeeze(np.unique(np.array(idx))), 
                               np.squeeze(np.array([idx])))):
            raise ValueError("need non-repeating indices in " +
                             "the range 0 - {:d}!".
                             format(nphys-1))

        # Reshape arrays if necessary
        if (np.array(photprop).shape[-1] != self.__nphot) and \
           (self.__nphot == 1):
            photprop = photprop.reshape(photprop.shape+(1,))
        if photerr is not None:
            if (np.array(photerr).shape[-1] != self.__nphot) and \
               (self.__nphot == 1):
                photerr = photerr.reshape(photerr.shape+(1,))

        # Set up the grid of outputs, and the versions of it that we
        # will be passing to the c code
        if grid is not None:
            grid_out = np.array(grid)
            qmin_tmp \
                = np.array(
                    np.copy(grid_out[(Ellipsis,)+(0,)*(grid_out.ndim-1)]),
                    dtype=np.double)
            qmax_tmp \
                = np.array(
                    np.copy(grid_out[(Ellipsis,)+(-1,)*(grid_out.ndim-1)]),
                    dtype=np.double)
            ngrid_tmp = np.array(grid_out.shape[1:], dtype=c_ulong)
        else:
            if qmin is None:
                qmin = np.amin(self.__dataset[:,idx], axis=0)
            if qmax is None:
                qmax = np.amax(self.__dataset[:,idx], axis=0)
            griddims = []
            if hasattr(idx, '__len__'):
                nidx = len(idx)
            else:
                nidx = 1
            if nidx > 1:
                # Case for multiple indices
                griddims = []
                if hasattr(ngrid, '__len__'):
                    ngrid_tmp = np.array(ngrid, dtype=c_ulong)
                else:
                    ngrid_tmp = np.array([ngrid]*len(idx), dtype=c_ulong)
                for i in range(len(idx)):
                    griddims.append(qmin[i] + np.arange(ngrid_tmp[i]) * 
                                    float(qmax[i]-qmin[i])/(ngrid_tmp[i]-1))
                grid_out = np.squeeze(np.array(np.meshgrid(*griddims,
                                                           indexing='ij')))
                out_shape = grid_out[0, ...].shape
                qmin_tmp \
                    = np.copy(grid_out[(Ellipsis,)+(0,)*(grid_out.ndim-1)])
                qmax_tmp \
                    = np.copy(grid_out[(Ellipsis,)+(-1,)*(grid_out.ndim-1)])
            else:
                # Case for a single index
                ngrid_tmp = np.array([ngrid], dtype=c_ulong)
                qmin_tmp = np.array([qmin], dtype=np.double).reshape(1)
                qmax_tmp = np.array([qmax], dtype=np.double).reshape(1)
                grid_out = qmin + \
                           np.arange(ngrid) * \
                           float(qmax-qmin)/(ngrid-1)
                out_shape = grid_out.shape

        # Figure out how many distinct photometric values we've been
        # given, and how many sets of photometric errors
        nphot = np.array(photprop).size//self.__nphot
        nphot_err = np.array(photerr).size//self.__nphot

        # Set up a grid to hold the outputs
        if photerr is not None:
            pdf = np.zeros(
                np.broadcast(np.array(photprop)[..., 0],
                             np.array(photerr)[..., 0]).shape +
                out_shape)
        else:
            pdf = np.zeros(np.array(photprop)[..., 0].shape +
                           out_shape)

        # Prepare data for c library
        if hasattr(idx, '__len__'):
            nidx = len(idx)
        else:
            nidx = 1
        dims = np.zeros(nidx+self.__nphot, dtype=c_ulong)
        dims[:nidx] = idx
        dims[nidx:] = self.__nphys+np.arange(self.__nphot, dtype=c_ulong)
        ndim = c_ulong(nidx + self.__nphot)
        phottmp = np.array(photprop, dtype=c_double)

        # Separate cases with single / no photometric errors from
        # cases with multiple sets of photometric errors
        if nphot_err <= 1:

            # Case with at most one set of photometric errors

            # Set the bandwidth based on the photometric errors if we
            # were given some
            if photerr is None:
                kd_tmp = self.__kd
            else:
                kd_tmp = self.__change_bw_err(photerr)

            # Call the PDF computation routine; note that we call
            # kd_pdf_int_vec if we're actually marginalizing over any
            # physical properties (which is the case if len(idx) <
            # nphys), but we invoke kd_pdf_vec if we're not actually
            # marginalizing (len(idx)==nphys) because then we don't
            # need to do any integration
            if nidx < self.__nphys:
                self.__clib.kd_pdf_int_reggrid(
                    kd_tmp, np.ravel(phottmp),
                    dims[nidx:], self.__nphot, nphot,
                    qmin_tmp, qmax_tmp, ngrid_tmp, dims[:nidx], nidx,
                    self.reltol, self.abstol, np.ravel(pdf))
            else:
                self.__clib.kd_pdf_reggrid(
                    kd_tmp, np.ravel(phottmp),
                    dims[nidx:], self.__nphot, nphot,
                    qmin_tmp, qmax_tmp, ngrid_tmp, dims[:nidx], nidx,
                    self.reltol, self.abstol, np.ravel(pdf))

        else:

            # Case with multiple sets of photometric errors

            # Loop over photometric errors
            kd_tmp = None
            for i in np.ndindex(*photerr.shape[:-1]):

                # Set bandwidth based on photometric error for this
                # iteration
                kd_tmp = self.__change_bw_err(photerr[i], kd_cur=kd_tmp)

                # Grab the corresponding portions of the arrays going
                # to and from the c code
                if photprop.size < photerr.size:
                    phot_sub = np.ravel(phottmp)
                else:
                    phot_sub = np.ravel(phottmp[i])
                pdf_sub = np.zeros(np.array(pdf[i]).shape)

                # Call kernel density estimate with this bandwidth
                if nidx < self.__nphys:
                    self.__clib.kd_pdf_int_reggrid(
                        kd_tmp, phot_sub,
                        dims[nidx:], self.__nphot, 
                        phot_sub.size//self.__nphot,
                        qmin_tmp, qmax_tmp, ngrid_tmp, 
                        dims[:nidx], nidx,
                        self.reltol, self.abstol, np.ravel(pdf_sub))
                else:
                    self.__clib.kd_pdf_reggrid(
                        kd_tmp, phot_sub,
                        dims[nidx:], self.__nphot, 
                        phot_sub.size//self.__nphot,
                        qmin_tmp, qmax_tmp, ngrid_tmp, dims[:nidx], nidx,
                        self.reltol, self.abstol, np.ravel(pdf_sub))
                pdf[i] = pdf_sub

        # Set the bandwidth back to its default if necessary
        if photerr is not None:
            self.__restore_bw_err(kd_tmp)

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
                             tuple(range(np.array(photprop).ndim-1, pdf.ndim)))

            # Normalize
            pdf = np.transpose(np.transpose(pdf)/normfac)

        # Return
        return grid_out, pdf


    ##################################################################
    # Function to return the marginal distribution of one or more of
    # the photometric properties for a specified set of physical
    # properties; this is just like mpdf, but in reverse
    ##################################################################
    def mpdf_phot(self, idx, physprop, ngrid=128,
                  qmin=None, qmax=None, grid=None, norm=True):
        """
        Returns the marginal probability for one or mode photometric
        quantities corresponding to an input set of physical
        properties. Output quantities are computed on a grid of
        values, in the same style as meshgrid.

        Parameters:
           idx : int or listlike containing ints
              index of the photometric quantity whose PDF is to be
              computed, starting at 0; if this is an iterable, the
              joint distribution of the indicated quantities is returned
           physprop : arraylike, shape (nphys) or (..., nphys)
              physical properties to be used; if this is an array of
              nphys elements, these give the physical properties; if
              it is a multidimensional array, the operation is
              vectorized over the leading dimensions
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

        # Check if we're marginalising over any dimensions, and if so
        # whether we have or should create a cache for those
        # dimensions
        if len(np.atleast_1d(idx)) < self.nphot:
            
            # Figure out which dimensions we're marginalising over
            keepdim = list(np.atleast_1d(idx)+self.nphys)
            margindim = list(set(range(self.nphys,self.ndim))
                             - set(keepdim))
            
            # If we're in lazy mode, add a cache for this set of
            # marginalised dimensions
            if self.caching == 'lazy':
                self.make_cache(margindim)
            
            # Look for this set of dimensions in our existing cache
            for c in self.__cache:
                if list(margindim) == c['margindims']:

                    # Found a match; call the bp object that matches
                    # and return its output
                    return c['bp'].mpdf_phot(
                        np.array(keepdim)-self.nphys,
                        physprop,
                        ngrid=ngrid,
                        qmin=qmin, qmax=qmax,
                        grid=grid, norm=norm)
                
        # Safety check
        if (np.array(physprop).shape[-1] != self.__nphys) and \
           (self.__nphot > 1):
            raise ValueError("need " + str(self.__nphys) +
                             " physical properties!")
        if (np.amax(np.asarray(idx)) > self.__nphot) or \
           (np.amin(np.asarray(idx)) < 0) or \
           (not np.array_equal(np.squeeze(np.unique(np.asarray(idx))), 
                               np.squeeze(np.array([idx])))):
            raise ValueError("need non-repeating indices in " +
                             "the range 0 - {:d}!".
                             format(self.__nphot-1))

        # Reshape arrays if necessary
        if (np.array(physprop).shape[-1] != self.__nphys) and \
           (self.__nphys == 1):
            physprop = physprop.reshape(physprop.shape+(1,))

        # Set up the grid of outputs, and the versions of it that we
        # will be passing to the c code
        if grid is not None:
            grid_out = np.array(grid)
            qmin_tmp \
                = np.array(
                    np.copy(grid_out[(Ellipsis,)+(0,)*(grid_out.ndim-1)]),
                    dtype=np.double)
            qmax_tmp \
                = np.array(
                    np.copy(grid_out[(Ellipsis,)+(-1,)*(grid_out.ndim-1)]),
                    dtype=np.double)
            ngrid_tmp = np.array(grid_out.shape[1:], dtype=c_ulong)
        else:
            if qmin is None:
                qmin = np.amin(self.__dataset[:,idx+self.__nphys], axis=0)
            if qmax is None:
                qmax = np.amax(self.__dataset[:,idx+self.__nphys], axis=0)
            griddims = []
            if hasattr(idx, '__len__'):
                nidx = len(idx)
            else:
                nidx = 1
            if nidx > 1:
                # Case for multiple indices
                griddims = []
                if hasattr(ngrid, '__len__'):
                    ngrid_tmp = np.array(ngrid, dtype=c_ulong)
                else:
                    ngrid_tmp = np.array([ngrid]*len(idx), dtype=c_ulong)
                for i in range(len(idx)):
                    griddims.append(qmin[i] + np.arange(ngrid_tmp[i]) * 
                                    float(qmax[i]-qmin[i])/(ngrid_tmp[i]-1))
                grid_out = np.squeeze(np.array(np.meshgrid(*griddims,
                                                           indexing='ij')))
                out_shape = grid_out[0, ...].shape
                qmin_tmp \
                    = np.copy(grid_out[(Ellipsis,)+(0,)*(grid_out.ndim-1)])
                qmax_tmp \
                    = np.copy(grid_out[(Ellipsis,)+(-1,)*(grid_out.ndim-1)])
            else:
                # Case for a single index
                ngrid_tmp = np.array([ngrid], dtype=c_ulong)
                qmin_tmp = np.array([qmin], dtype=np.double).reshape(1)
                qmax_tmp = np.array([qmax], dtype=np.double).reshape(1)
                grid_out = qmin + \
                           np.arange(ngrid) * \
                           float(qmax-qmin)/(ngrid-1)
                out_shape = grid_out.shape
        

        # Prepare things we'll be passing to the c library
        # Prepare data for c library
        if hasattr(idx, '__len__'):
            nidx = len(idx)
        else:
            nidx = 1
        dims = np.zeros(nidx+self.__nphys, dtype=c_ulong)
        dims[:nidx] = self.__nphys + idx
        dims[nidx:] = np.arange(self.__nphys, dtype=c_ulong)
        ndim = c_ulong(nidx + self.__nphys)

        # Prepare output holder
        pdf = np.zeros(np.array(physprop)[..., 0].shape +
                       out_shape)

        # Prepare inputs for passing to c
        nphys = np.array(physprop).size//self.__nphys
        phystmp = np.array(physprop, dtype=c_double)

        # Call c library; note that we call kd_pdf_int_reggrid if
        # we're actually marginalizing over any photometric
        # properties (which is the case if len(idx) <
        # self.__nphot), but we invoke kd_pdf_reggrid if we're not
        # actually marginalizing (len(idx)==self.__nphot) because
        # then we don't need to do any integration
        if nidx < self.__nphot:
            self.__clib.kd_pdf_int_reggrid(
                self.__kd, np.ravel(phystmp), dims[nidx:],
                self.__nphys, nphys, qmin_tmp, qmax_tmp,
                ngrid_tmp, dims[:nidx], nidx, self.reltol,
                self.abstol, np.ravel(pdf))
        else:
            self.__clib.kd_pdf_reggrid(
                self.__kd, np.ravel(phystmp), dims[nidx:],
                self.__nphys, nphys, qmin_tmp, qmax_tmp,
                ngrid_tmp, dims[:nidx], nidx,
                self.reltol, self.abstol, np.ravel(pdf))

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
                             tuple(range(np.array(physprop).ndim-1, pdf.ndim)))

            # Normalize
            pdf = np.transpose(np.transpose(pdf)/normfac)

        # Return
        return grid_out, pdf
                

    ##################################################################
    # This is the most general verison of mpdf. It allows the user to
    # fix or marginalize any number of physical or photometric
    # properties, returning the distribution in the remaining
    # properties
    ##################################################################
    def mpdf_gen(self, fixeddim, fixedprop, margindim, ngrid=128,
                 qmin=None, qmax=None, grid=None, norm=True):
        """
        Returns the marginal probability for one or more physical or
        photometric properties, keeping other properties fixed and
        marginalizing over other quantities. This is the most general
        marginal PDF routine provided.

        Parameters:
           fixeddim : int | arraylike of ints | None
              The index or indices of the physical or photometric
              properties to be held fixed; physical properties are
              numbered 0 ... nphys-1, and phtometric ones are numbered
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

        # If we're marginalising out some dimensions, and we have
        # cached data for that setup or are going to make it, call our
        # cached bp object to do this work for us
        if margindim is not None:

            # If we're in lazy mode, add a cache for this set of
            # marginalised dimensions
            if self.caching == 'lazy':
                self.make_cache(margindim)

            # Look for this set of dimensions in our existing cache
            for c in self.__cache:
                if list(margindim) == c['margindims']:

                    # Found a match; update the dimensional indexing
                    # on fixeddim, then call
                    if fixeddim is None:
                        fixeddim_ = None
                    else:
                        fixeddim_ = np.array(deepcopy(fixeddim))
                        diff = np.zeros(fixeddim_.size, dtype=np.int)
                        for d in margindim:
                            diff[fixeddim_ > d] += 1
                        fixeddim_ -= diff
                    return c['bp'].mpdf_gen(fixeddim_, fixedprop,
                                            None, ngrid=ngrid,
                                            qmin=qmin, qmax=qmax,
                                            grid=grid, norm=norm)

        # Safety check: make sure the dimensions all match up
        if fixeddim is None:
            nfixed = 0
        else:
            nfixed = np.atleast_1d(np.array(fixeddim)).shape[0]
            if len(set(np.atleast_1d(fixeddim))) != nfixed:
                raise ValueError("fixeddim must not have any "
                                 "repeated elements!")
            if np.array(fixedprop).shape[-1] != nfixed:
                raise ValueError("final dimension of fixedprop must" +
                                 " have {:d} dimensions".format(nfixed) +
                                 " to match fixeddim!")
        if margindim is None:
            nmargin = 0
        else:
            nmargin = np.atleast_1d(np.array(margindim)).shape[0]
            if len(set(np.atleast_1d(margindim))) != nmargin:
                raise ValueError("margindim must not have any "
                                 "repeated elements!")
        if qmin is not None:
            nq = np.atleast_1d(np.array(qmin)).shape[0]
            if nq + nmargin + nfixed != \
               self.__nphys + self.__nphot:
                raise ValueError("qmin must have {:d}".
                                 format(self.__nphys + self.__nphot -
                                        nfixed - nmargin) +
                                 " dimensions!")        
        else:
            nq = self.__nphys + self.__nphot - nfixed - nmargin
        if qmax is not None:
            if nq != np.atleast_1d(np.array(qmax)).shape[0]:
                raise ValueError("qmin and qmax must have matching" +
                                 " number of dimensions")
        if grid is not None:
            nq = len(grid)
            if nq + nmargin + nfixed != \
               self.__nphys + self.__nphot:
                raise ValueError("grid must have {:d}".
                                 format(self.__nphys + self.__nphot -
                                        nfixed - nmargin) +
                                 " elements!")
        if fixeddim is not None and margindim is not None:
            if len(set(margindim) & set(fixeddim)) > 0:
                raise ValueError("margindim and fixeddim must not "
                                 "have any elements in common!")

        # Make versions of all inputs arrays suitable for passing to c
        if margindim is not None:
            margindim_tmp = np.atleast_1d(
                np.array(margindim, dtype=c_ulong))
        if fixeddim is not None:
            fixeddim_tmp = np.atleast_1d(
                np.array(fixeddim, dtype=c_ulong))
            xfixed = np.atleast_1d(np.array(fixedprop, dtype=c_double))
            
        # Figure out the indices we're returning
        idx = set(range(self.__nphys+self.__nphot))
        if margindim is not None:
            idx = idx - set(margindim)
        if fixeddim is not None:
            idx = idx - set(fixeddim)
        idx = np.array(list(idx), dtype=c_ulong)
            
        # Set up output grid; make a version to pass to the c code
        # that has the axes tranposed, since the c code needs a
        # different data ordering
        nidx, qmin_tmp, qmax_tmp, grid_out, out_shape \
            = self.__make_outgrid(idx, ngrid, qmin, qmax, grid)
        grid_out_c = np.transpose(grid_out)

        # Set up array to hold the final pdf
        if fixedprop is not None:
            pdf = np.zeros(np.array(fixedprop)[..., 0] + out_shape,
                           dtype=c_double)
        else:
            pdf = np.zeros(out_shape, dtype=c_double)

        # Call the correct c library routine, depending on whether we
        # have fixed points, marginalized dimensions, both, or neither
        if fixeddim is None:

            if margindim is None:

                # No fixed points, no marginalized dimensions, so just
                # call kd_pdf_vec
                self.__clib.kd_pdf_vec(
                    self.__kd, np.ravel(grid_out_c), None,
                    grid_out_c.size//self.ndim,
                    self.reltol, self.abstol,
                    np.ravel(pdf))

            else:

                # Marginalizing over some dimensions, but no fixed
                # points, so call kd_pdf_int_vec
                self.__clib.kd_pdf_int_vec(
                    self.__kd, np.ravel(grid_out_c), None,
                    idx, len(idx), grid_out_c.size//len(idx),
                    self.reltol, self.abstol, np.ravel(pdf))

        else:

            if margindim is None:

                # We have fixed points, but no marginalized
                # dimensions, so call kd_pdf_reggrid
                self.__clib.kd_pdf_reggrid(
                    self.__kd, np.ravel(xfixed), fixeddim_tmp,
                    nfixed, xfixed.size//nfixed,
                    qmin_tmp, qmax_tmp, ngrid_tmp, idx,
                    idx.size, self.reltol, self.abstol,
                    np.ravel(pdf))

            else:

                # We have both fixed points and marginalized
                # dimensions, so call kd_pdf_int_reggrid
                self.__clib.kd_pdf_int_reggrid(
                    self.__kd, np.ravel(xfixed), fixeddim_tmp,
                    nfixed, xfixed.size//nfixed, qmin_tmp,
                    qmax_tmp, idx, idx.size, self.reltol,
                    self.abstol, np.ravel(pdf))

        # Normalize if requested
        if norm:
            pdf = self.__mpdf_normalize(nidx, grid_out, pdf)

        # Return
        return grid_out, pdf
    
    
    ##################################################################
    # Utility method to set up grids for outputs of marginal PDF
    # routines that are appropriate for passing to the c library. Used
    # by all the mpdf routines.
    ##################################################################
    def __make_outgrid(self, idx, ngrid, qmin, qmax, grid):

        # Were we given a grid?
        if grid is not None:
            # Yes, so just make it an array, and extract the limits
            # from it
            grid_out = np.array(grid, dtype=c_double)
            out_shape = grid_out.shape
            qmin_tmp \
                = np.array(
                    np.copy(grid_out[(Ellipsis,)+(0,)*(grid_out.ndim-1)]),
                    dtype=np.double)
            qmax_tmp \
                = np.array(
                    np.copy(grid_out[(Ellipsis,)+(-1,)*(grid_out.ndim-1)]),
                    dtype=np.double)
            ngrid_tmp = np.array(grid_out.shape[1:], dtype=c_ulong)
        else:
            # No input grid was given

            # Were we given upper or lower limits? If not, read them
            # from the library. For photometric dimensions, exclude
            # values of 99, which indicate bad data.
            if qmin is None:
                qmin = np.amin(self.__dataset[:,idx], axis=0)
            if qmax is None:
                qmax = np.zeros(len(idx))
                for i in range(len(idx)):
                    if i < self.__nphys:
                        qmax[i] = np.amax(self.__dataset[:,i])
                    else:
                        qmax[i] = np.amax(self.__dataset[:,i][
                            self.__dataset[:,i] < 99.0])
            griddims = []

            # Figure out how many dimensions the output has
            if hasattr(idx, '__len__'):
                nidx = len(idx)
            else:
                nidx = 1
            if nidx > 1:
                # Case for multiple indices
                griddims = []
                if hasattr(ngrid, '__len__'):
                    ngrid_tmp = np.array(ngrid, dtype=c_ulong)
                else:
                    ngrid_tmp = np.array([ngrid]*len(idx), dtype=c_ulong)
                for i in range(len(idx)):
                    griddims.append(qmin[i] + np.arange(ngrid_tmp[i]) * 
                                    float(qmax[i]-qmin[i])/(ngrid_tmp[i]-1))
                grid_out = np.squeeze(
                    np.array(np.meshgrid(*griddims,
                                         indexing='ij'), dtype=c_double))
                out_shape = grid_out[0, ...].shape
                qmin_tmp \
                    = np.copy(grid_out[(Ellipsis,)+(0,)*(grid_out.ndim-1)])
                qmax_tmp \
                    = np.copy(grid_out[(Ellipsis,)+(-1,)*(grid_out.ndim-1)])
            else:
                # Case for a single index
                ngrid_tmp = np.array([ngrid], dtype=c_ulong)
                qmin_tmp = np.array([qmin], dtype=np.double).reshape(1)
                qmax_tmp = np.array([qmax], dtype=np.double).reshape(1)
                grid_out = np.array(
                    qmin + np.arange(ngrid) *
                    float(qmax-qmin)/(ngrid-1),
                    dtype=c_double)
                out_shape = grid_out.shape

        # Return what we've built
        return nidx, qmin_tmp, qmax_tmp, grid_out, out_shape
    

    ##################################################################
    # Utility method to normalize marginal PDFs
    ##################################################################
    def __mpdf_normalize(self, nidx, grid_out, pdf):
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
        normfac = np.sum(pdf*np.abs(cellsize), axis = 
                         tuple(range(pdf.ndim-nidx, pdf.ndim)))

        # Normalize
        pdf = np.transpose(np.transpose(pdf)/normfac)

        # Return
        return pdf

    
    ##################################################################
    # Method to return log likelihood function at a specified set of
    # physical properties for a particular set of photometric
    # variables; this is set up for use by emcee, and is not intended
    # for use by humans
    ##################################################################
    def __logL(self, *args):
        x = np.zeros(self.__nphys+self.__nphot)
        x[:self.__nphys] = args[0]
        x[self.__nphys:] = args[1:-1]
        kd_tmp = args[-1]
        if not self.__diag_mode:
            return np.log(self.__clib.kd_pdf(kd_tmp, x, self.reltol,
                                             self.abstol))
        else:
            nodecheck = c_ulong(0)
            leafcheck = c_ulong(0)
            termcheck = c_ulong(0)
            return np.log(self.__clib.kd_pdf(kd_tmp, x, self.reltol,
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
        if (np.array(photprop).shape[-1] != self.__nphot) and \
           (self.__nphot > 1):
            raise ValueError("need " + str(self.__nphot) +
                             " photometric properties!")
        if photerr is not None:
            if (np.array(photerr).shape[-1] != self.__nphot) and \
               (self.__nphot > 1):
                raise ValueError("need " + str(self.__nphot) +
                                 " photometric errors!")
 
        # Reshape arrays if necessary
        if (np.array(photprop).shape[-1] != self.__nphot) and \
           (self.__nphot == 1):
            photprop = photprop.reshape(photprop.shape+(1,))
        if photerr is not None:
            if (np.array(photerr).shape[-1] != self.__nphot) and \
               (self.__nphot == 1):
                photerr = photerr.reshape(photerr.shape+(1,))

        # Prepare storage for output
        samples = []

        # Make dummy photometric errors if necessary
        if photerr is None:
            photerr = np.zeros(self.__nphot)

        # Loop over photometric errors
        kd_tmp = None
        for i in np.ndindex(*np.array(photerr).shape[:-1]):

            # Set bandwidth based on photometric error for this
            # iteration
            kd_tmp = self.__change_bw_err(photerr[i], kd_cur=kd_tmp)
 
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
                    kd_neighbors(kd_tmp, ph_tmp,
                                 dims.ctypes.data_as(POINTER(c_ulong)),
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
                                              args=[ph_tmp, kd_tmp])
                sampler.run_mcmc(pos, mc_steps)

                # Store the result
                samples.append(sampler.chain[:,mc_burn_in:,:].
                               reshape((-1,self.__nphys)))


        # Set bandwidth back to default if necessary
        if photerr is not None:
            self.__restore_bw_err(kd_tmp)

        # Reshape the samples
        samples = np.squeeze(
            np.array(samples).reshape(
                np.broadcast(np.array(photprop)[..., 0],
                             np.array(photerr)[..., 0]).shape +
                samples[0].shape))

        # Return
        return samples


    ##################################################################
    # Function to return the N best matches in the library to an input
    # set of photometric properties
    ##################################################################
    def bestmatch(self, phot, photerr=None, nmatch=1, 
                  bandwidth_units=False):
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

        Returns:
           matches : array, shape (..., nmatch, nphys + nfilter)
              best matches to the input photometry; shape in the
              leading dimensions will be the same as for phot, and if
              nmatch == 1 then that dimension will be omitted
           dist : array, shape (..., nmatch)
              distances between the matches and the input photometry
        """

        # Safety check
        if (np.array(phot).shape[-1] != self.__nphot) and \
           (self.__nphot > 1):
            raise ValueError("need " + str(self.__nphot) +
                             " photometric properties!")
        if photerr is not None:
            if (np.array(photerr).shape[-1] != self.__nphot) and \
               (self.__nphot > 1):
                raise ValueError("need " + str(self.__nphot) +
                                 " photometric errors!")
            if np.array(photerr).shape != np.array(phot).shape:
                raise ValueError("phot and photerr must have the same shape!")

        # Reshape arrays if necessary
        if (np.array(phot).shape[-1] != self.__nphot) and \
           (self.__nphot == 1):
            phot = phot.reshape(phot.shape+(1,))
        if photerr is not None:
            if (np.array(photerr).shape[-1] != self.__nphot) and \
               (self.__nphot == 1):
                photerr = photerr.reshape(photerr.shape+(1,))

        # Figure out how many distinct photometric values we've been
        # given, and how many sets of photometric errors
        nphot = np.array(phot).size//self.__nphot
        nphot_err = np.array(photerr).size//self.__nphot

        # Figure out what shape the output should have
        if nphot == 1 and nphot_err == 1:
            outshape = [self.__nphys + self.__nphot]
        elif photerr is None:
            outshape = list(phot.shape[:-1]) + \
                       [self.__nphys + self.__nphot]
        else:
            outshape = list(
                np.broadcast(np.array(phot)[..., 0],
                             np.array(photerr)[..., 0]).shape) + \
                [self.__nphys + self.__nphot]
        if nmatch > 1:
            outshape.insert(-1, nmatch)

        # Create array to hold output
        matches = np.zeros(outshape)
        if len(outshape) > 1:
            wgts = np.zeros(outshape[:-1])
            d2 = np.zeros(outshape[:-1])
        else:
            wgts = np.zeros([1])
            d2 = np.zeros([1])

        # Do we have errors?
        if photerr is None:

            # No, so just call the c neighbor-finding routine
            dims = np.arange(self.__nphys, self.__nphot+self.__nphys,
                             dtype=c_ulong)
            self.__clib.kd_neighbors_vec(
                self.__kd, np.ravel(phot), 
                dims.ctypes.data_as(POINTER(c_ulong)),
                self.__nphot, np.array(phot).size // self.__nphot,
                nmatch, bandwidth_units, np.ravel(matches),
                wgts.ctypes.data_as(POINTER(c_double)), 
                np.ravel(d2))

        elif nphot_err == 1:

            # Yes, we have errors, but only one set, so no need to
            # loop

            # Change the bandwidth to match the input errors
            kd_tmp = self.__change_bw_err(
                photerr, err_only = not bandwidth_units)

            # Now call c neighbor-finding routine
            dims = np.arange(self.__nphys, self.__nphot+self.__nphys,
                             dtype=c_ulong)
            self.__clib.kd_neighbors_vec(
                kd_tmp, np.ravel(phot), 
                dims.ctypes.data_as(POINTER(c_ulong)),
                self.__nphot, np.array(phot).size//self.__nphot,
                nmatch, True, np.ravel(matches),
                wgts.ctypes.data_as(POINTER(c_double)), 
                np.ravel(d2))

            # Restore bandwidth to previous value
            self.__restore_bw_err(kd_tmp)

        else:

            # We have multiple sets of errors, so loop
            ptr = 0
            kd_tmp = None
            for i in np.ndindex(*photerr.shape[:-1]):

                # Set bandwidth based on photometric error for this
                # iteration
                kd_tmp = self.__change_bw_err(
                    photerr[i], kd_cur = kd_tmp,
                    err_only = not bandwidth_units)

                # Call c neighbor-finding routine
                offset1 = ptr*nmatch
                offset2 = offset1*(self.__nphot+self.__nphys)
                dims = np.arange(self.__nphys, self.__nphot+self.__nphys,
                                 dtype=c_ulong)                
                self.__clib.kd_neighbors_vec(
                    kd_tmp, np.ravel(phot[i]), 
                    dims.ctypes.data_as(POINTER(c_ulong)),
                    self.__nphot, np.array(phot[i]).size//self.__nphot,
                    nmatch, True, np.ravel(matches)[offset2:],
                    wgts.ctypes.data_as(POINTER(c_double)), 
                    np.ravel(d2)[offset1:])
                ptr = ptr+1

            # Restore bandwidth to previous value
            self.__restore_bw_err(kd_tmp)

        # Return
        return matches, np.sqrt(d2)

    ##################################################################
    # Function to return the N best matches in the library to an input
    # set of physical properties
    ##################################################################
    def bestmatch_phys(self, phys, nmatch=1, bandwidth_units=False):
        """
        Searches through the simulation library and returns the closest
        matches to an input set of photometry.

        Parameters:
           phot : arraylike, shape (nphys) or (..., nphys)
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

        Returns:
           matches : array, shape (..., nmatch, nphys + nfilter)
              best matches to the input properties; shape in the
              leading dimensions will be the same as for phot, and if
              nmatch == 1 then that dimension will be omitted
           dist : array, shape (..., nmatch)
              distances between the matches and the input physical
              properties
        """

        # Safety check
        if (np.array(phys).shape[-1] != self.__nphys) and \
           (self.__nphot > 1):
            raise ValueError("need " + str(self.__nphys) +
                             " physical properties!")


        # Figure out how many points we were given
        phys_tmp = np.array(phys)
        npt = phys_tmp.size // self.__nphys

        # Prepare the output arrays
        outshape = list(phys_tmp.shape[:-1]) + \
                   [self.__nphys + self.__nphot]
        if nmatch > 1:
            outshape.insert(-1, nmatch)

        # Create array to hold output
        matches = np.zeros(outshape)
        if len(outshape) > 1:
            wgts = np.zeros(outshape[:-1])
            d2 = np.zeros(outshape[:-1])
        else:
            wgts = np.zeros([1])
            d2 = np.zeros([1])

        # Call c routine
        dims = np.arange(self.__nphys, dtype=c_ulong)
        self.__clib.kd_neighbors_vec(
            self.__kd, np.ravel(phys_tmp), 
            dims.ctypes.data_as(POINTER(c_ulong)),
            self.__nphys, npt, nmatch, bandwidth_units, 
            np.ravel(matches),
            wgts.ctypes.data_as(POINTER(c_double)), 
            d2)

        # Return
        return matches, np.sqrt(d2)


    ##################################################################
    # Functions to return an approximate kernel density representation
    # of a given point; that is, for an input physical point or
    # photometric point (there are two versions below), this routine
    # returns a set of points and weights that can be used to compute
    # an approximation to the PDF for it anywhere in space, much
    # faster than could be done using the full data set
    ##################################################################
    def make_approx_phot(self, phys, squeeze=True, filter_ignore=None):
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

        Returns:
           x : array, shape (M, nphot), or a list of such arrays
              an array containing the list of points to be used for
              the approximation, where nphot is the number of
              photometric filters being returned
           wgts : array, shape (M), or a list of such arrays
              an array containing the weights of the points

        Notes:
           if the requested relative tolerance cannot be reached for
           numerical reasons (usually because the input point is too
           far from the library to allow accurate computation), x and
           wgts will be return as None, and a warning will be issued
        """

        # If given only one set of physical variables, make a
        # temporary list we can loop over
        if hasattr(phys[0], '__iter__'):
            phystmp = phys
        else:
            phystmp = [phys]

        # Specify dimensions to return, and set bandwidth for them
        if filter_ignore is None:
            dim_return_ptr = None
            ndim_return = self.__nphot
            bw = self.__bandwidth[self.__nphys:]
        else:
            dim_return = np.where(np.logical_not(
                np.array(filter_ignore)))[0] + self.__nphys
            dim_return_ptr = dim_return.ctypes.data_as(POINTER(c_ulong))
            ndim_return = len(dim_return)
            bw = self.__bandwidth[dim_return]

        # Loop over input physical variables
        x = []
        wgts = []
        for ph in phystmp:

            # Safety check
            if len(ph) != self.__nphys:
                raise ValueError("need " + str(self.__nphys) + 
                                 " physical properties!")

            # Call c library routine to make the representation
            xout = POINTER(c_double)()
            wgtsout = POINTER(c_double)()
            npts = self.__clib.kd_rep(self.__kd, np.array(ph), 
                                      np.arange(self.__nphys, dtype=c_ulong),
                                      self.__nphys, self.reltol, 
                                      dim_return_ptr, ndim_return,
                                      ctypes.byref(xout), 
                                      ctypes.byref(wgtsout))

            # Check if the c routine encountered an error due to
            # insufficient precision; if so, issue warning and return
            # None in place of xpt and wgts
            if npts == 0:
                warn("bp.make_approx_phot: requested precision cannot be reached")
                return None

            # Call c library routine to squeeze the representation
            if squeeze:
                npts = self.__clib.squeeze_rep(
                    npts, ndim_return, bw, self.reltol,
                    ctypes.byref(xout), ctypes.byref(wgtsout))

            # Convert the returned values into numpy arrays; copy the
            # data so that we can free the c buffers
            xsave = np.copy(npct.as_array(xout, shape=(npts, ndim_return)))
            wgtssave = np.copy(npct.as_array(wgtsout, shape=(npts,)))

            # Free the c buffers
            self.__clib.free_kd_rep(ctypes.byref(xout), 
                                    ctypes.byref(wgtsout))

            # Append to the lists we'll be returning
            x.append(xsave)
            wgts.append(wgtssave)

        # If we were given a single object, return something in the
        # same shape
        if not hasattr(phys[0], '__iter__'):
            x = x[0]
            wgts = wgts[0]

        # Return
        return x, wgts

    def make_approx_phys(self, phot, photerr=None, squeeze=True,
                         phys_ignore=None, tol=None):
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

        # If given only one set of photometric variables, make a
        # temporary list we can loop over; same for photometric errors
        if hasattr(phot[0], '__iter__'):
            phottmp = phot
        else:
            phottmp = [phot]
        if photerr is not None:
            if hasattr(photerr[0], '__iter__'):
                photerrtmp = photerr
            else:
                photerrtmp = [photerr]
        else:
            photerrtmp = [None]

        # Specify dimensions to return, and set bandwidth for them
        if phys_ignore is None:
            dim_return_ptr = None
            ndim_return = self.__nphys
            bw = self.__bandwidth[:self.__nphys]
        else:
            dim_return = np.where(np.logical_not(
                np.array(phys_ignore)))[0]
            dim_return_ptr = dim_return.ctypes.data_as(POINTER(c_ulong))
            ndim_return = len(dim_return)
            bw = self.__bandwidth[dim_return]

        # Set tolerance
        if tol is None:
            tol = self.reltol

        # Loop over input photometric variables
        x = []
        wgts = []
        for ph in phottmp:

            # Safety check
            if len(ph) != self.__nphot:
                raise ValueError("need " + str(self.__nphot) + 
                                 " photometric properties!")

            # Loop over photometric errors
            kd_tmp = None
            for pherr in photerrtmp:

                # Change the bandwidth if necessary
                if pherr is not None:
                    if len(pherr) != self.__nphot:
                        raise ValueError("need " + str(self.__nphot) + 
                                         " photometric errors!")
                    kd_tmp = self.__change_bw_err(pherr, kd_cur=kd_tmp)

                # Call c library routine to make the representation
                xout = POINTER(c_double)()
                wgtsout = POINTER(c_double)()
                npts = self.__clib.kd_rep(
                    kd_tmp, np.array(ph), 
                    np.arange(self.__nphot, dtype=c_ulong)+self.__nphys,
                    self.__nphot, tol, 
                    dim_return_ptr, ndim_return,
                    ctypes.byref(xout), 
                    ctypes.byref(wgtsout))

                # Check if the c routine encountered an error due to
                # insufficient precision; if so, issue warning and
                # return None in place of xpt and wgts
                if npts == 0:
                    warn("bp.make_approx_phys: requested precision cannot be reached")
                    if kd_tmp is not None:
                        self.__restore_bw_err(kd_tmp)
                    return None

                # Call c library routine to squeeze the representation
                if squeeze:
                    npts = self.__clib.squeeze_rep(
                        npts, ndim_return, bw, self.reltol,
                        ctypes.byref(xout), ctypes.byref(wgtsout))

                # Convert the returned values into numpy arrays; copy the
                # data so that we can free the c buffers
                xsave = np.copy(npct.as_array(xout, 
                                              shape=(npts, ndim_return)))
                wgtssave = np.copy(npct.as_array(wgtsout, shape=(npts,)))

                # Free the c buffers
                self.__clib.free_kd_rep(ctypes.byref(xout), 
                                        ctypes.byref(wgtsout))

                # Append to the lists we'll be returning
                x.append(xsave)
                wgts.append(wgtssave)

        # Reset the bandwidth
        if photerr is not None:
            self.__restore_bw_err(kd_tmp)

        # If we were given a single object, return something in the
        # same shape
        if not hasattr(phot[0], '__iter__'):
            if photerr is None:
                x = x[0]
                wgts = wgts[0]
            elif not hasattr(photerr[0], '__iter__'):
                x = x[0]
                wgts = wgts[0]

        # Return
        return x, wgts
            
    ##################################################################
    # Method to squeeze a representation that has already been created
    ##################################################################
    def squeeze_rep(self, x, wgts, dims=None):
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

        Returns:
           Nothing
        """

        # Get the bandwidth
        if dims is None:
            bw = self.__bandwidth[:x.shape[1]]
        else:
            bw = self.__bandwidth[dims]

        # Prepare data to pass to the c routine
        npt = x.shape[0]
        ndim = x.shape[1]
        xtmp = POINTER(c_double)(x.ctypes.data_as(POINTER(c_double)))
        wgtstmp = POINTER(wgts.ctypes.data_as(POINTER(c_double)))

        # Call the c squeeze routine
        npt_out = self.__clib.squeeze_rep(
            npt, ndim, bw, self.reltol,
            ctypes.byref(xtmp), ctypes.byref(wgtstmp))

        # Discard the old metadata for x and wgts, and rebuild them
        # using the new metadata
        x = npct.as_array(xtmp, shape=(npt_out, ndim))
        wgts = npct.as_array(wgtstmp, shape=(ndim,))


    ##################################################################
    # Routines to use the approximate representations returned by
    # make_approx_phys and make_approx_phot to compute marginal
    # posterior probabilities
    ##################################################################

    def mpdf_approx(self, x, wgts, dims='phys', dims_return=None,
                    ngrid=64, qmin='all', qmax='all', grid=None,
                    norm=True):
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
              specified in dms; if not, then dimreturn must be a
              subset of dim, and a marginal PDF in certain dimensions
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

        # Set up input and output dimensions
        if dims == 'phys':
            dim_in = np.arange(self.__nphys)
        elif dims == 'phot':
            dim_in = np.arange(self.__nphot) + self.__nphys
        else:
            dim_in = dims
        if dims_return is None:
            dim_out = dim_in
        elif not hasattr(dims_return, '__iter__'):
            dim_out = [dims_return]
            if not dim_out in dim_in:
                raise ValueError("dims_return must be a subset of dims!")
        elif set(dims_return) <= set(dim_in):
            dim_out = dims_return
        else:
            raise ValueError("dims_return must be a subset of dims!")
        nidx = len(dim_out)

        # Set up the output grid
        if grid is not None:
            grid_out = np.meshgrid(grid)
            qmin \
                = np.array(
                    np.copy(grid_out[(Ellipsis,)+(0,)*(grid_out.ndim-1)]),
                    dtype=np.double)
            qmax \
                = np.array(
                    np.copy(grid_out[(Ellipsis,)+(-1,)*(grid_out.ndim-1)]),
                    dtype=np.double)
            ngrid_tmp = np.squeeze(grid_out).shape
            for i in range(len(dim_out)):
                griddims.append(np.linspace(qmin[i], qmax[i], 
                                            ngrid_tmp[i]))
        else:
            if qmin == 'zoom' or qmin == 'all' or \
               qmax == 'zoom' or qmax == 'all':
                if qmin is 'zoom' or qmin == 'all':
                    get_qmin = True
                    qmin = []
                else:
                    get_qmin = False
                if qmax == 'zoom' or qmax == 'all':
                    qmax = []
                    get_qmax = True
                else:
                    get_qmax = False
                for d in dim_out:
                    xtmp = np.copy(x[:,d])
                    wgttmp = np.copy(wgts)
                    idx = np.argsort(xtmp)
                    xtmp = xtmp[idx]
                    wgttmp = wgttmp[idx]
                    wgtsum = np.cumsum(wgttmp)
                    if get_qmin:
                        if qmin == 'all':
                            qmin.append(xtmp[0])
                        else:
                            qmin.append(xtmp[np.argmax(wgtsum > 0.001)])
                    if get_qmax:
                        if wgtsum[-1] > 0.999 and qmax is not 'all':
                            qmax.append(xtmp[np.argmax(wgtsum >
                                                       0.999)])
                        else:
                            qmax.append(xtmp[-1])
                if get_qmin:
                    qmin = np.array(qmin)
                if get_qmax:
                    qmax = np.array(qmax)
            griddims = []
            if nidx > 1:
                # Case for multiple indices
                griddims = []
                if hasattr(ngrid, '__len__'):
                    ngrid_tmp = np.array(ngrid)
                else:
                    ngrid_tmp = np.array([ngrid]*len(dim_out))
                for i in range(len(dim_out)):
                    griddims.append(np.linspace(qmin[i], qmax[i], 
                                                ngrid_tmp[i]))
                grid_out = np.squeeze(np.array(np.meshgrid(*griddims,
                                                           indexing='ij')))
            else:
                # Case for a single index
                grid_out = qmin + \
                           np.arange(ngrid) * \
                           float(qmax-qmin)/(ngrid-1)
                griddims = [grid_out]

        # Compute result; in the multi-dimensional case this involves
        # some tricky array indexing
        if nidx == 1:
            pdf = np.einsum('i,ij', wgts, 
                           np.exp(-np.subtract.outer(
                               x[:,dim_out[0]], grid_out)**2 / 
                                  (2*self.__bandwidth[dim_out[0]]**2)))
                         
        else:
            compsum = np.exp(
                -np.subtract.outer(x[:,dim_out[0]], griddims[0])**2 /
                (2.0*self.__bandwidth[dim_out[0]]**2))
            for d, grd in zip(dim_out[1:], griddims[1:]):
                comp = np.exp(-np.subtract.outer(x[:,d], grd)**2 / \
                              (2.0*self.__bandwidth[d]**2))
                compsum = np.einsum('...j,...k->...jk', compsum, comp)
            pdf = np.einsum('i,i...', wgts, compsum)

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
            normfac = np.sum(pdf*cellsize)

            # Normalize
            pdf = pdf/normfac

        # Return
        return grid_out, pdf


    ##################################################################
    # Routines to draw from the kernel density PDF, constraining 
    # dimensions or not
    ##################################################################
    def draw_sample(self, photerr=None, nsample=1):
        """
        Returns a randomly-drawn sample of physical and photometric
        properties from the kernel density function

        Parameters:
           nsample : int
              number of random samples to draw for each set of
              physical properties; must be positive
           photerr : arraylike, shape (nphot)
              photometric errors to apply to the output photometry;
              these are added in quadrature with the kernel density
              estimation bandwidth

        Returns:
           samples : array, shape (nsample, nphys+nphot)
              random sample drawn from the kernel density object; for
              the final dimension in the output, the first nphys
              elements are the physical quantities, the next nphot are
              the photometric quantities
        """

        # Safety check
        if (nsample < 1):
            raise ValueError("need to draw at least 1 sample")
        
        # Allocate a random number generator from gsl if we do not
        # already have one
        if not hasattr(self, '__rng'):
            self.__rng = self.__clib.rng_init(0)

        # Allocate array to hold outputs
        samples = np.zeros((nsample, self.__nphys+self.__nphot),
                           dtype=c_double)

        # Apply photometric errors if requested
        if photerr is None:
            kd_tmp = self.__kd
        else:
            kd_tmp = self.__change_bw_err(photerr)

        # Call library function
        self.__clib.kd_pdf_draw(kd_tmp, None, None, 0,
                                nsample, 1, self.__rng,
                                np.ravel(samples))

        # Restore kernel density bandwidth
        if photerr is not None:
            self.__restore_bw_err(kd_tmp)

        # Return samples
        return samples
            
    
    def draw_phot(self, physprop, physidx=None, photerr=None,
                  nsample=1):
        """
        Returns a randomly-drawn sample of photometric properties for
        one or more input sets of physical properties.

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

        Returns:
           samples : arraylike
              sample of photometric properties; the shape of the
              output is (..., nsample, nphot), where nphot is the
              number of photometric properties, and the leading
              dimension(s) match the leading dimension(s) of physprop
        """

        # Safety check and sanitization on inputs
        if (nsample < 1):
            raise ValueError("need to draw at least 1 sample")
        if physidx is None:
            physidx_ = np.arange(self.__nphys, dtype=c_ulong)
        else:
            physidx_ = np.flatten(np.array(physidx_))
            if (len(physidx_) > self.__nphys) or \
               (len(physidx_) != len(np.unique(physidx))) or \
               (np.amin(physidx_) < 0) or \
               (np.amax(physidx_) >= self.__nphys):
                raise ValueError("physidx must be an array of "
                                 "unique indices from 0 to {:d}".
                                 format(self.__nphys-1))
            
        physprop_ = np.array(physprop)
        if physprop_.shape[-1] != physidx_.size:
            raise ValueError("physidx must have as many elements as "
                             "trailing dimension of physprop")

        # Set up array to hold outputs
        if physprop_.ndim == 1:
            physprop_ = physprop_.reshape((1,)+physprop_.shape)
        sample_dim = physprop_.shape[:-1]
        samples = np.zeros(sample_dim + (nsample, self.__nphot))

        # Apply photometric errors if requested
        if photerr is None:
            kd_tmp = self.__kd
        else:
            kd_tmp = self.__change_bw_err(photerr)

        # Allocate a random number generator from gsl if we do not
        # already have one
        if not hasattr(self, '__rng'):
            self.__rng = self.__clib.rng_init(0)

        # Loop over leading dimensions of input physical properties
        for i in np.ndindex(*physprop_.shape[:-1]):

            # Call library function
            self.__clib.kd_pdf_draw(
                kd_tmp,
                physprop_[i].ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                physidx_.ctypes.data_as(ctypes.POINTER(c_ulong)),
                physidx_.size,
                nsample, 0,
                self.__rng,
                np.ravel(samples[i]))

        # Restore kernel density bandwidth
        if photerr is not None:
            self.__restore_bw_err(kd_tmp)

        # Reshape output if needed
        if physprop_.ndim == 1:
            samples = samples[0]

        # Return samples
        return samples
