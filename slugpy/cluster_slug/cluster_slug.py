"""
This defines a class that can be used to estimate the PDF of star
cluster properties (mass, age, extinction) from a set of input
photometry in various bands.
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
import numpy.ctypeslib as npct
import warnings
from copy import deepcopy
import emcee
from ..read_cluster_prop import read_cluster_prop
from ..read_cluster_phot import read_cluster_phot

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

class cluster_slug(object):
    """
    A class that can be used to estimate the PDF of star cluster
    properties (mass, age, extinction) from a set of input photometry
    in various bands.
    """

    ##################################################################
    # Initializer method
    ##################################################################
    def __init__(self, libname=None, bandwidth=0.1):
        """
        Initialize a cluster_slug object.

        Parameters
           libname : string
              name of the SLUG model to load; if left as None, the default
              is $SLUG_DIR/cluster_slug/CLUSTER_SLUG
           bandwidth : float
             bandwidth of the kernel to use in density estimates

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

        # Only keep cases where the luminosity is non-zero; zero
        # luminosity cases correspond to clusters where the most
        # massive star produced was below the minimum mass included in
        # our evolutionary tracks
        idx = np.where(np.amin(phot.phot, axis=1) > 0)[0]

        # Store filters
        self.__dataset_filters = phot.filter_names

        # Build dataset array; 1st 3 dimensions are log mass, log age,
        # and A_V remaining ones are photometric values
        self.__dataset = np.zeros((len(idx),
                                   3+len(phot.filter_names)))
        self.__dataset[:,0] = np.log10(prop.actual_mass[idx])
        self.__dataset[:,1] = np.log10(prop.time[idx] - prop.form_time[idx])
        self.__dataset[:,2] = prop.A_V[idx]
        self.__dataset[:,3:] = phot.phot_ex[idx,:]

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
            self.__dataset[:,3+i] = phot.phot[idx,i]

        # Take log of photometric values if they are recorded in a
        # linear system
        for i, f in enumerate(phot.filter_units):
            if 'mag' not in f:
                self.__dataset[:,3+i] = np.log10(self.__dataset[:,3+i])

        # Record the bandwidth
        self.__bandwidth = bandwidth

        # Set up the interface to the cluster_slug c library; note
        # that certain arrays have to be decleared as
        # POINTER(c_double) instead of array_1d_double, because they
        # can be passed either an array or a NULL
        self.__clib = npct.load_library("cluster_slug", 
                                        osp.realpath(__file__))
        self.__clib.build_kd.restype = c_void_p
        self.__clib.build_kd.argtypes \
            = [ array_1d_double,   # x
                c_uint,            # ndim
                c_uint,            # npt
                ctypes.
                POINTER(c_double), # wgt
                c_uint,            # leafsize
                c_double,          # bandwidth
                c_int ]            # ktype
        self.__clib.free_kd.restype = None
        self.__clib.free_kd.argtypes = [ c_void_p ]
        self.__clib.kd_pdf.restype = c_double
        self.__clib.kd_pdf.argtypes \
            = [ array_1d_double,   # x
                c_void_p ]         # kd
        self.__clib.kd_pdf_vec.restype = None
        self.__clib.kd_pdf_vec.argtypes \
            = [ array_1d_double,   # x
                c_uint,            # npt
                c_void_p,          # kd
                array_1d_double ]  # pdf
        self.__clib.kd_pdf_int.restype = c_double
        self.__clib.kd_pdf_int.argtypes \
            = [ array_1d_double,   # q
                array_1d_uint,     # qdim
                c_uint,            # nqdim
                c_void_p ]         # kd
        self.__clib.reweight_kd.restype = None
        self.__clib.reweight_kd.argtypes \
            = [ POINTER(c_double), # wgt
                c_void_p ]         # kd
        self.__clib.find_neighbors.restype = None
        self.__clib.find_neighbors.argtypes \
            = [ array_1d_double,   # xpt
                POINTER(c_uint),  # dims
                c_uint,            # ndim
                c_uint,            # nneighbor
                c_void_p,          # kd
                array_1d_double,   # pos
                POINTER(c_double), # dptr
                array_1d_double ]  # dist2

        # Initialize the stored filter list, kernel_density objects,
        # and priors
        self.__stored_filters = None
        self.__dataset_mass_pdf = None
        self.__dataset_age_pdf = None
        self.__dataset_AV_pdf = None
        self.__dataset_wgts = None
        self.__kd = None
        self.__prior_mass = None
        self.__prior_age = None
        self.__prior_AV = None


    ##################################################################
    # De-allocation method
    ##################################################################
    def __del__(self):
        if self.__kd is not None:
            self.__clib.free_kd(
                ctypes.byref(c_void_p(self.__kd)))


    ##################################################################
    # Return a copy of the list of available filters
    ##################################################################
    def filters(self):
        return deepcopy(self.__dataset_filters)


    ##################################################################
    # Method to set up kernel density estimation
    ##################################################################
    def setup_kd(self, filter_names, prior_mass=None,
              prior_age=None, prior_AV=None, prior_mass_args=None, 
              prior_age_args=None, prior_AV_args=None):
        """
        This function sets up kernel density estimation

        Parameters:
           filter_names : iterable of stringlike
              names of the photometric filters being used
           prior_mass : callable
              a callable specifying the prior probability distribution
              on mass; must accept an array of mass values in Msun,
              and return an array of the same shape giving the prior
              probability at each mass. If left as None, the prior
              used is left as whatever distribution was present in the
              input data.
           prior_age : callable
              a callable specifying the prior probability distribution
              on age; must accept an array of age values in yr, 
              and return an array of the same shape giving the prior
              probability at each age. If left as None, the prior used
              is left as whatever distribution was present in the
              input data.
           prior_AV : callable
              a callable specifying the prior probability distribution
              on A_V; must accept an array of age values in mag, 
              and return an array of the same shape giving the prior
              probability at each A_V. If left as None, the prior used
              is left as whatever distribution was present in the
              input data.
           prior_mass_args : dict
              a dict of additional arguments to be passed to the
              prior_mass function
           prior_age_args : dict
              a dict of additional arguments to be passed to the
              prior_mass function
           prior_AV_args : dict
              a dict of additional arguments to be passed to the
              prior_mass function

        Returns
           Nothing
        """

        # If we are going to use a prior on mass, and we haven't yet
        # figured out the mass PDF of the input data set, do so now
        if prior_mass is not None and self.__dataset_mass_pdf is None:

            # Store the list of input masses
            mass_pdf_data = np.copy(self.__dataset[:,0])

            # Build a kernel density estimator for the PDF of masses
            # in the input data set
            mass_pdf_kd \
                = self.__clib.build_kd(mass_pdf_data, 1,
                                       mass_pdf_data, None, 32, 
                                       self.__bandwidth, 0)

            # Use the kernel density estimator to compute the PDF of
            # input masses; store this for future use
            mass_list = np.copy(self.__dataset[:,0])
            self.__dataset_mass_pdf = np.zeros(mass_list.shape[0])
            self.__clib.kd_pdf_vec(mass_list, mass_list.shape[0],
                                   mass_pdf_kd,
                                   self.__dataset_mass_pdf)

            # Free the kernel density object used to make this
            # estimate
            self.__clib.free_kd(ctypes.byref(c_void_p(mass_pdf_kd)))

        # Repeat the procedure if we're going to need the PDF of age
        # or AV in the input data set
        if prior_age is not None and self.__dataset_age_pdf is None:
            age_pdf_data = np.copy(self.__dataset[:,1])
            age_pdf_kd \
                = self.__clib.build_kd(age_pdf_data, 1,
                                       age_pdf_data, None, 32, 
                                       self.__bandwidth, 0)
            age_list = np.copy(self.__dataset[:,1])
            self.__dataset_age_pdf = np.zeros(age_list.shape[0])
            self.__clib.kd_pdf_vec(age_list, age_list.shape[0],
                                   age_pdf_kd,
                                   self.__dataset_age_pdf)
            self.__clib.free_kd(ctypes.byref(c_void_p(age_pdf_kd)))
        if prior_AV is not None and self.__dataset_AV_pdf is None:
            AV_pdf_data = np.copy(self.__dataset[:,2])
            AV_pdf_kd \
                = self.__clib.build_kd(AV_pdf_data, 1,
                                       AV_pdf_data, None, 32, 
                                       self.__bandwidth, 0)
            AV_list = np.copy(self.__dataset[:,2])
            self.__dataset_AV_pdf = np.zeros(AV_list.shape[0])
            self.__clib.kd_pdf_vec(AV_list, AV_list.shape[0],
                                   AV_pdf_kd,
                                   self.__dataset_AV_pdf)
            self.__clib.free_kd(ctypes.byref(c_void_p(AV_pdf_kd)))

        # If using a different prior than the last time we were
        # called, compute the new weights that result from that
        # prior. Note that, since the user function returns p(mass),
        # but we're working with log(mass), we have to multiply the
        # PDF returned by the user's function by mass to get it into
        # logarithmic units. The same goes for age, but not for A_V
        # which is already logarithmic.
        reweight_pdf = False
        if prior_mass != self.__prior_mass:
            reweight_pdf = True
            self.__prior_mass = prior_mass
            if prior_mass is not None:
                self.__prior_mass_wgt \
                    = self.__dataset[:,0] * \
                    prior_mass(self.__dataset[:,0]) \
                    / self.__dataset_mass_pdf
            else:
                self.__prior_mass_wgt = None
        if prior_age != self.__prior_age:
            reweight_pdf = True
            self.__prior_age = prior_age
            if prior_age is not None:
                self.__prior_age_wgt \
                    = self.__dataset[:,1] * \
                    prior_age(self.__dataset[:,1]) \
                    / self.__dataset_age_pdf
            else:
                self.__prior_age_wgt = None
        if prior_AV != self.__prior_AV:
            reweight_pdf = True
            self.__prior_AV = prior_AV
            if prior_AV is not None:
                self.__prior_AV_wgt \
                    = prior_AV(self.__dataset[:,2]) \
                    / self.__dataset_AV_pdf
            else:
                self.__prior_AV_wgt = None

        # Compute new weight for overall PDF
        if reweight_pdf:
            if self.__prior_mass_wgt is not None or \
               self.__prior_age_wgt is not None or \
               self.__prior_AV_wgt is not None:
                self.__dataset_wgts = np.ones(self.__dataset.shape[0])
                if self.__prior_mass_wgt is not None:
                    self.__dataset_wgts \
                        = self.__dataset_wgts * self.__prior_mass_wgt
                if self.__prior_age_wgt is not None:
                    self.__dataset_wgts \
                        = self.__dataset_wgts * self.__prior_age_wgt
                if self.__prior_AV_wgt is not None:
                    self.__dataset_wgts \
                        = self.__dataset_wgts * self.__prior_AV_wgt
            else:
                self.__dataset_wgts = None

        # If we have not yet built the PDF of the full data set for
        # this list of filters, do so now
        if self.__stored_filters != filter_names:
            self.__stored_filters = filter_names
            self.__dataset_tmp = np.zeros((self.__dataset.shape[0],
                                           3+len(filter_names)))
            self.__dataset_tmp[:,:3] = self.__dataset[:,:3]
            for i, f in enumerate(filter_names):
                if f not in self.__dataset_filters:
                    raise ValueError("unknown filter "+str(f))
                self.__dataset_tmp[:,3+i] \
                    = self.__dataset[:, 
                                     3+self.__dataset_filters.index(f)]
            self.__kd \
                = self.__clib.build_kd(np.ravel(self.__dataset_tmp),
                                       3 + len(filter_names),
                                       self.__dataset_tmp.shape[0],
                                       self.__dataset_wgts,
                                       32, self.__bandwidth, 0)
            reweight_pdf = False

        # If we have a built PDF, but we need to reweight it, do so
        # now
        if reweight_pdf:
            self.__clib.reweight_kd(self.__dataset_wgts, self.__kd)


    ##################################################################
    # Method to return log likelihood function at a specified mass,
    # age, AV combination for a particular set of photometric
    # variables
    ##################################################################
    def logL(self, *args):
        x = np.zeros(3+len(self.__stored_filters))
        x[:3] = args[0]
        x[3:] = args[1:]
        return np.log(self.__clib.kd_pdf(x, self.__kd))


    ##################################################################
    # Function to compute an MCMC sampler for a particular set of
    # photometric values
    ##################################################################
    def pdf(self, phot, filter_names, prior_mass=None,
            prior_age=None, prior_AV=None, prior_mass_args=None, 
            prior_age_args=None, prior_AV_args=None, mc_walkers=100,
            mc_steps=500, mc_burn_in=50):
        """
        This function returns the posterior probability distribution
        for 

        Parameters:
           phot : array, shape (nfilter) or (ncluster, nfilter)
              array giving photometric values for 1 or more clusters
           filter_names : iterable of stringlike
              names of the photometric filters being used
           prior_mass : callable
              a callable specifying the prior probability distribution
              on mass; must accept an array of mass values in Msun,
              and return an array of the same shape giving the prior
              probability at each mass. If left as None, the prior
              used is left as whatever distribution was present in the
              input data.
           prior_age : callable
              a callable specifying the prior probability distribution
              on age; must accept an array of age values in yr, 
              and return an array of the same shape giving the prior
              probability at each age. If left as None, the prior used
              is left as whatever distribution was present in the
              input data.
           prior_AV : callable
              a callable specifying the prior probability distribution
              on A_V; must accept an array of age values in mag, 
              and return an array of the same shape giving the prior
              probability at each A_V. If left as None, the prior used
              is left as whatever distribution was present in the
              input data.
           prior_mass_args : dict
              a dict of additional arguments to be passed to the
              prior_mass function
           prior_age_args : dict
              a dict of additional arguments to be passed to the
              prior_mass function
           prior_AV_args : dict
              a dict of additional arguments to be passed to the
              prior_mass function
           mc_walkers : int
              number of walkers to use in the MCMC
           mc_steps : int
              number of steps in the MCMC
           mc_burn_in : int
              number of steps to consider "burn-in" and discard

        Returns
           Nothing
        """

        # Set up kernel density for this filter set and set of priors
        self.setup_kd(filter_names, 
                      prior_mass=prior_mass,
                      prior_age=prior_age,
                      prior_AV=prior_AV, 
                      prior_mass_args=prior_mass_args, 
                      prior_age_args=prior_age_args,
                      prior_AV_args=prior_AV_args)

        # If only 1 cluster, reshape array to be 1 x n_filters
        if phot.ndim == 1:
            phot_tmp = phot.reshape(1, phot.shape[0])
        else:
            phot_tmp = phot

        # Prepare storage for output
        samples = []

        # Loop over clusters
        for i in range(phot_tmp.shape[0]):

            # Grab photometric values for this cluster
            ph_tmp = phot_tmp[i,:]

            # Search the data set for the sample cluster closest to
            # the observed luminosities; this will be the starting
            # point for our search
            nfilter = len(filter_names)
            dims = np.arange(3, 3+nfilter, dtype=np.uint32)
            nearpt = np.zeros(3+nfilter)
            wgt = np.zeros(1)
            dist2 = np.zeros(1)
            self.__clib. \
                find_neighbors(ph_tmp, 
                               dims.ctypes.data_as(
                                   ctypes.POINTER(ctypes.c_uint)),
                               nfilter,
                               1, self.__kd, nearpt, 
                               wgt.ctypes.data_as(
                                   ctypes.POINTER(ctypes.c_double)),
                               dist2)

            # Generate a set of starting points by scattering walkers
            # around the starting position
            pos = [nearpt[:3] + self.__bandwidth*np.random.randn(3) 
                   for i in range(mc_walkers)]

            # Run the MCMC
            sampler=emcee.EnsembleSampler(mc_walkers, 3, 
                                          self.logL, 
                                          args=ph_tmp)
            sampler.run_mcmc(pos, mc_steps)

            # Store the result
            samples.append(sampler.chain[:,mc_burn_in:,:].
                           reshape((-1,3)))

        # Remove extraneous dimension if given just one set of
        # photometry
        if phot.ndim == 1:
            samples = samples[0]

        # Return
        return samples
