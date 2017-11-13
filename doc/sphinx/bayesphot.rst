.. highlight:: rest

.. _sec-bayesphot:

bayesphot: Bayesian Inference for Stochastic Stellar Populations
================================================================

What Does bayesphot Do?
-----------------------

Bayesphot is a package for performing Bayesian inference for the physical properties of a stellar system using its measured photometric properties, in a case where the photometric properties vary non-deterministically with the physical properties. Formally, bayesphot answers the following question: consider a stellar system characterized by a vector of :math:`\mathbf{x} = (x_1, x_2, \ldots x_N)` physical properties. We have a physical model that lets us sample the expected photometric properties as a function of physical properties, i.e., that for some sample of :math:`K` systems with physical properties :math:`\mathbf{x}_k` we are able to compute the corresponding photometric properties :math:`\mathbf{y}_k = \mathbf{y} = (y_1, y_2, \ldots y_M)_k`. Now suppose that we observe such a system, and we observe it to have photometric properties :math:`\mathbf{y}_{\mathrm{obs}}`, with some set of photometric errors :math:`\mathbf{\sigma}_{\mathbf{y}} = (\sigma_{y_1}, \sigma_{y_2}, \ldots \sigma_{y_M})`, which are assumed to be Gaussian-distributed. What should we infer about the posterior probability distribution of the physical properties, i.e., given a set of prior probabilities :math:`p(\mathbf{x})`, plus our measurements, what is :math:`p(\mathbf{x} \mid \mathbf{y}_{\mathrm{obs}}, \mathbf{\sigma}_{\mathrm{y}})`?

The kernel density estimation algorithm that bayesphot uses to answer this question is described and derived in the slug methods paper. Bayesphot is implemented in two parts: a shared object library that is implemented in c, and that is compiled at the same time that slug is built, and a python wrapper class called ``bp`` that is included in the slugpy.bayesphot module. The following sections describe how to use ``bp`` objects to generate posterior PDFs.

Creating ``bp`` Objects
-----------------------

The ``bp`` class can be imported via::

  from slugpy.bayesphot import *

or::

  from slugpy.bayesphot import bp

Once imported, a ``bp`` object can be instantiated. The call signature for the ``bp`` constructor class is::

  def __init__(self, dataset, nphys, filters=None, bandwidth='auto',
               ktype='gaussian', priors=None, sample_density=None,
               reltol=1.0e-3, abstol=1.0e-10, leafsize=16):

A full description of all options is included in the :ref:`ssec-slugpy-bayesphot`, but the essential features are summarized here.

The argument ``dataset`` is an array of shape (N, M) that contains the library of N models that represents the training set for the Bayesian analysis. Each model consists of M properties; the first ``nphys`` of these are physical properties that are the quantities to be inferred from the observations, while the remaining ones are photometric properties. An important point is that the ``dataset`` object is NOT copied, so altering it after the ``bp`` object is created will result in erroneous results.

The ``priors`` and ``sample_density`` arguments are used to compute the weighting to apply to the input models. The ``priors`` argument specifies the prior probability to assign to each model; it can be either an array giving a prior probability directly, or a callable that can take the physical properties of models as an input and return the prior probability as an output. Similarly, the ``sample_density`` argument specifies the probablity distribution from which the physical models were selected; as with ``priors``, it can be an array or a callable.

The ``bandwidth`` argument specifies the bandwidth to use in the kernel density estimation; this need not be the same in each dimension. The ``bandwidth`` can be specified as a float, in which case it is the same for every dimension, or as an array of M elements giving the bandwidth for every dimension. Finally, it can be set to the string ``auto``, in which case the ``bp`` will attempt to make a reasonable choice of bandwidth autonomously. However, this autonomous choice will probably perform less well than something that is hand-chosen by the user based on their knowledge of the library. As a rule of thumb, bandwidths should be chosen so that, for typical input photometric values, there are ~10 simulations within the 1 kernel size.

Note that both ``priors`` and ``bandwidth`` are properties of the ``bp`` class, and can be altered after the ``bp`` object is created. This makes it possible to alter the priors and bandwidth without incurring the computational or memory cost of generating an entirely new ``bp`` object.


Using ``bp`` Objects
--------------------

Once a ``bp`` object is instantiated, it can be used to compute likelihood functions, marginal probabilities, and MCMC sample ensembles, and to search the library for the best matches to an input set of photometry.

The likelihood function is implemented via the ``bp.logL`` method, which has the call signature::

  def logL(self, physprop, photprop, photerr=None):

The argument ``physprop`` is a set of physical properties, the argument ``photprop`` is a set of photometric properties, and the argument ``photerr`` is an (optional) set of photometric errors. All of these must be arrays, the size of whose trailing dimension matches the number of physical properties (for ``physprop``) or the number of photometric properties (for ``photprop`` and ``photerr``); the leading dimensions of these arrays are broadcast together using normal broadcasting rules. The quantity returned is the log of the joint probability distribution of physical and photometric properties. Specifically, the quantity returned for each input set of physical and photometric properties is

.. math:: \log p(\mathbf{x}, \mathbf{y}, \sigma_{\mathbf{y}}) = \log A \sum_{i=1}^N w_i G(\mathbf{x}, \mathbf{y}; \mathbf{h}')

where :math:`A` is a normalization constant chosen to ensure that the PDF integrated over all space is unity, :math:`\mathbf{x}` is the vector of physical properties, :math:`\mathbf{y}` is the vector of photometric properties, :math:`\sigma_\mathbf{y}` is the vector of photomtric errors, :math:`w_i` is the weight of the ith model as determined by the priors and sample density,

.. math:: G\left(\mathbf{x}, \mathbf{y}; \mathbf{h}'\right) \propto \exp\left[-\left(\frac{x_1^2}{2h_{x_1}'^2} + \cdots + \frac{x_N^2}{2h_{x_N}'^2} + \frac{y_1^2}{2h_{y_1}'^2} + \cdots + \frac{y_M^2}{2h_{y_M}'^2} \right)\right]


is the N-dimensional Gaussian function, and

.. math:: \mathbf{h'} = \sqrt{\mathbf{h}^2 + \sigma_{\mathbf{y}}^2}

is the modified bandwidth, which is equal to the bandwidth used for kernel density estimation added in quadrature sum with the errors in the photometric quantities (see the slug method paper for details).

Estimation of marginal PDFs is done via the ``bp.mpdf`` method, which has the call signature::

  def mpdf(self, idx, photprop, photerr=None, ngrid=128,
           qmin=None, qmax=None, grid=None, norm=True):

The argument ``idx`` is an int or a list of ints between 0 and nphys-1, which specifies for which physical quantity or physical quantities the marginal PDF is to be computed. These indices refer to the indices in the ``dataset`` array that was input when the ``bp`` object was instantiated. The arguments ``photprop`` and ``photerr`` give the photometric measurements and their errors for which the marginal PDFs are to be computed; they must be arrays whose trailing dimension is equal to the number of photometric quantities. The leading dimensions of these arrays are broadcast together following the normal broadcasting rules. By default each physical quantity will be estimated on a grid of 128 points, evenly spaced from the lowest value of that physical property in the model library to the highest value. The parameters ``qmin``, ``qmax``, ``ngrid``, and ``grid`` can be used to override this behavior and set the grid of evaluation points manually. The function returns a tuple ``grid_out, pdf``; here ``grid_out`` is the grid of points on which the marginal PDF has been computed, and ``pdf`` is the value of the marginal PDF evaluated at those gridpoints.

MCMC calculations are implemented via the method ``bp.mcmc``; this method relies on the `emcee <http://dan.iel.fm/emcee/current/>`_ python module, and will only function if it is installed. The call signature is::

  def mcmc(self, photprop, photerr=None, mc_walkers=100,
           mc_steps=500, mc_burn_in=50):

The quantities ``photprop`` and ``photerr`` have the same meaning as for ``bp.mpdf``, and the quantities ``mc_walkers``, ``mc_steps``, and ``mc_burn_in`` are passed directly to ``emcee``, and are described in `emcee's documentation <http://dan.iel.fm/emcee/current/>`_. The quantity returned is an array of sample points computed by the MCMC; its format is also described in `emcee's documentation <http://dan.iel.fm/emcee/current/>`_. Note that, although ``bp.mcmc`` can be used to compute marginal PDFs of the physical quantities, for marginal PDFs of 1 quantity or joint PDFs of 2 quantities it is almost always faster to use ``bp.mpdf`` than ``bp.mcmc``. This is because ``bp.mpdf`` takes advantage of the fact that integrals of cuts through N-dimensional Gaussians can be integrated analytically to compute the marginal PDFs directly, though needing to evaluate the likelihood function point by point. In contrast, the general MCMC algorithm used by ``emcee`` effectively does the integral numerically.

The ``bp.bestmatch`` method searches through the model library and finds the N library entries that are closest to an input set of photometry. The call signature is::

  def bestmatch(self, phot, nmatch=1, bandwidth_units=False):

Here ``phot`` is the set of photometric properties, which is identical to the ``photprop`` parameter used by ``logL``, ``mpdf``, and ``mcmc``. The argument ``nmatch`` specifies how many matches to return, and the argument ``bandwidth_units`` specifies whether distances are to be measured using an ordinary Euclidean metric, or in units of the kernel bandwidth in a given direction. The function returns, for each input set of photometry, the physical and photometric properties of the ``nmatch`` models in the library that are closest to the input photometric values. This can be used to judge if a good match to the input photometry is present in the library.


Caching
-------

The ``bp`` class is built to compute posterior PDFs of some
quantities, marginalising over others. To speed this process, it uses
an internal KD tree representation of the data. By default the KD tree
spans all physical and photometric dimensions of the underlying data
set. When marginalising over some dimensions, however, it can be much
more efficient to have a KD tree that only spans the dimensions that
are not being marginalised, particularly if there are many of
them. For example, suppose that we have a library of sample star
clusters with a range of masses, ages, extinctions, metallicities, and
photometric measurements. We are interested in the posterior PDF of
mass given a set of photometry, marginalising over age, extinction,
and metallicity. In this case it is more efficient to have a KD tree
that does not split in the age, extinction, or metallicity
dimensions. Since evaluating the marginal posterior PDF of mass using
such a tree is much faster, if we are going to compute marginal
posterior PDFs of mass many times (for example for many photometric
measurements) it is advantageous to pay the one-time cost of
constructing the better-optimised KD tree and then use it for all
subsequent calculations.

To handle this, ``bp`` has an optional caching capability that will
cache KD tree representations of the data that are optimal for
computing posterior PDFs marginalised over certain dimensions. One can
construct such a cached tree by invoking the ``bp.make_cache``
method. The syntax is simple::

  bp.make_cache(margindims)

where ``margindims`` is a listlike containing the dimensions that will
be marginalized out. Once created, the cache will be used for all
subsequent evaluations using ``bp.mpdf`` and related methods
(``mpdf_phot`` and ``mpdf_gen``) that marginalise over those
dimensions.

Caching can also be done automatically. The ``bp`` constructor accepts
the keyword ``caching`` which specifies the type of caching to
use. The default value, ``none``, uses no caching. The value ``lazy``
causes a cached KD tree to be build whenever one of the marginal PDF
methods is invoked, and is stored for future use. (Warning: ``lazy``
mode is not thread-safe -- see :ref:`ssec-bayesphot-threading`.)
Finally, ``aggressive`` mode constructs caches for all possible
one-dimensional marginalisations of physical variables given observed
photometry, and all possible one-dimensional marginalisations of
variables by themselves, when the ``bp`` object is first constructed.


.. _ssec-bayesphot-threading:

Parallelism in bayesphot
------------------------

The ``bp`` class supports parallel calculations of posterior PDFs and
related quantities, through the python `multiprocessing module
<https://docs.python.org/2.7/library/multiprocessing.html>`_. This
allows efficient use of multiple cores on a shared memory machine,
circumventing the python global interpreter lock, without the need for
every process to read a large simulation library or store it in
memory. The recommended method for writing threaded code using ``bp``
objects is to use have a master process create the ``bp`` object, and
then use a `Process
<https://docs.python.org/2.7/library/multiprocessing.html#multiprocessing.Process>`_
or `Pool
<https://docs.python.org/2.7/library/multiprocessing.html#module-multiprocessing.pool>`_
object to create child processes the perform computations using ``bp``
methods such as ``bp.logL`` or ``bp.mpdf``. It is often most efficient
to combine this with shared memory objects such as `RawArray
<https://docs.python.org/2.7/library/multiprocessing.html#module-multiprocessing.sharedctypes>`_
to hold the outputs.

An example use case for computing 1D marginal PDFs on a large set of photometric values is::

  # Import what we need
  from slugpy.bayesphot import bp
  from multiprocessing import Pool, RawArray
  from ctypes import c_double
  import numpy as np

  # Some code here to create / read the data set to be used by
  # bayesphot and store it in a variable called dataset

  # Create the bayesphot object
  my_bp = bp(dataset, nphys)

  # Some code here to create / read the photometric data we want to
  # process using bayesphot and store it in an array called phot,
  # which is of shape (nphot, nfilter). There is also an array of
  # photometric errors, called photerr, of the same shape.

  # Create a holder for the output
  pdf_base = RawArray(c_double, 128*nphot)
  pdf = np.frombuffer(pdf_base, dtype=c_double). \
        reshape((nphot,128))
  grd_base = RawArray(c_double, 128)
  grd = np.frombuffer(grd_base, dtype=c_double)

  # Define the function that will be responsible for computing the
  # marginal PDF
  def mpdf_apply(i):
      grd[:], pdf[i] = my_bp.mpdf(0, phot[i], photerr[i])

  # Main thread starts up a process pool and starts the computation
  if __name__ == '__main__':
      pool = Pool()
      pool.map(mpdf_apply, range(nphot))

      # At this point the grd and PDF contain the same as they would
      # if we had done
      #     grd, pdf = my_bp.mpdf(0, phot, photerr)
      # but the results will be computed much faster this way

For an example of a more complex use case, see the `LEGUS cluster pipeline <https://bitbucket.org/krumholz/legus-cluster-pipeline/overview>`_.

The full list of ``bp`` methods that are thread-safe is:

* ``bp.logL``
* ``bp.mpdf``
* ``bp.mcmc``
* ``bp.bestmatch``
* ``bp.make_approx_phot``
* ``bp.make_approx_phys``
* ``bp.squeeze_rep``
* ``bp.mpdf_approx``

Thread safety involves a very modest overhead in terms of memory and speed, but for non-threaded computations this can be avoided by specifying::

  bp.thread_safe = False

or by setting the ``thread_safe`` keyword to ``False`` when the ``bp``
object is constructed.

Finally, note that this parallel paradigm avoids duplicating the large
library only on unix-like operating systems that support copy-on-write
semantics for `fork
<https://en.wikipedia.org/wiki/Fork_(system_call)>`_. Code such as the
above example should still work on windows (though this has not been
tested), but each worker process will duplicate the library in
physical memory, thereby removing one of the main advantages of
working in parallel.


.. _ssec-slugpy-bayesphot:

Full Documentation of slugpy.bayesphot
--------------------------------------

.. automodule:: slugpy.bayesphot.bp
   :members:
   :special-members:

