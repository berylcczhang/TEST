.. highlight:: rest
	       
.. _sec-library-mode:

Using SLUG as a Library
=======================

In addition to running as a standalone program, SLUG can be
compiled as a library that can be called by external programs. This is
useful for including stellar population synthesis calculations within
some larger code, e.g., a galaxy simulation code in which star
particles represent individual star clusters, where the stars in them
are treated stochastically. 

.. _ssec-library-mode:

Compiling in Library Mode
-------------------------

To compile in library mode, simply do::

  make lib

in the main directory. This will cause a dynamically linked library
file ``libslug.x`` to be created in the ``src`` directory, where ``x``
is whatever the standard extension for dynamically linked libraries on
your system is (``.so`` for unix-like systems, ``.dylib`` for MacOS).

Alternately, if you prefer a statically-linked version, you can do::

  make libstatic

and a statically-linked archive ``libslug.y`` will be created instead,
where ``y`` is the standard statically-linked library extension on
your system (generally ``.a``).

In addition to ``lib`` and ``libstatic``, the makefile supports
``lib-debug`` and ``libstatic-debug`` as targets as well. These
compile the same libraries, but with optimization disabled and
debugging symbols enabled.

Finally, if you want MPI functionality, you can compile with::

  make lib MPI=ENABLE_MPI

See :ref:`ssec-compiling` for more on compiling with MPI enabled.


.. _ssec-predefined-objects:

Predefined Objects
------------------

In order to make it more convenient to use slug as a library, the
library pre-defines some of the most commonly-used classes, in order
to save users the need to construct them. These predefined objects can
be accessed by including the file ``slug_predefined.H`` in your source
file. This function defines the class ``slug_predef``, which
pre-defines all the IMFs, evolutionary tracks, spectral synthesizers,
and yields that ship with slug, without forcing the user to interact
with the parameter parsing structure.

The ``slug_predef`` class provides the methods ``imf``, ``tracks``,
``specsyn``, and ``yields``. These methods take as arguments a string
specifying one of the predefined names of an IMF, set of tracks, or
spectral synthesizer, and return an object of that class that can then
be passed to ``slug_cluster`` to produce a cluster object. For
example, the following sytax creates a ``slug_cluster`` with ID number
1, a mass of 100 solar masses, age 0, a Chabrier IMF, Padova solar
metallicity tracks, starburst99-style spectral synthesis, and slug's
default nuclear yields::

  #include "slug_predefined.H"
  #include "slug_cluster.H"
  
  slug_cluster *cluster =
     new slug_cluster(1, 100.0, 0.0, slug_predef.imf("chabrier"),
                      slug_predef.tracks("modp020.dat"),
		      slug_predef.specsyn("sb99"),
		      nullptr, nullptr, nullptr,
		      slug_predef.yields());

.. _ssec-mpi-support:

Using SLUG as a Library with MPI-Enabled Codes
----------------------------------------------

In large codes where one might wish to use slug for subgrid stellar
models, it is often necessary to pass information between processors
using MPI. Since slug's representation of stellar populations is
complex, and much information is shared between particles rather than
specific to individual particles (e.g., tables of yields and
evolutionary tracks), passing slug information between processors is
non-trivial.

To facilitate parallel implementations, slug provides routines that
wrap the base MPI routines and allow seamless and efficient exchange
of the slug_cluster class (which slug uses to represent simple stellar
populations) between processors. The prototypes for these functions
are found in the ``src/slug_MPI.H`` header file, and the functions are
available if the library was compiled with MPI support enabled (see
:ref:`ssec-library-mode`).

Here is an example of MPI usage, in which one processor creates a
cluster and then sends it to another one::

  #include "slug_cluster.H"
  #include "slug_MPI.H"
  #include "mpi.h"
  #include <vector>
  #include <cstdio>

  int main(int argc, char *argv[]) {

    // Start MPI
    MPI_Init(&argc, &argv);

    // Get rank
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Rank 0 creates a cluster and prints out the masses of the stars
    slug_cluster *cluster;
    if (rank == 0) {
      cluster =
         new slug_cluster(1, 100.0, 0.0, slug_predef.imf("chabrier"),
                          slug_predef.tracks("modp020.dat"),
	    	          slug_predef.specsyn("sb99"),
		          nullptr, nullptr, nullptr,
		          slug_predef.yields());
      const std::vector<double> stars = cluster->get_stars();
      for (int j=0; j<stars.size(); j++)
	std::cout << "rank 0, star " << j
		  << ": " << stars[j] << std::endl;
    }
    
    // Barrier to make sure rank 0 outputs come first
    MPI_Barrier(MPI_COMM_WORLD);

    // Rank 0 sends cluster, rank 1 receives it
    if (rank == 0) {
      MPI_send_slug_cluster(*cluster, 1, 0, MPI_COMM_WORLD);
    } else if (rank == 1) {
      cluster = MPI_recv_slug_cluster(0, 1, MPI_COMM_WORLD,
                                      slug_predef.imf("chabrier"),
                                      slug_predef.tracks("modp020.dat"),
	    	                      slug_predef.specsyn("sb99"),
		                      nullptr, nullptr, nullptr,
		                      slug_predef.yields());
    }

    // Rank 1 prints the masses of the stars; the resulting masses
    // should be identical to that produced on rank 0
    if (rank == 1) {
      const std::vector<double> stars = cluster->get_stars();
      for (int j=0; j<stars.size(); j++)
	std::cout << "rank 1, star " << j
		  << ": " << stars[j] << std::endl;
    }
  }
