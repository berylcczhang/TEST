.. highlight:: rest

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


.. _ssec-mpi-support:

Support for MPI Parallelism
---------------------------

In large codes where one might wish to use SLUG for subgrid stellar
models, it is often necessary to pass information between processors
using MPI. Since SLUG's representation of stellar populations is
complex, and much information is shared between particles rather than
specific to individual particles (e.g., tables of yields and
evolutionary tracks), passing SLUG information between processors is
non-trivial.

To facilitate parallel implementations, SLUG provides routines that
wrap the base MPI routines and allow seamless and efficient exchange
of the slug_cluster class (which SLUG uses to represent simple stellar
populations) between processors. The prototypes for these functions
are found in the ``src/slug_MPI.H`` header file.

By default MPI support is not included in the library. To enable MPI
support, compile the library as follows::

  make lib MPI=ENABLE_MPI

This will enable MPI support. In addition, you may need to specify the
names of your preferred MPI C++ compiler by setting the variable
``MACH_MPICXX`` in your machine-specific makefile -- see
:ref:`ssec-machine-makefiles`. The Makefiles contain reasonable
guesses, but since MPI compiler names are much less standardized than
general compiler names, you may need to supply yours rather than
relying on the default.
