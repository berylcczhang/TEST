.. highlight:: rest

Compiling and Installing SLUG
=============================

Dependencies
------------

The core SLUG program requires

* The `Boost C++ libraries <http://www.boost.org/>`_
* The `GNU scientific library <http://www.gnu.org/software/gsl/>`_ (version 2.x preferred, code can be compiled with version 1.x -- see below)
* The `cfitsio library <http://heasarc.gsfc.nasa.gov/fitsio/fitsio.html>`_ (optional, only required for FITS capabilities)

Compilation will be easiest if you install these libraries such that the header files are included in your ``CXX_INCLUDE_PATH`` (for Boost) and ``C_INCLUDE_PATH`` (for GSL) and the compiled object files are in your ``LIBRARY_PATH``. Alternately, you can manually specify the locations of these files by editing the Makefiles -- see below. The cfitsio library is optional, and is only required if you want the ability to write FITS output. To compile without it, use the flag ``FITS=DISABLE_FITS`` when calling ``make`` (see below). Note that SLUG uses some Boost libraries that must be built separately (see the Boost documentation on how to build and install Boost libraries).

In addition to the core dependencies, slugpy, the python helper library requires:

* `numpy <http://www.numpy.org/>`_
* `scipy <http://www.scipy.org/>`_
* `astropy <http://www.astropy.org/>`_ (optional, only required for FITS capabilities)

Finally, the cloudy coupling capability requires:

* `cloudy <http://nublado.org>`_

This is only required performing cloudy runs, and is not required for any other part of SLUG.

Compiling
---------

If you have Boost in your ``CXX_INCLUDE_PATH``, GSL in your ``C_INCLUDE_PATH``, and (if you're using it) cfitsio in your ``C_INCLUDE_PATH``, and the compiled libraries for each of these in your ``LIBRARY_PATH`` environment variables, and your system is running either MacOSX or Linux, you should be able to compile simply by doing::

   make

from the main ``slug`` directory.

To compile in debug mode, do::

   make debug

instead. 

If you are compiling using GSL version 1.x or without cfitsio, you must specify these options when compiling. If you are using version 1.x of the GSL, do::

  make GSLVERSION=1

To compile without cfitsio, do::

   make FITS=DISABLE_FITS

Alternately, you can manually specify the compiler flags to be used by creating a file named ``Make.mach.MACHINE_NAME`` in the ``src`` directory, and then doing::

   make MACHINE=MACHINE_NAME

An example machine-specific file, ``src/Make.mach.ucsc-hyades`` is included in the repository. You can also override or reset any compilation flag you want by editing the file ``src/Make.config.override``.

Finally, note that SLUG is written in C++11, and requires some C++11 features, so it may not work with older C++ compilers. The following compiler versions are known to work: gcc >= 4.8 (4.7 works on most but not all platforms), clang/llvm >= 3.3, icc >= 14.0. Earlier versions may work as well, but no guarantees.

Using SLUG as a Library
-----------------------

In addition to using SLUG as a standalone program, SLUG can be compiled as a library that can be called by external programs. This is useful for including stellar population synthesis calculations within some larger code, e.g., a galaxy simulation code. To compile in library mode, simply do::

  make lib

in the main directory. This will cause a dynamically linked library file ``libslug.x`` to be created in the ``src`` directory, where ``x`` is whatever the standard extension for dynamically linked libraries on your system is (``.so`` for unix-like systems, ``.dylib`` for MacOS).

Alternately, if you prefer a statically-linked version, you can do::

  make libstatic

and a statically-linked archive ``libslug.y`` will be created instead, where ``y`` is the standard statically-linked library extension on your system (generally ``.a``).

In addition to ``lib`` and ``libstatic``, the makefile supports ``lib-debug`` and ``libstatic-debug`` as targets as well. These compile the same libraries, but with optimization disabled and debugging symbols enabled.
