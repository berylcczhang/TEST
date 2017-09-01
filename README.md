### What is SLUG ###

SLUG is the Stochastically Lighting Up Galaxies Code. SLUG is a stellar population synthesis code, and in many respects is similar to other commonly-used SPS codes such as [starburst99](http://www.stsci.edu/science/starburst99/docs/default.htm) and [FSPS](https://code.google.com/p/fsps/). Given an input star formation history, stellar initial mass function, and a set of evolutionary tracks and stellar atmospheres, all of these codes can predict the spectra and photometric properties of a stellar population, and some can also predict chemical yields. The main difference between SLUG and conventional SPS codes is that, instead of the usual approach of assuming that all stellar masses and ages are fully populated, SLUG is capable of stochastically sampling from the stellar initial mass function and age distribution, and thereby predicting not just the mean spectrum, photometry, and yields, but also the full distribution of these quantities that results from stochastic sampling. This capability is critical in the regime of low star formation rates and total stellar masses, where finite sampling can lead to a distribution of properties that is extremely broad, and the mean values produced by other SPS codes are therefore of limited predictive power.

### Full documenation ###

Full documentation of SLUG is available at <http://slug2.readthedocs.io/en/latest/>. Documentation is also available in the source repository; see *doc/html/index.html* for the html version, and *doc/latex/slug.pdf* for a PDF version.

### Version history ###

This repository contains SLUG v2. SLUG v1 is available at <http://www.slugsps.com/downloads/v1>. The code here is a complete re-implementation of SLUG, with a number of major improvements.

* The ability to predict full spectra as well as photometry
* A vastly larger array of photometric filters, with the ability to add more without needing to alter the source code; there is also wide range of choices for photometric system
* Full support for both thread-based or MPI-based parallelism
* Predictions of element yields as well as light output
* The ability to handle a very wide range of functional forms for the initial mass function, cluster mass function, cluster lifetime function, and star formation history
* Many sampling methods available for mass-limited sampling
* The ability to simulate simple stellar populations with finite masses, as well as full star formation histories
* More accurate numerical methods for interpolation and integration
* A cleaner, easier-to-use control file structure
* Greater flexibility in output and output formatting (ASCII, binary, or FITS)
* Python helper routines for parsing output (replacing the older IDL routines used by SLUG v1)
* The ability to apply stochastically-selected extinctions to the predicted light output, including diffential stellar and nebular extinction
* Quick estimation of the nebular contribution to spectra and photometry
* An automated interface to couple the code to [cloudy](http://nublado.org) in order to provide much more accurate predictions of nebular emission
* A suite of Bayesian estimator tools that can be used to derive posterior probabilities on star formation rates and star cluster properties from libraries of SLUG simulations.

SLUG v1 will no longer be maintained, so all users are encouraged to migrate to SLUG v2 as soon as possible.

### Layout of the repository ###

The SLUG repository contains several subdirectories, whose contents are as follows:

* *bin*: this contains the slug source code when it is compiled, as well as the parallel wrapper python scripts
* *cloudy_slug*: code for the interface to [cloudy](http://nublado.org)
* *cluster_slug*: example code and data for cluster_slug
* *doc*: documentation for the code
    - *doc/html*: documentation in html format
    - *doc/latex*: documentation in latex/pdf format
    - *doc/sphinx*: source for documentation in sphinx/rst format
* *lib*: data files required by the code
    - *lib/atmospheres*: stellar atmosphere models
    - *lib/atomic*: atomic data files
    - *lib/avdist*: predefined distributions of visual extinction
    - *lib/clf*: predefined star cluster lifetime functions
    - *lib/cmf*: predefined star cluster mass fuctions
    - *lib/extinct*: predefined extinction curves
    - *lib/filters*: predefined photometric filters
    - *lib/imf*: predefined initial mass functions
    - *lib/tracks*: stellar evolution tracks
    - *lib/yields*: chemical yield tables
* *param*: default location for parameter files that control slug simulations; example files are provided
* *sfr_slug*: example code and data for sfr_slug
* *slugpy*: this contains the python module slugpy, which includes routines for parsing and manipulating simulation outputs
    - *slugpy/bayesphot*: python routines to do Bayesian estimation of physical properties from photometry
      - *slugpy/bayesphot/bayesphot_c*: source code for c extensions to slugpy/bayesphot
    - *slugpy/cloudy*: parsing routines related to the [cloudy](http://nublado.org) interface
    - *slugpy/cluster_slug*: tools for using a suite of SLUG simulations to do Bayesian inference of star cluster properties
    - *slugpy/sfr_slug*: tools for using a suite of SLUG simulations to do Bayesian inference of star formation rates
* *src*: the main source code for SLUG
* *test*: test scripts
* *tools*: utility tools for manipulating input and output files

### Installing and configuring ###

#### Downloading ####

The easiest way to get SLUG is via git. From the command line, just do

`git clone https://bitbucket.org/krumholz/slug2.git`

You can also download the entire repository from the [SLUG downloads page](https://bitbucket.org/krumholz/slug2/downloads), but git is the preferred method.

#### Dependencies ####

The core SLUG c++ code requires:

* The [GNU scientific library](http://www.gnu.org/software/gsl/) (version 2.x preferred, 1.x can be used with compiler flag)
* The [BOOST C++ libraries](http://www.boost.org/)
* The [cfitsio libary](http://heasarc.gsfc.nasa.gov/fitsio/fitsio.html) (optional, only required for FITS output capability)

The slugpy python routines require:

* [numpy](http://www.numpy.org/)
* [scipy](http://www.scipy.org/)
* [astropy](http://www.astropy.org/) (optional, only required for FITS handling capability)

Support for MPI-based parallelism requires:

* [MPI](http://mpi-forum.org/) (an MPI implementation that supports the MPI 3 standard required for full functionality; limited functionality with earlier versions)

The interface to cloudy (obviously) requires:

* [cloudy](http://nublado.org) (not required for core functionality)

#### Compiling ####

SLUG comes with a Makefile in the main directory, and it should be possible to build the code in most environments simply by typing "make". Compliation requires that the GSL and BOOST header files be included in the include path, and that their shared object libraries be included in the link path. If FITS capability is desired, the cfitsio library and header must be included as well. If MPI capability is desired, and an MPI library is installed, type "make MPI=ENABLE_MPI".

#### Running ####

Once the code is compiled, running a SLUG simulation is fairly straightforward. The required steps are:

1. Set the environment variable SLUG_DIR to the directory where you have installed SLUG. This is not strictly necessary, but it avoids the need to specify manually the locations of all the data files that SLUG needs. If you plan to use the [cloudy](http://nublado.org) interface, you should also set the environment variable CLOUDY_DIR to the directory where the cloudy executable is located.
2. Create a parameter file to set up the desired simulation. The files *param/example.param* and *param/example_cluster.param* can provide useful starting points for creating your own parameter files. Alternately, you can use the menu-driven parameter file generator script to create a parameter file for you by doing `python ./tools/write_param.py`.
3. Run the simulation with the command line `./bin/slug param/mysimulation.param`. Alternately, you can run the simulation in parallel by doing `python ./bin/slug.py param/mysimulation.param` (for thread-based parallelism) or `mpirun -np 1 python ./bin/slug_mpi.py param/mysimulation.param` (for MPI-based parallelism).
4. Once the simulation is done, you can examine the output by eye (if you chose to write the output in ASCII format) or, more likely, read the data using the provided python routines. The two basic output reading routines can be invoked by doing the following in a python program of the command prompt of an interactive session:

`from slugpy import *`

`int_data = read_integrated('RUN_NAME')`

`cluster_data = read_cluster('RUN_NAME')`

### Questions, bugs, and getting involved ###

If you have questions about SLUG, have discovered any bugs, or want to contribute to ongoing development, please contact [Mark Krumholz](http://www.mso.anu.edu.au/~krumholz/), mark.krumholz@anu.edu.au.

### Acknowledgements ###

SLUG was the product of many people's work. Version 1 was written by Robert da Silva and Michele Fumagalli, with contributions and improvements from Jonathan Parra. Version 2 of SLUG, as well as slugpy, were written by primarily by Mark Krumholz, with contributions from Michele Fumagalli, Teddy Rendahl, Evan Demers, and Greg Ashworth. See the [documentation](http://slug2.readthedocs.io/en/latest/) for a full list of contributors, including both those who wrote parts of the code and who provided tabulated data.

### License ###

SLUG is distributed under the terms of the [GNU General Public License](http://www.gnu.org/copyleft/gpl.html) version 3.0. The text of the license is included in the main directory of the repository as GPL-3.0.txt.