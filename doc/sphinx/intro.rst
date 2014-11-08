.. highlight:: rest

Introduction to SLUG
====================

This is a guide for users of the SLUG software package. SLUG is distributed under the terms of the `GNU General Public License v. 3.0 <http://www.gnu.org/licenses/gpl.html>`_. A copy of the license notification is included in the main SLUG directory. If you use SLUG in any published work, please cite the SLUG method paper, `da Silva, R. L., Fumagalli, M., & Krumholz, M. R., 2012, The Astrophysical Journal, 745, 145 <http://adsabs.harvard.edu/abs/2012ApJ...745..145D>`_. A second method paper, describing the upgraded version 2 code and a set of ancillary tools, is in preparation at this time.

What Does SLUG Do?
------------------

SLUG is a stellar population synthesis (SPS) code, meaning that, for a specified stellar initial mass function (IMF), star formation history (SFH), cluster mass function (CMF), cluster lifetime function (CLF), and (optionally) distribution of extinctions (A_V), it predicts the spectra and photometry of both individual star clusters and the galaxies (or sub-regions of galaxies) that contain them. In this regard, SLUG operates much like any other SPS code. The main difference is that SLUG regards the functions describing the stellar population as probability distributions, and the resulting stellar population as being the result of a draw from them. SLUG performs a Monte Carlo simulation to determine the PDF of the light produced by the stellar populations that are drawn from these distributions. The remainder of this section briefly describes the major conceptual pieces of a SLUG simulation. For a more detailed description, readers are referred to `da Silva, Fumagalli, & Krumholz (2012) <http://adsabs.harvard.edu/abs/2012ApJ...745..145D>`_.

Cluster Simulations and Galaxy Simulations
------------------------------------------

SLUG can simulate either a simple stellar population (i.e., a group of stars all born at one time) or a composite stellar population, consisting of stars born at a distribution of times. We refer to the former case as a "cluster" simulation, and the latter as a "galaxy" simulation, since one can be thought of as approximating the behavior of a single star cluster, and the other as approximating a whole galaxy.

.. _ssec-slugpdfs:

Probability Distribution Functions: the IMF, SFH, CMF, CLF, A_V distribution
----------------------------------------------------------------------------

As mentioned above, SLUG regards the IMF, SFH, CMF, CLF, and extinction A_V as probability distribution functions. These PDFs can be described by a very wide range of possible functional forms; see :ref:`sec-pdfs` for details on the exact functional forms allowed, and on how they can be specified in the code. When SLUG runs a cluster simulation, it draws stars from the specified IMF in an attempt to produce a cluster of a user-specified total mass. There are a number of possible methods for performing such mass-limited sampling, and SLUG gives the user a wide menu of options; see :ref:`sec-pdfs`. SLUG will also, upon user request, randomly draw a visual extinction A_V to be applied to the light.

For a galaxy simulation, the procedure involves one extra step. In this case, SLUG assumes that some fraction :math:`f_c` of the stars in the galaxy are born in star clusters, which, for the purposes of SLUG, means that they all share the same birth time. The remaining fraction :math:`1-f_c` of stars are field stars. When a galaxy simulation is run, SLUG determines the total mass of stars :math:`M_*` that should have formed since the start of the simulation (or since the last output, if more than one output is requested) from the star formation history, and then draws field stars and star clusters in an attempt to produce masses :math:`(1-f_c)M_*` and :math:`f_c M_*`. For the field stars, the stellar masses are drawn from the IMF, in a process completely analogous to the cluster case, and each star is given its own randomly-generated extinction. For star clusters, the masses of the clusters are drawn from the CMF, and each cluster is then populated from the IMF as in the cluster case. Again, each cluster gets its own extinction. For both the field stars and the star clusters, the time of their birth is drawn from the PDF describing the SFH.

Finally, star clusters can be disrupted independent of the fate of their parent stars. When each cluster is formed, it is assigned a lifetime drawn from the CLF. Once that time has passed, the cluster ceases to be entered in the lists of individual cluster spectra and photometry (see next section), although the individual stars continue to contribute to the integrated light of the galaxy.

.. _ssec-spec-phot:

Spectra, Photometry, and Extinction
-----------------------------------

Once SLUG has drawn a population of stars, its final step is to compute the light they produce. SLUG does this in several steps. First, it computes the physical properties of all the stars present user-specified times using a set of stellar evolutionary tracks. Second, it uses these physical properties to compute the composite spectra produced by the stars, using a user-specified set of stellar atmosphere models. Formally, the quantity computed is the specific luminosity per unit wavelength :math:`L_\lambda`. Third, it computes photometry for the stellar population by integrating the computed spectra over a set of specified photometric filters. If extinction is enabled, photometric values are computed for both the unextincted and the extincted spectrum. Depending on the options specified by the user and the filter under consideration, the photometric value output will be one of the following:

* The frequency-averaged luminosity across the filter, defined as

.. math:: \langle L_\nu\rangle_R = \frac{\int L_\nu \, d\ln\nu}{\int R_\nu (\nu/\nu_c)^\beta \, d\ln\nu},

where :math:`L_\nu` is the specific luminosity per unit frequency, :math:`R_\nu` is the filter response function per photon at frequency :math:`\nu`, :math:`\nu_c` is the central wavelength of the filter, and :math:`\beta` is a constant that is defined by convention for each filter, and is either 0, 1, or 2; usually it is 0 for optical and UV filters.

* The wavelength-averaged luminosity across the filter, defined as

.. math::\langle L_\lambda\rangle_R = \frac{\int L_\lambda \, d\ln\lambda}{\int R_\lambda (\lambda/\lambda_c)^{-\beta} \, d\ln\lambda},

where :math:`L_\lambda` is the specific luminosity per unit wavelength, :math:`R_\lambda` is the filter response function per photon at wavelength :math:`\lambda`, and :math:`\lambda_c` is the central wavelength of the filter.

* The AB magnitude, defined by

.. math:: M_{\rm AB} = -2.5 \log_{10} \left[\frac{\langle L_\nu\rangle_R}{4\pi\left(10\,\mathrm{pc}\right)^2}\right] - 48.6,

where :math:`\langle L_\nu\rangle_R` is in units of :math:`\mathrm{erg\,s}^{-1}\,\mathrm{Hz}^{-1}`.

* The ST magnitude, defined by

.. math:: M_{\rm ST} = -2.5 \log_{10} \left[\frac{\langle L_\lambda\rangle_R}{4\pi\left(10\,\mathrm{pc}\right)^2}\right] - 21.1,

where :math:`\langle L_\lambda\rangle_R` is in units of :math:`\mathrm{erg\, s}^{-1}\,\mathrm{Angstrom}^{-1}`.

* The Vega magnitude, defined by

.. math:: M_{\rm Vega} = M_{\rm AB} - M_{\rm AB}(\mbox{Vega}),

where :math:`M_{\rm AB}(\mbox{Vega})` is the AB magnitude of Vega. The latter quantity is computed on the fly, using a stored Kurucz model spectrum for Vega. 

* The photon flux above some threshold :math:`\nu_0`, defined as

.. math:: Q(\nu_0) = \int_{\nu_0}^\infty \frac{L_\nu}{h\nu} \, d\nu.

* The bolometric luminosity,

.. math:: L_{\rm bol} = \int_0^\infty L_\nu \, d\nu.

Fourth and finally, if extinction is enabled, SLUG also computes an extincted spectrum

.. math:: L_{\lambda,\mathrm{ex}} = L_{\lambda} e^{-\tau_\lambda}

where the optical depth :math:`\tau_\lambda = (\kappa_\lambda / \kappa_V) (A_V/1.086)`, :math:`A_V` is the visual extinction in mag, the factor 1.086 is the conversion between magnitudes and the true dimensionless optical depth, :math:`\kappa_\lambda` is a user-specified input extinction at wavelength :math:`\lambda`, and the V-band mean opacity is defined by

.. math:: \kappa_V = \frac{\int \kappa_\nu R_\nu(V) \, d\nu}{\int R_\nu(V) \, d\nu}

where :math:`R_\nu(V)` is the filter response function as frequency :math:`\nu` for the Johnson V filter. Photometric values for the extincted spectrum are computed exactly as for the unextincted one.

For a cluster simulation, this procedure is applied to the star cluster being simulated at a user-specified set of output times. For a galaxy simulation, the procedure is much the same, but it can be done both for all the stars in the galaxy taken as a whole, and individually for each star cluster that is still present (i.e., that has not been disrupted).

Monte Carlo Simulation
----------------------

The steps described in the previous two section are those required for a single realization of the stellar population. However, the entire point of SLUG is to repeat this procedure many times in order to build up the statistics of the population light output. Thus the entire procedure can be repeated as many times as the user desires.
