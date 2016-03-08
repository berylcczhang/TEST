.. highlight:: rest

Introduction to SLUG
====================

This is a guide for users of the SLUG software package. SLUG is distributed under the terms of the `GNU General Public License v. 3.0 <http://www.gnu.org/licenses/gpl.html>`_. A copy of the license notification is included in the main SLUG directory. If you use SLUG in any published work, please cite the SLUG method papers, `da Silva, R. L., Fumagalli, M., & Krumholz, M. R., 2012, The Astrophysical Journal, 745, 145 <http://adsabs.harvard.edu/abs/2012ApJ...745..145D>`_ and `Krumholz, M. R., Fumagalli, M., da Silva, R. L., Rendahl, T., & Parra, J. 2015, Monthly Notices of the Royal Astronomical Society, 452, 1447 <http://adsabs.harvard.edu/abs/2015MNRAS.452.1447K>`_.

What Does SLUG Do?
------------------

SLUG is a stellar population synthesis (SPS) code, meaning that, for a specified stellar initial mass function (IMF), star formation history (SFH), cluster mass function (CMF), cluster lifetime function (CLF), and (optionally) distribution of extinctions (A_V), it predicts the spectra and photometry of both individual star clusters and the galaxies (or sub-regions of galaxies) that contain them. It also predicts the yields of various isotopes. In this regard, SLUG operates much like any other SPS code. The main difference is that SLUG regards the functions describing the stellar population as probability distributions, and the resulting stellar population as being the result of a draw from them. SLUG performs a Monte Carlo simulation to determine the PDF of the light and yields produced by the stellar populations that are drawn from these distributions. The remainder of this section briefly describes the major conceptual pieces of a SLUG simulation. For a more detailed description, readers are referred to `da Silva, Fumagalli, & Krumholz (2012) <http://adsabs.harvard.edu/abs/2012ApJ...745..145D>`_.

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

Spectra and Photometry
----------------------

Once SLUG has drawn a population of stars, its final step is to compute the light they produce. SLUG does this in several steps. First, it computes the physical properties of all the stars present user-specified times using a set of stellar evolutionary tracks. Second, it uses these physical properties to compute the composite spectra produced by the stars, using a user-specified set of stellar atmosphere models. Formally, the quantity computed is the specific luminosity per unit wavelength :math:`L_\lambda`. Third, if nebular emission is enabled, the code calculates the spectrum :math:`L_{\lambda,\mathrm{neb}}` that emerges after the starlight passes through the HII region aruond the star -- see :ref:`ssec-nebula`. Fourth, if extinction is enabled, SLUG computes the extincted stellar and nebula-processed spectra :math:`L_{\lambda,\mathrm{ex}}` and :math:`L_{\lambda,\mathrm{neb,ex}}` -- see :ref:`ssec-extinction`. Fifth and finally, SLUG computes photometry for the stellar population by integrating all computed spectra over a set of specified photometric filters. Depending on the options specified by the user and the filter under consideration, the photometric value output will be one of the following:

* The frequency-averaged luminosity across the filter, defined as

.. math:: \langle L_\nu\rangle_R = \frac{\int L_\nu R_\nu \, d\ln\nu}{\int R_\nu (\nu/\nu_c)^\beta \, d\ln\nu},

where :math:`L_\nu` is the specific luminosity per unit frequency, :math:`R_\nu` is the filter response function per photon at frequency :math:`\nu`, :math:`\nu_c` is the central wavelength of the filter, and :math:`\beta` is a constant that is defined by convention for each filter, and is either 0, 1, or 2; usually it is 0 for optical and UV filters.

* The wavelength-averaged luminosity across the filter, defined as

.. math:: \langle L_\lambda\rangle_R = \frac{\int L_\lambda R_\lambda \, d\ln\lambda}{\int R_\lambda (\lambda/\lambda_c)^{-\beta} \, d\ln\lambda},

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

If nebular processing and/or extinction are enabled, photometric quantities are computed separately for each available version of the spectrum, :math:`L_\lambda`, :math:`L_{\lambda,\mathrm{neb}}`, :math:`L_{\lambda,\mathrm{ex}}`, and :math:`L_{\lambda,\mathrm{neb,ex}}`.

For a cluster simulation, this procedure is applied to the star cluster being simulated at a user-specified set of output times. For a galaxy simulation, the procedure is much the same, but it can be done both for all the stars in the galaxy taken as a whole, and individually for each star cluster that is still present (i.e., that has not been disrupted).

Monte Carlo Simulation
----------------------

The steps described in the previous two section are those required for a single realization of the stellar population. However, the entire point of SLUG is to repeat this procedure many times in order to build up the statistics of the population light output. Thus the entire procedure can be repeated as many times as the user desires.

.. _ssec-nebula:

Nebular Processing
------------------

SLUG includes methods for post-processing the output starlight to compute the light that will emerge from the HII region around star clusters, and to further apply extinction to that light.

Nebular emission is computed by assuming that, for stars / star clusters younger than 10 Myr, all the ionizing photons are absorbed in a uniform-density, uniform-temperature HII region around each star cluster / star, and then computing the resulting emission at non-ionizing energies. The calculation assumes that the HII region is in photoionization equilibrium, and consists of hydrogen that is fully ionized and helium that is singly ionized. Under these assumptions the volume :math:`V`, electron density :math:`n_e`, and hydrogen density :math:`n_{\mathrm{H}}` are related to the hydrogen ionizing luminosity :math:`Q(\mathrm{H}^0)` via

.. math:: \phi Q(\mathrm{H}^0) = \alpha_{\mathrm{B}}(T) n_e n_{\mathrm{H}} V

Here :math:`\phi` is the fraction of ionizing photons that are absorbed by hydrogen rather than dust grains, and :math:`\alpha_{\mathrm{B}}(T)` is the temperature-dependent case B recombination rate coefficient. SLUG approximates :math:`\alpha_{\mathrm{B}}(T)` using the analytic approximation given by equation 14.6 of `Draine (2011, Physics of the Interstellar and Intergalactic Medium, Princeton University Press) <http://adsabs.harvard.edu/abs/2011piim.book.....D>`_, while :math:`\phi` is a user-chosen parameter. The temperature can either be set by the user directly, or can be looked up automatically based on the age of the stellar population.

The relation above determines :math:`n_e n_{\mathrm{H}} V`, and from this SLUG computes the nebular emission including the following processes:

* :math:`\mathrm{H}^+` and :math:`\mathrm{He}^+` free-free emission
* :math:`\mathrm{H}` and :math:`\mathrm{He}` bound-free emission
* Hydrogen 2-photon emission
* Hydrogen recombination lines from all lines with upper levels :math:`n_u \leq 25`
* Non-hydrogen line emission based on a tabulation (see below)

Formally, the luminosity per unit wavelength is computed as

.. math:: L_{\lambda,\mathrm{neb}} = \left[\gamma_{\mathrm{ff}}^{(\mathrm{H})} + \gamma_{\mathrm{bf}}^{(\mathrm{H})} + \gamma_{\mathrm{2p}}^{(\mathrm{H})} + \sum_{n,n' \leq 25, n<n'} \alpha_{nn'}^{\mathrm{eff,B,(H)}} E_{nn'}^{(\mathrm{H})} +  x_{\mathrm{He}} \gamma_{\mathrm{ff}}^{(\mathrm{He})} +  x_{\mathrm{He}} \gamma_{\mathrm{bf}}^{(\mathrm{He})} + \sum_i \gamma_{i,\mathrm{line}}^{(\mathrm{M})}\right] n_e n_{\mathrm{H}}{V}

Here :math:`n_e n_{\mathrm{H}} V = \phi_{\mathrm{dust}} Q(\mathrm{H}^0)/ \alpha_{\mathrm{B}}(T)` from photoionization equilibrium, :math:`E_{nn'}` is the energy difference between hydrogen levels :math:`n` and :math:`n'`, and the remaining terms and their sources appearing in this equation are:

* :math:`\gamma_{\mathrm{ff}}^{(\mathrm{H})}` and :math:`\gamma_{\mathrm{ff}}^{(\mathrm{He})}`: HII and HeII free-free emission coefficients; these are computed from eqution 10.1 of `Draine (2011) <http://adsabs.harvard.edu/abs/2011piim.book.....D>`_, using the analytic approximation to the Gaunt factor given by equation 10.8 of the same source 

* :math:`\gamma_{\mathrm{bf}}^{(\mathrm{H})}` and :math:`\gamma_{\mathrm{bf}}^{(\mathrm{He})}`: HI and HeI bound-free emission coefficients; these are computed using the tabulation and interpolation method given in `Ercolano & Storey (2006, MNRAS, 372, 1875) <http://adsabs.harvard.edu/abs/2006MNRAS.372.1875E>`_

* :math:`\alpha_{nn'}^{\mathrm{eff,B,(H)}}` is the effective emission rate coefficient for the :math:`n` to :math:`n'` H recombination line, taken from the tabulation of `Storey & Hummer (1995, MNRAS, 272, 41) <http://adsabs.harvard.edu/abs/1995MNRAS.272...41S>`_

* :math:`\gamma_{i,\mathrm{line}}^{(\mathrm{M})}` is the emissivity for the brightest non-hydrogen lines, computed using a set of pre-tabulated values, following the procedure described in the `SLUG 2 method paper <http://adsabs.harvard.edu/abs/2015MNRAS.452.1447K>`_

* :math:`\gamma_{\mathrm{2p}}^{(\mathrm{H})}`: hydrogen two-photon emissivity, computed as

.. math:: \gamma_{\mathrm{2p}}^{(\mathrm{H})} = \frac{hc}{\lambda^3} I(\mathrm{H}^0) \alpha_{2s}^{\mathrm{eff,(H)}} \frac{1}{1 + \frac{n_{\mathrm{H}} q_{2s-2p,p} + (1+x_{\mathrm{He}}) n_{\mathrm{H}} q_{2s-2p,e}}{A_{2s-1s}}} P_\nu 

Here

  * :math:`I(\mathrm{H}^0)` is the hydrogen ionization potential
  * :math:`\alpha_{2s}^{\mathrm{eff,(H)}}` is the effective recombination rate to the 2s state, taken from the tabulation of `Storey & Hummer (1995, MNRAS, 272, 41) <http://adsabs.harvard.edu/abs/1995MNRAS.272...41S>`_
  * :math:`q_{2s-2p,p}` and :math:`q_{2s-2p,e}` are the collisional rate coefficients for transitions from the 2s to the 2p state induced by collisions with protons and electrons, respectively, taken from `Osterbrock (1989, University Science Books, table 4.10) <http://adsabs.harvard.edu/abs/1989agna.book.....O>`_
  * :math:`A_{2s-1s}` is the Einstein coefficient for the hydrogen 2s-1s two-photon emission process, taken from `Draine (2011, section 14.2.4) <http://adsabs.harvard.edu/abs/2011piim.book.....D>`_
  * :math:`P_\nu` is the frequency distribution for two-photon emission, computed from the analytic approximation of `Nussbaumer & Schmutz (1984, A&A, 138, 495) <http://adsabs.harvard.edu/abs/1984A%26A...138..495N>`_

.. _ssec-extinction:

Extinction
----------

If extinction is enabled, SLUG applies extinction to the stellar spectra and, if nebular processing is enabled as well, to the spectrum that emerges from the nebula. Note that the nebular plus extincted spectrum computation is not fully self-consistent, in that the dust absorption factor :math:`\phi_{\mathrm{dust}}` used in the nebular emission calculation (see :ref:`ssec-nebula`) is not affected by the value of :math:`A_V` used in the calculation.

SLUG computes the extincted spectrum as

.. math:: L_{\lambda,\mathrm{ex}} = L_{\lambda} e^{-\tau_\lambda}

where the optical depth :math:`\tau_\lambda = (\kappa_\lambda / \kappa_V) (A_V/1.086)`, :math:`A_V` is the visual extinction in mag, the factor 1.086 is the conversion between magnitudes and the true dimensionless optical depth, :math:`\kappa_\lambda` is a user-specified input extinction at wavelength :math:`\lambda`, and the V-band mean opacity is defined by

.. math:: \kappa_V = \frac{\int \kappa_\nu R_\nu(V) \, d\nu}{\int R_\nu(V) \, d\nu}

where :math:`R_\nu(V)` is the filter response function as frequency :math:`\nu` for the Johnson V filter. The extinction curve :math:`\kappa_\lambda` can be specified via a user-provided file, or the user may select from a set of pre-defined extinction curves; see :ref:`ssec-ism-keywords` for details.

The computation for :math:`L_{\lambda,\mathrm{neb,ex}}` is analogous.

.. _ssec-yields:

Chemical Yields
---------------

In addition to computing the light output by a stellar population, SLUG can also predict the yield of isotopes. At present SLUG includes yields for core collapse supernovae only. These are based on the yield tables provided by `Sukhbold et al. (2016, ApJ, in press) <http://adsabs.harvard.edu/abs/2015arXiv151004643S/abstract>`_, which provide a finely-spaced set of yields for progenitors of mass 9 - 120 :math:`M_\odot`. Yields from other stellar types will be added in later work.

SLUG operates by separately tracking the isotopes produced by stars winds and by their end-of-life explosions, as a function of progenitor mass. To obtain the yield of any star, SLUG interpolates on the wind and explosive yields using `Steffen interpolation <http://adsabs.harvard.edu/abs/1990A%26A...239..443S>`_ to guarantee that no artifical extrema are produced. The code reports production of 302 distinct isotopes, including many stable species and select unstable ones; for unstable species, radioactive decay over the course of the simulation is properly handled.

At present there is a slight inconsistency in that the evolutionary tracks used by Sukhbold et al. do not precisely match those available for the purposes of omputing spectra and photometr. The death times of stars are computed using the tracks used output the photometric quantities, and thus are slightly inconsistent with the tracks used for the yields. This inconsistency also means that it is not possible to accurately place the wind yields in time; SLUG simply assumes that all yields due to winds are released when the stars end their lives, just as for the explosive yields.
