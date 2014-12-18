.. highlight:: rest

.. _sec-tests:

===============
 Test Problems
===============

This section describes a set of problems that can be used to test and explore the different capabilities of SLUG. SLUG ships a 
set of problems ``problemname`` that are specified by a parameter file ``param/problemname.param``. Problems that require 
multiple simulations are described instead by multiple paramater files, each with unique ID XX:  ``param/problemnameXX.param``. 
Users can reproduce the output of the test problems with the provided executable scripts  ``test/run_problemname.sh``. 
For each problem, a script for analysis is distributed  in ``test/problemname.py``. Details for each test problem are given below. Throughout this section, it is assumed that the ``SLUG_DIR`` has been properly set.
These test problems are designed to work with outputs in FITS format, but that can be easily changed in the 
``.param`` files. To run all the problems and the analysis scripts in one go, the user can simply 
run ``test/run_alltest.sh``. It will take around 15 minutes 
for the script to complete on a standard laptop. About 700MB of data are generated. 
If SLUG is correctly installed and working, the first part of the script (i.e. the SLUG
simulations) should run flawlessly. The second part of the script relies instead on external python procedures, 
including slugpy, numpy, and matplotlib. While these packages are fairly standard, the user needs to ensure that 
they are properly installed and visible to the script. This script has been written for and tested with Python 2.7.
 

Problem ``example_galaxy``: basic galaxy simulation
===================================================

This problem illustrates the basic usage of \slug\ in ``galaxy`` mode by running 48 realizations of a galaxy with constant 
:math:`\mathrm{SFR}=0.001\; M_\odot\;\mathrm{yr}^{-1}`, up to a maximum time of :math:`2\times 10^8` yr. By issuing the 
command ``test/run_example_galaxy.sh`` the output files ``SLUG_GALAXY_EXAMPLE*`` are generated. Once the models are ready, 
``python test/plot_example_galaxy.py`` produces a multi-panel figure ``test/SLUG_GALAXY_EXAMPLE_f1.pdf``. 

The top-left panel shows the actual mass produced by SLUG for each of the 48 models at different time steps as a 
function of the targeted mass. One can see that SLUG realizations only approximate the desired mass, which is a consequence 
of SLUG core algorithm. The 1:1 relation is shown by a red dashed line. 
The remaining panels show examples of integrated photometry (as labeled) of all simulated galaxies 
at different time steps, as a function of the actual mass. Due to its stochastic nature, SLUG produces 
distributions rather than single values for each time step. The expected rate of ionizing 
photon and the bolometric luminosities for a deterministic model with a
continuous star formation rate of :math:`\mathrm{SFR}=0.001\; M_\odot\;\mathrm{yr}^{-1}` are shown 
by red dashed lines in the relevant panels. 


Problem ``example_cluster``: basic cluster simulation
=====================================================

This problem illustrates the basic usage of SLUG in ``cluster`` mode by running 1000 realizations of a cluster 
with mass 500 :math:`M_\odot`, up to a maximum time of 10 Myr. By issuing the command 
``test/run_example_cluster.sh`` the output files ``SLUG_CLUSTER_EXAMPLE*`` are 
generated. Once the models are ready, ``python test/plot_example_cluster.py`` produces a multi-panel 
figure ``test/SLUG_CLUSTER_EXAMPLE_f1.pdf``. 

This figure is divided in two columns: the left one shows outputs at the first time step, 1 Myr, while 
the second one shows outputs at the last time step, 10 Myr.  The top row shows the actual cluster mass for an 
input mass of :math:`500\;M_\odot`.
In ``cluster`` mode, all clusters are generated at the first time step and they evolve 
passively after that. Thus, the mass does not change. As a consequence of the 
random drawing from the IMF, masses are distributed around the input mass. 
As the wanted mass is large enough to allow for many stars to be drawn, the 
actual mass distribution is narrow. 

The second row shows instead the distribution of the maximum mass of all stars that are still 
alive at a given time step. At 1 Myr, this distribution is a good approximation of the 
input distribution, which is the result of random draws from the IMF. At 10 Myr, which is the 
typical lifetime of a 15-20 :math:`M_\odot` star, the most massive stars have died, and 
SLUG stops following them. The distribution of luminosities, and particularly those 
most sensitive to the presence of massive stars, change accordingly 
(third and fourth row for :math:`Q_{H_0}` and FUV).

.. _probsampl-label:

Problem ``constsampl``: importance of constrained sampling
==========================================================

This problem illustrates in more detail the effects of constrained sampling on SLUG simulations. 
This is the first key ingredient in the core algorithm of SLUG. With the command ``test/run_constsampl.sh``, 
three different ``cluster`` simulations are run, each with 1000 trials, but with masses of :math:`50\;M_\odot`, 
:math:`250\;M_\odot`, and :math:`500\;M_\odot`. A single timestep of :math:`10^6` yr is generated. 
The analysis script ``python test/plot_constsampl.py`` produces a multi-panel 
figure ``test/SLUG_CONSTSAMPL_f1.pdf``. 

This figure shows the maximum mass of the stars in these realizations (top row), the 
rate of ionizing photons :math:`Q_{H_0}` (central row), and the FUV luminosity (bottom row). 
Histograms refer, form left to right, to clusters with :math:`50\;M_\odot`, :math:`250\;M_\odot`, 
and :math:`500\;M_\odot`.

Due to the small timestep, the distributions of stellar masses shown in the top panels reflect 
to good approximation the distribution of the maximum stellar masses that are drawn from the IMF by 
SLUG in each realization. For a cluster of :math:`50\;M_\odot`, the vast majority of the 
stars are drawn below  :math:`20-50\;M_\odot`. This is an obvious consequence of the 
fact that a cluster cannot contain stars much more massive than its own mass. However, stars 
more massive then the targeted mass are not impossible realizations for the default 
sampling algorithm (see below). For instance, if the first star to be drawn has 
mass :math:`60\;M_\odot`, then SLUG would add it to the cluster and stop. Leaving this star out
would indeed be a worse approximation than overshooting the targeted cluster mass by only 
:math:`10\;M_\odot`.  From left to right, one can see that, as the targeted cluster mass increases, the 
histogram shifts to progressively higher masses. In the limit of an infinite cluster, 
all stellar masses would be represented, and the histogram would peak at :math:`120\;M_\odot`.
Essentially, this constrained sampling introduces a stochastic (and not deterministic)
variation in the IMF. An IMF truncated above :math:`60\;M_\odot` would roughly 
approximate the results of the left column; however, a deterministic cut-off 
would not correctly reproduce the non-zero tail at higher masses, thus artificially 
reducing the scatter introduced by random sampling. 

The second and third row simply reflect what said above: for large clusters that can host 
stars at all masses, the luminosity peaks around what is expected according to a deterministic 
stellar population synthesis codes. At lower cluster masses, ionizing and UV fluxes 
are instead suppresses, due to the lack of massive stars. However, tails to high values exist 
in all cases. 
  

Problem ``sampling``: different sampling techniques
===================================================

As highlighted in the previous section, the method with which stars are sampled from the 
IMF has a great influence on the final output. Starting from v2, SLUG has the capability of 
specifying the desired sampling algorithm for a given PDF. 
The command  ``test/run_sampling.sh`` runs four ``cluster`` simulations, each with 1000 trials
of masses of :math:`50\;M_\odot`, and a Kroupa (2002) IMF. 
The following four sampling methods are chosen for each simulation: 1) ``stop_nearest``, 
which is the default in SLUG; 2) ``stop_before``; 3) ``stop_after``; 4) ``sorted_sampling``.
A description of each method is provided in Section :ref:`sampling_metod_label`. 
The analysis script ``python test/plot_sampling.py`` produces a multi-panel 
figure ``test/SLUG_SAMPLING_f1.pdf``. 

By comparing the panels in each column, one can understand the fundamental differences
induced by the sampling technique. The top row shows the maximum stellar mass drawn from the
IMF in each realization. The targeted cluster mass is also shown with red vertical lines.   
In the default mode, SLUG is allowed to overshoot the targeted mass if that constitutes 
a good approximation for the total cluster mass. Thus, a tail at stellar masses above the 
targeted cluster mass is visible. This tail is accentuated when the stop after method 
is selected (third column). In this case, SLUG always overshoots the cluster mass, and thus
extreme realizations above :math:`100\;M_\odot`  are possible. Conversely, in the 
stop after method (second column), SLUG always under-fills the clusters, and (in this case) 
the cluster mass becomes a limit to the maximum stellar mass that can be drawn. A similar effect 
is seen when sorted sampling is enable (fourth column). However, the correspondence between the 
cluster mass and the maximum stellar mass is not trivially established, as it depends on the 
shape of the IMF. The second and third row show how the sampling techniques affect the output 
photometry. 


.. _probimf-label:

Problem ``imfchoice``: different IMF implementations
====================================================

This problem highlights how SLUG can handle different IMF implementations by running 
three simulations with a Kroupa, a Salpeter, and a Chabrier IMF. However, SLUG is not 
restricted to these choices, as the user can in fact easily input an arbitrary IMF. 
The command  ``test/run_imfchoice.sh`` runs three ``cluster`` simulations, each with 1000 trials
of masses of :math:`500\;M_\odot` and different IMF. The analysis script 
``python test/plot_imfchoice.py`` produces a multi-panel figure ``test/SLUG_IMFCHOICE_f1.pdf``. 
Each column shows different statistics for the three IMF. From top to bottom, these are:
the maximum stellar mass in a cluster, the number of stars that SLUG treats stochastically, 
and the distributions of :math:`Q_{H_0}`  and bolometric luminosities. 
As expected for a steep lower-end of the IMF, in the Salpeter case SLUG prefers to fill the 
clusters with a higher number of low mass stars. 


Problem ``clfraction``: cluster fraction at work
================================================

With the exception of the first example, these test problems have focused on how SLUG handles 
cluster simulations, and how these clusters are filled with stars drawn from the IMF. 
This new problem highlights instead the presence of additional stochasticity induced by a 
second level in the hierarchy of ``galaxy`` simulations: how clusters are drawn from the CMF to satisfy the 
targeted galaxy mass. Although it may not appear obvious at first, 
the fraction of stars that are formed in clusters, :math:`f_c`, is a very important parameter that regulates 
the stochastic behavior of SLUG. This can be understood by considering two limiting cases.
In the limit :math:`f_c \rightarrow 0`, SLUG fills a galaxy by drawing stars from the 
IMF. Thus, because the mass of a galaxy is typically much larger than the mass of the upper 
end of the IMF, the effects of mass-constrained sampling highlighted in :ref:`probsampl-label` are simply
not relevant anymore. In this case, stochasticity is minimal.  
Conversely, in the limit :math:`f_c \rightarrow 1`, not only the IMF sampling contributes to the
stochastic behavior of SLUG, but also clusters themselves contribute to additional stochasticity,
as clusters are now drawn from the CMF to fill the targeted galaxy mass following the similar rules 
to those specified for the IMF draws. Thus, in this case, constrained mass sampling applies to both 
stars in clusters and clusters in galaxies, and stochasticity is amplified.  

The command  ``test/run_clfraction.sh`` runs three ``galaxy`` simulations, each with 500 trials
of continuous  SFR :math:`=0.001\rm\;M_\odot\;yr^{-1}` which are evolved for a 
single timestep of  :math:`2\times 10^6\rm\;yr`. A Chabrier IMF and a cluster mass function 
:math:`\propto M^{-2}` are adopted. Cluster disruption is disabled. The three simulations
differ only for the fraction of stars formed in clusters, respectively :math:`f_c=1,0.5,0.01`.
The analysis script ``python test/plot_clfraction.py`` produces a multi-panel figure 
``test/SLUG_CLFRACTION_f1.pdf``. Each column shows properties of simulations for different 
fractions of stars formed in clusters. 

The top row shows the maximum stellar mass in clusters. Clearly, :math:`f_c` has no effect on the way 
clusters are filled up with stars, but the normalization changes. Thus,  the least probable realizations 
in the tail of the distribution simply do not appear for :math:`f_c \rightarrow 0`. The second row 
shows the number of stars in clusters. Obviously, this scales directly with  :math:`f_c`, as it does the number 
of field stars in the third row. This is expected as, by definition, :math:`f_c` regulates the number of stars in 
clusters versus the field. However, as discussed, :math:`f_c` also affects the stochastic behavior of the 
simulation. The fourth row shows histograms of the actual galaxy mass versus the targeted mass (red line).
As :math:`f_c` increases, one can see that the spread around the targeted mass increase. This is again 
a consequence of the mass-constrained sampling and the stop-nearest condition. For :math:`f_c \rightarrow 0`,
the code tries to fill a galaxy of mass :math:`0.001\rm\;M_\odot\;yr^{-1} \times 2\times 10^6\rm\;yr`
with stars. Thus, since the targeted mass is at least a factor of 10 larger than the mass of the 
building block, SLUG can approximate the desired mass very well (to better than :math:`120\rm\;M_\odot`, in fact).
Conversely, for :math:`f_c \rightarrow 1`, SLUG is using clusters as building blocks. As the typical 
mass of the building blocks is now more comparable to the targeted galaxy mass, the problem of the 
mass constrained sampling becomes a relevant one. Not only :math:`f_c` affects the precision with which 
SLUG builds galaxies, but, as shown in the bottom row, it also affects photometry. One can see that 
:math:`Q_{H_0}` increases as :math:`f_c` decreases (the red lines indicate medians). 
The reason for this behavior should now be clear: 
in the case of clustered star formation (:math:`f_c \rightarrow 1`), the mass of the most massive stars 
is subject to the mass constrained sampling of the IMF at the cluster level, reducing the occurrence of 
very massive stars and thus suppressing the flux of ionizing radiation. Conversely, for non clustered star formation 
(:math:`f_c \rightarrow 0`), the sampling of the IMF is constrained only at the galaxy mass level, and since this 
is typically much greater than the mass of the most massive stars, one recovers higher fluxes on average.  


Problem ``cmfchoice``: different CMF implementations
====================================================

Given the ability of SLUG v2 to handle generic PDFs, the user can specify arbitrary CMF, 
similarly to what shown in  :ref:`probimf-label`.
The command  ``test/run_cmfchoice.sh`` runs three ``galaxy`` simulations, each with 500 trials
of continuous  SFR :math:`=0.001\rm\;M_\odot\;yr^{-1}` which are evolved for a 
single timestep of  :math:`2\times 10^6\rm\;yr`. A Chabrier IMF and :math:`f_c=1`
are adopted. Cluster disruption is disabled. The three simulations
differ only for the cluster mass function, which are: 
1) the default powerlaw :math:`M^{-2}` between :math:`20-10^{7}~\rm M_\odot`; 
2) a truncated powerlaw :math:`M^{-2}` between :math:`20-100~\rm M_\odot`;
3) a mass-independent CMF :math:`M^{0}` between :math:`20-10^3~\rm M_\odot`.
The analysis script ``python test/plot_cmfchoice.py`` produces a multi-panel figure 
``test/SLUG_CMFCHOICE_f1.pdf``. Each column shows properties of simulations for the different 
cluster mass functions.

The top row shows the maximum stellar mass in clusters. Compared to the default case, 
the histogram of the truncated CMF is steeper towards low masses. Given that the upper end of the 
CMF is comparable to the maximum stellar mass of the chosen IMF, low stellar masses are typically 
preferred  as a result of the stop-nearest condition. A flat CMF
prefers instead more massive clusters on average, which in turn results in higher probabilities 
of drawing massive stars. In this case, the residual slope of the distribution towards 
low stellar masses is a result of the shape of the IMF. A reflection of the effects induced by the 
shape of the CMF are also apparent in the bottom row, which shows the distribution of 
ionizing photons from these simulations. The second row shows instead the difference 
between the targeted galaxy mass (red line), and the distribution of actual masses.
The spread is minimal for the truncated CMF because, as discussed above, SLUG is using 
small building blocks, and it can approximate the targeted galaxy mass very well. 
Larger spread is visible in the case of the flat CMF, as this choice allows for clusters with masses 
up to :math:`10^3~\rm M_\odot`, without imposing an excess of probability at the low 
mass end. The largest scatter is visible for the default case, as this CMF is virtually 
a pure powerlaw without cutoff at the high mass end, and thus clusters as massive as the entire galaxy 
are accessible to SLUG.  


Problem ``sfhsampling``: realizations of SFH
============================================

The algorithm at the heart of SLUG is quite simple: for a given star formation history 
:math:`\dot\psi(t)` a stellar population with mass :math:`\dot\psi(t)\times \Delta t`
is generated at each timestep, according to the constraints set by IMF, CMF and other 
controlling parameters. As discussed in the previous examples, SLUG builds a best 
approximation for the targeted mass :math:`\dot\psi(t)\times \Delta t`. This means that 
the input SFH and the output SFHs are not identical. SLUG receives an input SFH which 
is used to constrain the rate with which clusters and stars are drawn to achieve the 
desired targeted mass in each timestep. However, the output SFHs are only realizations
and not exact copies  of the input SFH. This problem is designed to illustrate this behavior.  

The command  ``test/run_sfhsampling.sh`` runs two ``galaxy`` simulations, each with 100 trials
of continuous  SFR :math:`=0.0001\rm\;M_\odot\;yr^{-1}` which are evolved for a 
10 timesteps of  :math:`5\times 10^6\rm\;yr`. A Chabrier IMF and a :math:`M^{-2}`
CMF are adopted. Cluster disruption is disabled. The two simulations
differ only for the fraction of stars in clusters, :math:`f_c = 1` and :math:`f_c = 0` respectively. 
The analysis script ``python test/plot_sfhsampling.py`` produces a two-panel figure 
``test/SLUG_SFHSAMPLING_f1.pdf``, showing the box plot for the output SFH of the two simulations
(:math:`f_c = 1` top, and :math:`f_c = 0` bottom).

In each panel, the median SFH over 100 trials is represented by the red lines, while the red squares 
show the mean. The box sizes represent instead the first and third quartile, with the 
ends of the whiskers representing the 5th and 95th percentiles. One can see that the input 
SFH at :math:`\dot\psi(t)=10^{-4}\rm\;M_\odot\;yr^{-1}` is recovered on average, albeit with 
significant variation in each realization. The reason for this variation lies in the fact that, 
at low SFRs, SLUG samples the input SFH with coarse sampling points, which are clusters and stars. 
One can also notice a widely different scatter between the :math:`f_c = 1` and :math:`f_c = 0` 
case. In the former case, the basic elements used by SLUG to sample the targeted mass in  a
given interval are clusters. In the latter case, they are stars. Given that the typical mass of a 
cluster is of the same order of the targeted mass in each interval, the output SFH for 
the :math:`f_c = 1` case are more sensitive to the history of drawings from the CMF. 
Conversely, for  :math:`f_c = 0`, the sampling elements are less massive than the 
targeted mass in a given interval, resulting in an output SFH distribution which is 
better converged towards the input value. Clearly, a comparable amplitude in the scatter 
will be present in the output photometry, especially for the traces that are more sensitive
to variations in the SFHs on short timescales. 


Problem ``cldisrupt``: cluster disruption at work
=================================================

One additional ingredient in SLUG is the lifetime distribution for clusters. Since v2, SLUG is flexible in 
controlling the rate with which clusters are disrupted. This problem shows a comparison between 
two simulations with and without cluster disruption. 

The command  ``test/run_cldisrup.sh`` runs two ``galaxy`` simulations, each with 100 trials
which are evolved in timesteps of  :math:`5\times 10^5\rm\;yr` up to a maximum age of
:math:`1\times 10^7\rm\;yr`. Both simulations are characterized by a burst of star formation 
:math:`=0.001\rm\;M_\odot\;yr^{-1}` within the first Myr. A Chabrier IMF and a :math:`M^{-2}`
CMF are adopted, and :math:`f_c = 1`. For the first simulation, cluster disruption is 
disabled. In the second simulation, cluster disruption operates at times :math:`>1\rm\;Myr`,
with a cluster lifetime function which is a powerlaw of index -1.9. 
The analysis script ``python test/plot_cldisrup.py`` produces the figure ``test/SLUG_CLDISRUP_f1.pdf``.
The two columns show results with (right) and without (left) cluster disruption. 

The first row shows the median stellar mass of the 100 trials as a function of time.
The blue dashed lines show the mass inside the galaxy, while the black solid lines show the 
median mass in clusters. The red band shows the first and fourth quartile of the distribution. 
One can see that in both cases the galaxy mass rises in the first Myr up to the desired 
targeted mass of :math:`=1000\rm\;M_\odot` given the input SFH. After 1Myr, star formation 
stops and the galaxy mass does not evolve with time. Conversely, the cluster mass (black line, red
regions) evolves differently. In the case without cluster disruption, because :math:`f_c = 1`, 
the cluster mass tracks the galaxy mass at all time. When cluster disruption is enabled (right), 
one can see that the mass in clusters rise following the galaxy mass in the first Myr. Past that time, 
clusters start being disrupted and the mass in clusters declines. 
The same behavior is visible in the second row, which shows the median number of alive (black) and 
disrupted (black) clusters. To the left, without cluster disruption, the number of clusters alive
tracks the galaxy mass. Conversely, this distribution declines with time to the right when cluster disruption is
enabled. The complementary quantity (number of disrupted clusters) rises accordingly. 
The last two rows show instead the integrated fluxes in FUV and bolometric luminosity. 
Again, medians are in black and the first and third quartiles in red. One can see a nearly identical distribution 
in the left and right panels. In these simulations, the controlling factors of the integrated photometry 
are the SFH and the sampling techniques, which do not depend on the cluster disruption rate. Clearly, the 
photometry of stars in cluster would exhibit instead a similar dependence to what shown in the top panels. 



Problem ``spectra``: full spectra
=================================

Since v2, SLUG is able to generate spectra for star clusters and for galaxies, which can also be computed for 
arbitrary redshifts. This problem highlights the new features. 
It also demonstrates how SLUG can handle dust extinction, both in a deterministic and stochastic way. 

The command  ``test/run_spectra.sh`` runs four ``galaxy`` simulations, each with 500 trials
of continuous SFR :math:`=0.001\rm\;M_\odot\;yr^{-1}` which are evolved for a 
single timestep of  :math:`2\times 10^6\rm\;yr`. A Chabrier IMF and a :math:`M^{-2}`
CMF are adopted, cluster disruption is disabled, and :math:`f_c = 1`.
The simulations differ in the following way:
1) the reference model, computed without extinction and at :math:`z = 0`;
2) same as the reference model, but at :math:`z = 3`;
3) same as the reference model at :math:`z = 0`, but with a deterministic extinction of :math:`A_V = 0.5` and
a Calzetti+2000 starburst attenuation curve;
4) same as model number 3, but with stochastic extinction.
The analysis script ``python test/plot_spectra.py`` produces the figure ``test/SLUG_SPECTRA_f1.pdf``, 
which shows a gallery of galaxy SEDs for each model. The median SED is shown in black, the blue region 
corresponds to the first and third quartile of the distribution, and the red shaded region 
marks the 5 and 95 percentiles. 

The top panel shows the default model, where stochasticity occurs as detailed in the previous examples. 
The second panel from the top shows instead a model with deterministic extinction. This is simply 
a scaled-down version of the reference model, according to the input dust law and normalization 
coefficient :math:`A_V`. As the dust law extends only to 915 Angstrom the output SED is truncated. 
The third panel shows that, once SLUG handles dust  in a stochastic way, the intrinsic scatter is 
amplified. This is a simple consequence of applying dust extinction with varying normalizations, which 
enhances the final scatter about the median. Finally, the bottom panel shows the trivial case in which 
the spectrum is shifted in wavelength by a constant factor :math:`(1+z)`. Obviously, redshift enhances
the stochasticity in the optical due to a simple shift of wavelengths. 
