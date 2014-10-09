.. highlight:: rest

.. _sec-tests:

===============
 Test Problems
===============

This section describes a set of problems that can be used to test and explore the different capabilities of SLUG. SLUG ships a 
set of problems ``problemname`` that are specified by a parameter file ``param/problemname.param``. Problems that require 
multiple simulations are described instead by multiple paramater files, each with unique ID XX:  ``param/problemnameXX.param``. 
Users can reproduce the output of the test problems with the provided executable scripts  ``test/run_problemname.sh``. 
For each problem, a script for analysis is distributed  in ``test/problemname.py``. Details for each test problem are given below.  
Throughout this section, it is assumed that the ``SLUG_DIR`` has been properly set. 

Problem ``example``: basic galaxy simulation
============================================

This problem illustrates the basic usage of \slug\ in ``galaxy`` mode by running 48 realizations of a galaxy with constant 
:math:`\mathrm{SFR}=0.001\; M_\odot\;\mathrm{yr}^{-1}`, up to a maximum time of :math:`2\times 10^8` yr. By issuing the 
command ``test/run_example.sh`` the output files ``SLUG_EXAMPLE*`` are generated. Once the models are ready, 
``python test/plot_example.py`` produces a multi-panel figure ``test/SLUG_EXAMPLE_f1.pdf``. 

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

The second raw shows instead the distribution of the maximum mass of all stars that are still 
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
This new problem highlight instead the additional stochasiticy induced by the second level in the 
hierarchy of ``galaxy`` simulations: how clusters are drawn from the CMF to satisfy the 
target mass. Although may not seem obvious at first, a very important parameter that regulates the 
stochastic behavior of SLUG simulation is the fraction of stars that are formed in clusters, 
:math:`f_c`. In the limit :math:`f_c \rightarrow 0`, SLUG fills a galaxy by drawing stars from the 
IMF. Thus, because the mass of a galaxy is typically much larger than the mass of the upper 
end of the IMF, the effects of mass-constrained sampling highlighted in :ref:`probsampl-label` are simply
not relevant anymore. Conversely, in the limit :math:`f_c \rightarrow 1`, not only the IMF sampling 
plays a role, but clusters themselves contribute to the stochastic behavior. Indeed, similarly 
to the sampling of stars to fill clusters, now clusters are drawn from the CMF to fill the target
galaxy mass. Thus, the same problem of constrained mass sampling applies, inducing additional 
stochasticity.  




This problem highlights the flexible choice of cluster fraction during SLUG simulations.
Introduce importance of this parameter. 

[basic run with different fc]


Problem ``cmfchoice``: different CMF implementations
====================================================

This problem highlights the flexible choice of CMF implementations in SLUG.

[basic run with different cmfs]

Problem ``sfhsampling``: realizations of SFH
============================================

This problem illustrates the conceptual difference between an input SFH and the effective 
realizations produced by SLUG, in comparison to deterministic codes.

[show how a input SFH gets implemented in different realizations]

Problem ``cldisrupt``: cluster disruption at work
=================================================

This problem highlights the flexible choice of CLF implementations in SLUG.

[basic run with different clfs]

Problem ``spectra``: full spectra
=================================

This problem highlights the power of the new feature offered in SLUG v2: the ability to produce 
full spectra and implement dust extinction. 

[basic run with full spectra out: shows stochasticity applied to spectra]

Problem ``redshift``: trivial redshift example
==============================================

This problem shows a trivial example of the redshift capability in SLUG v2.

[basic run with full spectra out at a different redshift]

