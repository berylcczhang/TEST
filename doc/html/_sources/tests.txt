.. highlight:: rest

.. _sec-tests:

Test Problems
=============

This section describes a set of problems that can be used to test and explore the different capabilities of SLUG. SLUG ships a set of problems ``problemname`` that are specified by a parameter file ``param/problemname.param``. Problems that require multiple simulations are described instead by multiple paramater files, each with unique ID XX:  ``param/problemnameXX.param``. Users can reproduce the output of the test problems with the provided executable scripts  ``test/run_problemname.sh``. For each problem, a script for analysis is distributed  in ``test/problemname.py``. Details for each test problem are given below.  Throughout this section, it is assumed that the ``SLUG_DIR`` has been properly set. 

Problem ``example``: basic galaxy simulation
--------------------------------------------

This problem illustrates the basic usage of \slug\ in ``galaxy`` mode by running 48 realizations of a galaxy with constant :math:`\mathrm{SFR}=0.001\; M_\odot\;\mathrm{yr}^{-1}`, up to a maximum time of :math:`2\times 10^8` yr. By issuing the command ``test/run_example.sh`` the output files ``SLUG_EXAMPLE*`` are generated. Once the models are ready, ``python test/plot_example.py`` produces a multi-panel figure ``test/SLUG_EXAMPLE_f1.pdf``. 

The top-left panel shows the actual mass produced by SLUG for each of the 48 models at different time steps as a function of the targeted mass. One can see that SLUG realizations only approximate the desired mass, which is a consequence 
of SLUG core algorithm. The 1:1 relation is shown by a red dashed line. 
The remaining panels show examples of integrated photometry (as labeled) of all simulated galaxies 
at different time steps, as a function of the actual mass. Due to its stochastic nature, SLUG produces 
distributions rather than single values for each time step. The expected rate of ionizing 
photon and the bolometric luminosities for a deterministic model with a
continuous star formation rate of :math:`\mathrm{SFR}=0.001\; M_\odot\;\mathrm{yr}^{-1}` are shown 
by red dashed lines in the relevant panels. 


Problem ``example_cluster``: basic cluster simulation
-----------------------------------------------------

This problem illustrates the basic usage of SLUG in ``cluster`` mode by running 1000 realizations of a cluster with mass 500 :math:`M_\odot`,
up to a maximum time of 10 Myr. By issuing the command 
``test/run_example_cluster.sh`` the output files ``SLUG_CLUSTER_EXAMPLE*`` are 
generated. Once the models are ready, ``python test/plot_example_cluster.py`` produces a multi-panel figure ``test/SLUG_CLUSTER_EXAMPLE_f1.pdf``. 

This figure is divided in two columns: the left one shows outputs at the first time step, 1 Myr, while the second one shows outputs at the last time step, 10 Myr.  The top row shows the actual cluster mass for an input mass of :math:`500\;M_\odot`.
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


Problem ``constsampl``: importance of constrained sampling
-------------------------------------------------------------

This problem illustrates in more details the effects of constrained sampling on SLUG simulations, 
which is the first key ingredient in the core algorithm of SLUG.

 


[change cluster mass to show how max mass changes]


Problem ``sampling``: different sampling techniques
-----------------------------------------------------

This problem highlights the flexible choice of sampling techniques in SLUG, which is 
a new capability of v2.

[basic run with different sampling]


Problem ``imfchoice``: different IMF implementations
------------------------------------------------------

This problem highlights the flexible choice of IMF implementations in SLUG.

[basic run with different imfs]


Problem ``cmfchoice``: different CMF implementations
------------------------------------------------------

This problem highlights the flexible choice of CMF implementations in SLUG.

[basic run with different cmfs]

Problem ``sfhsampling``: realizations of SFH
----------------------------------------------

This problem illustrates the conceptual difference between an input SFH and the effective 
realizations produced by SLUG, in comparison to deterministic codes.

[show how a input SFH gets implemented in different realizations]

Problem ``cldisrupt``: cluster disruption at work
---------------------------------------------------

This problem highlights the flexible choice of CLF implementations in SLUG.

[basic run with different clfs]

Problem ``clfraction``: cluster fraction at work
--------------------------------------------------

This problem highlights the flexible choice of cluster fraction during SLUG simulations.

[basic run with different fc]

Problem ``spectra``: full spectra
-----------------------------------

This problem highlights the power of the new feature offered in SLUG v2: the ability to produce 
full spectra. 

[basic run with full spectra out: shows stochasticity applied to spectra]

Problem ``redshift``: trivial redshift example
------------------------------------------------

This problem shows a trivial example of the redshift capability in SLUG v2.

[basic run with full spectra out at a different redshift]


