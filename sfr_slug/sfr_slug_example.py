# Import the sfr_slug tool and other libraries
from slugpy.sfr_slug import *
import numpy as np
import matplotlib.pyplot as plt
import time

# Create a star formation rate estimator based on ionizing photon
# fluxes, using a flat prior on log SFR
sfrest=sfr_slug(filters='QH0', priors='flat')

# Print out when we're done reading, just to make the point about how
# quick computation is
print("Done reading library...")
tottime = 0.0

# Create array of input point mass-estimated SFRs from 10^-6 to 10^-1
# Msun/yr, spaced logarithmically in steps of 1 dex. These are the
# SFRs that one would infer from a non-stochastic estimate.
logSFR_phot = np.linspace(-6, -1, 7)

# Compute posterior probability distributions for true SFR, assuming
# negligible photometric errors
lasttime = time.clock()
logSFR_true, pdf_flat = sfrest.mpdf(logSFR_phot)
newtime = time.clock()
tottime = newtime - lasttime

# Same calculation assuming 0.5 dex errors
lasttime = time.clock()
logSFR_true, pdf_flat_err = sfrest.mpdf(logSFR_phot, 0.5)
newtime = time.clock()
tottime = tottime + newtime - lasttime

# Change the prior to Schechter
sfrest.priors = 'schechter'

# Now repeat the two calculations for Schechter function prior
lasttime = time.clock()
logSFR_true, pdf_sch = sfrest.mpdf(logSFR_phot)
logSFR_true, pdf_sch_err = sfrest.mpdf(logSFR_phot, 0.5)
newtime = time.clock()
tottime = tottime + newtime - lasttime

# Print out that we're done computing
endtime = time.clock()
print("PDF computation done in {:f} sec, average {:f} sec per input".
      format(tottime, tottime / (4*len(logSFR_phot))))

# Make figure
plt.figure(1, figsize=(7,7))

# Make plot for flat prior
plt.subplot(211)
p1,=plt.plot(logSFR_true-logSFR_phot[1], pdf_flat[1,:], 'b', lw=2)
p1_err,=plt.plot(logSFR_true-logSFR_phot[1], pdf_flat_err[1,:], 'b--', lw=2)
p2,=plt.plot(logSFR_true-logSFR_phot[3], pdf_flat[3,:], 'g', lw=2)
p2_err,=plt.plot(logSFR_true-logSFR_phot[3], pdf_flat_err[3,:], 'g--', lw=2)
p3,=plt.plot(logSFR_true-logSFR_phot[5], pdf_flat[5,:], 'r', lw=2)
p3_err,=plt.plot(logSFR_true-logSFR_phot[5], pdf_flat_err[5,:], 'r--', lw=2)
leg1=plt.legend([p1,p2,p3],
                [r"$\log\,\mathrm{SFR}_{Q(\mathrm{H}^0)} = -5$",
                 r"$\log\,\mathrm{SFR}_{Q(\mathrm{H}^0)} = -3$",
                 r"$\log\,\mathrm{SFR}_{Q(\mathrm{H}^0)} = -1$"],
                title=r"Flat, $\sigma = 0$", loc='upper left')
leg2=plt.legend([p1_err,p2_err,p3_err],
                [r"$\log\,\mathrm{SFR}_{Q(\mathrm{H}^0)} = -5$",
                 r"$\log\,\mathrm{SFR}_{Q(\mathrm{H}^0)} = -3$",
                 r"$\log\,\mathrm{SFR}_{Q(\mathrm{H}^0)} = -1$"],
                title=r"Flat, $\sigma = 0.5$ dex", loc='upper right')
plt.gca().add_artist(leg1)
plt.xlim([-6,6])
#plt.xlabel(r"$\log\,\mathrm{SFR}-\log\,\mathrm{SFR}_{Q(\mathrm{H}^0)}$")
plt.setp(plt.gca(), 'xticklabels', [])
plt.ylabel("PDF")

# Make plot for Schechter prior
plt.subplot(212)
p1,=plt.plot(logSFR_true-logSFR_phot[1], pdf_sch[1,:], 'b', lw=2)
p1_err,=plt.plot(logSFR_true-logSFR_phot[1], pdf_sch_err[1,:], 'b--', lw=2)
p2,=plt.plot(logSFR_true-logSFR_phot[3], pdf_sch[3,:], 'g', lw=2)
p2_err,=plt.plot(logSFR_true-logSFR_phot[3], pdf_sch_err[3,:], 'g--', lw=2)
p3,=plt.plot(logSFR_true-logSFR_phot[5], pdf_sch[5,:], 'r', lw=2)
p3_err,=plt.plot(logSFR_true-logSFR_phot[5], pdf_sch_err[5,:], 'r--', lw=2)
leg1=plt.legend([p1,p2,p3],
                [r"$\log\,\mathrm{SFR}_{Q(\mathrm{H}^0)} = -5$",
                 r"$\log\,\mathrm{SFR}_{Q(\mathrm{H}^0)} = -3$",
                 r"$\log\,\mathrm{SFR}_{Q(\mathrm{H}^0)} = -1$"],
                title=r"Shechter, $\sigma = 0$", loc='upper left')
leg2=plt.legend([p1_err,p2_err,p3_err],
                [r"$\log\,\mathrm{SFR}_{Q(\mathrm{H}^0)} = -5$",
                 r"$\log\,\mathrm{SFR}_{Q(\mathrm{H}^0)} = -3$",
                 r"$\log\,\mathrm{SFR}_{Q(\mathrm{H}^0)} = -1$"],
                title=r"Schechter, $\sigma = 0.5$ dex", loc='upper right')
plt.gca().add_artist(leg1)
plt.xlim([-6,6])
plt.xlabel(r"$\log\,\mathrm{SFR}-\log\,\mathrm{SFR}_{Q(\mathrm{H}^0)}$")
plt.ylabel("PDF")
plt.ylim([0, 1.78])

plt.subplots_adjust(hspace=0, top=0.95)
plt.savefig('sfr_slug.pdf')
