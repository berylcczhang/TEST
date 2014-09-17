#
#This script displays the output of the test problem: example_cloudy
#See section 10 of the slug manual 
#

import numpy as np
import matplotlib.pyplot as plt
import os
import os.path as osp
import copy
from os import listdir 

try:
    from slugpy import *    # If slugpy is already in our path
except ImportError:
    # If import failed, try to find slugpy in $SLUG_DIR
    if 'SLUG_DIR' in os.environ:
        import sys
        cur_path = copy.deepcopy(sys.path)
        sys.path.append(os.environ['SLUG_DIR'])
        from slugpy import *
        sys.path = cur_path
    else:
        raise ImportError("No module named slugpy")

# Get slug dir
if 'SLUG_DIR' in os.environ:
    slugdir = os.environ['SLUG_DIR']
else:
    slugdir = '.'

#set model name and output type
modname='SLUG_CLOUDY'

#load output 
data=read_integrated(modname)

#get the mean spectra
mean_emit = np.mean(data.cloudy_emit, axis=2)
mean_trans_emit = np.mean(data.cloudy_trans_emit, axis=2)

#get the 10th, 50th, and 90th percentiles
perc_emit = np.array(np.percentile(data.cloudy_emit, (10,50,90), axis=2))
perc_trans_emit = np.array(np.percentile(data.cloudy_trans_emit, 
                                         (10,50,90), axis=2))

#plot the spectra
plt.figure(1)
wl = data.cloudy_wl
pt,=plt.loglog(wl, mean_emit[:,0], 'b', lw=2)
pte,=plt.loglog(wl, mean_trans_emit[:,0], 'r', lw=2)
pt1,=plt.loglog(wl, perc_emit[1,:,0], 'b--')
pte1,=plt.loglog(wl, perc_trans_emit[1,:,0], 'r--')
plt.fill_between(wl, perc_emit[0,:,0], perc_emit[2,:,0], color='b',
                 alpha=0.25)
plt.fill_between(wl, perc_trans_emit[0,:,0], perc_trans_emit[2,:,0], 
                 color='r', alpha=0.25)
plt.xlim([5e2, 1e6])
plt.ylim([1e32, 1e37])
plt.xlabel(r'$\lambda$ [$\AA$]')
plt.ylabel(r'$F_\lambda$ [erg s$^{-1}$ $\AA^{-1}$]')
plt.legend([pt, pt1, pte, pte1],
           ["Nebular, mean",
            "Nebular, median",
            "Nebular+stellar, mean",
            "Nebular+stellar, median"],
           loc='upper right')

# Save
plt.savefig('test/'+modname+'_f1.pdf')

#extract Ha, Hbeta, NII, OIII line luminosities
Ha_idx = np.where(data.cloudy_linewl == 6563)
Hb_idx = np.where(data.cloudy_linewl == 4861)
NII_idx = np.where(data.cloudy_linewl == 6584)
OIII_idx = np.where(data.cloudy_linewl == 4959)
lHa = data.cloudy_linelum[Ha_idx,0,:].flatten()
lHb = data.cloudy_linelum[Hb_idx,0,:].flatten()
lNII = data.cloudy_linelum[NII_idx,0,:].flatten()
lOIII = data.cloudy_linelum[OIII_idx,0,:].flatten()

#plot the BPT diagram
plt.figure(2)
plt.scatter(np.log10(lNII/lHa), np.log10(lOIII/lHb))
plt.xlim([-2, 0])
plt.ylim([-2, 1])
plt.xlabel(r'log NII/H$\alpha$')
plt.ylabel(r'log OIII/H$\beta$')

# Save
plt.savefig('test/'+modname+'_fw.pdf')


